"""Updates of MSpektra from its GitHub releases.

Each release carries the full zip and a signed list of the parts of the program (update.json and
update.json.sig: the version, a hash of every part, the notes). The files of a part are stored one after the
other in the full zip, so a part is one byte range of it (offset, size, SHA-256 in the list). At start
MSpektra reads the list of the latest release; when it is newer, the update bar appears. Update downloads
only the byte ranges of the parts that differ from those installed (install.json in the program folder),
checks the signature of the list and the SHA-256 of every range, and unpacks them into _update\\staging. "Restart now" starts update_apply.py (from a copy in the temp folder, with the Python of
the folder), which swaps the parts after MSpektra has closed, starts it again and puts the previous version
back when the new one does not start. The previous version stays in _update\\old: About > "Back to ..."
returns to it.

The parts are the files of the folder (MSpektra.exe, README.txt, LICENSES ...), the files of _portable and
each folder of _portable (msengine, baf2sql, timsdata, wheels ...). The Python of the folder is not a part:
when a release changes it ("runtime" in the list), the update needs the full zip and the download page
opens instead. Nothing is sent but the requests for these files (no user data, no statistics)."""
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from decimal import Decimal, InvalidOperation

REPO = "its-kietsu/MSpektra"
RELEASES = "https://github.com/%s/releases" % REPO
# the public key of the MSpektra release key (the private key stays on the release computer)
PUBLIC_KEY = bytes.fromhex("a35709dee950d523bf6bd69d27071d592625621f6a8611fc40e8915f32397bd1")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
UPD = os.path.join(ROOT, "_update")
STAGING = os.path.join(UPD, "staging")
OLD = os.path.join(UPD, "old")
PENDING = os.path.join(UPD, "pending.json")
_STATE = {"started": False, "result": None, "listeners": []}


def _log(*a):
    try:
        print("[update]", *a)
    except Exception:
        pass


def version_key(v):
    """Versions compare as decimal numbers, as MSpektra counts them: 4.3 is newer than 4.24."""
    try:
        return Decimal(str(v).strip().lstrip("vV"))
    except (InvalidOperation, ValueError):
        return Decimal(0)


def current_version():
    try:
        import unidec_theme
        return unidec_theme.APP_VERSION
    except Exception:
        return "0"


def _read_json(path):
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


def installed():
    """The list of the parts installed (install.json, written with the release or by the last update)."""
    return _read_json(os.path.join(ROOT, "install.json"))


def _settings():
    try:
        import unidec_theme
        return unidec_theme._load() or {}
    except Exception:
        return {}


def enabled():
    """Updates only for a folder made from a release (install.json) that MSpektra may write to, and not when
    switched off ("check_updates": false in the settings or MSPEKTRA_NO_UPDATE set)."""
    if os.environ.get("MSPEKTRA_NO_UPDATE") or _settings().get("check_updates", True) is False:
        return False
    if installed() is None:
        return False
    try:
        os.makedirs(UPD, exist_ok=True)
        probe = os.path.join(UPD, ".write")
        with open(probe, "w") as f:
            f.write("1")
        os.remove(probe)
        return True
    except OSError:
        return False


def _url(name, version=None):
    if version is None:
        return "%s/latest/download/%s" % (RELEASES, name)
    return "%s/download/v%s/%s" % (RELEASES, version, name)


def _open(name, version=None, timeout=20):
    """A file of the latest release (or of version): from GitHub, or from the folder MSPEKTRA_UPDATE_FEED
    (tests: a folder with the files of a release)."""
    feed = os.environ.get("MSPEKTRA_UPDATE_FEED")
    if feed:
        return open(os.path.join(feed, name), "rb")
    import urllib.request
    req = urllib.request.Request(_url(name, version), headers={"User-Agent": "MSpektra/%s" % current_version()})
    return urllib.request.urlopen(req, timeout=timeout)


def _open_range(name, version, offset, size, timeout=60):
    """size bytes from offset of a file of the release (a part inside the full zip)."""
    feed = os.environ.get("MSPEKTRA_UPDATE_FEED")
    if feed:
        f = open(os.path.join(feed, name), "rb")
        f.seek(offset)
        return _Limited(f, size)
    import urllib.request
    req = urllib.request.Request(_url(name, version), headers={
        "User-Agent": "MSpektra/%s" % current_version(), "Range": "bytes=%d-%d" % (offset, offset + size - 1)})
    resp = urllib.request.urlopen(req, timeout=timeout)
    if resp.status != 206:  # (the whole file: no)
        resp.close()
        raise ValueError("the server did not send a part of %s" % name)
    return _Limited(resp, size)


class _Limited:
    def __init__(self, f, size):
        self.f, self.left = f, size

    def read(self, n):
        block = self.f.read(min(n, self.left)) if self.left > 0 else b""
        self.left -= len(block)
        return block

    def __enter__(self):
        return self

    def __exit__(self, *a):
        self.f.close()


def _unpack(path, staging):
    """The zip entries of a byte range of the full zip (local headers, one after the other) into staging,
    without the top folder ("MSpektra <version>/")."""
    import struct
    import zlib
    with open(path, "rb") as f:
        while True:
            head = f.read(30)
            if len(head) < 30:
                if head:
                    raise ValueError("a part ends in the middle of a file")
                return
            sig, _, flag, method, _, _, crc, csize, usize, nlen, xlen = struct.unpack("<4s5H3L2H", head)
            if sig != b"PK\x03\x04" or flag & 0x08:
                raise ValueError("a part is not a run of zip entries")
            raw = f.read(nlen)
            extra = f.read(xlen)
            name = raw.decode("utf-8" if flag & 0x800 else "cp437")
            if csize == 0xFFFFFFFF or usize == 0xFFFFFFFF:  # zip64: the sizes in the extra field
                i = 0
                while i + 4 <= len(extra):
                    tag, ln = struct.unpack("<2H", extra[i:i + 4])
                    if tag == 1:
                        vals = list(struct.unpack("<%dQ" % (ln // 8), extra[i + 4:i + 4 + ln]))
                        if usize == 0xFFFFFFFF:
                            usize = vals.pop(0)
                        if csize == 0xFFFFFFFF:
                            csize = vals.pop(0)
                    i += 4 + ln
            rel = name.split("/", 1)[1] if "/" in name else name
            dest = _inside(staging, rel)
            if name.endswith("/"):
                os.makedirs(dest, exist_ok=True)
                continue
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            dec = zlib.decompressobj(-15) if method == 8 else None
            if method not in (0, 8):
                raise ValueError("%s: unknown compression" % rel)
            left, c = csize, 0
            with open(dest, "wb") as out:
                while left > 0:
                    block = f.read(min(left, 1 << 20))
                    if not block:
                        raise ValueError("a part ends in the middle of %s" % rel)
                    left -= len(block)
                    if dec:
                        block = dec.decompress(block)
                    c = zlib.crc32(block, c)
                    out.write(block)
                if dec:
                    block = dec.flush()
                    c = zlib.crc32(block, c)
                    out.write(block)
            if c & 0xFFFFFFFF != crc:
                raise ValueError("%s is damaged (CRC differs)" % rel)


def check():
    """The newer release, if there is one: its list of parts (verified) with what this folder needs of it;
    None when there is none, or nothing could be read (offline: nothing is said)."""
    inst = installed()
    if inst is None:
        return None
    try:
        with _open("update.json") as f:
            data = f.read()
        with _open("update.json.sig") as f:
            sig = bytes.fromhex(f.read().decode("ascii").strip())
    except Exception as ex:
        _log("no update information (%s)" % ex)
        return None
    import ed25519
    if not ed25519.verify(PUBLIC_KEY, data, sig):
        _log("update.json is not signed with the MSpektra release key: ignored")
        return None
    try:
        man = json.loads(data.decode("utf-8"))
    except ValueError:
        return None
    if man.get("app") != "MSpektra" or version_key(man.get("version")) <= version_key(current_version()):
        return None
    skip = _read_json(os.path.join(UPD, "skip.json")) or {}
    if skip.get("version") == man.get("version"):
        _log("%s was put back after it did not start: not offered again" % man.get("version"))
        return None
    have = inst.get("items", {})
    man["_changed"] = [k for k, v in man.get("items", {}).items() if have.get(k, {}).get("hash") != v.get("hash")]
    man["_removed"] = [k for k in have if k not in man.get("items", {})]
    man["_bytes"] = sum(int(man["items"][k].get("size", 0)) for k in man["_changed"])
    # the Python of the folder ("runtime" until 4.25, "python" since 4.26: "runtime" of the newer lists sends
    # 4.25, which cannot read parts in the full zip, to the download page)
    man["_full"] = (man.get("python", man.get("runtime")) != inst.get("python", inst.get("runtime"))
                    or not man.get("zip"))
    man["_raw"] = data.decode("utf-8")
    return man


def staged():
    """An update downloaded before and not installed yet ("Restart to update")."""
    ready = _read_json(os.path.join(STAGING, "ready.json"))
    if ready and version_key(ready.get("version")) > version_key(current_version()):
        return ready
    return None


def _inside(base, rel):
    p = os.path.normpath(os.path.join(base, rel))
    if not (p + os.sep).startswith(os.path.normpath(base) + os.sep):
        raise ValueError("a path outside the folder: %s" % rel)
    return p


def download(man, progress=None, cancelled=lambda: False):
    """Downloads the parts this folder needs into _update\\staging; progress(done_bytes, total_bytes).
    True when everything was downloaded and checked; raises on a bad download."""
    if os.path.isdir(STAGING):
        shutil.rmtree(STAGING, ignore_errors=True)
    dl = os.path.join(UPD, "download")
    os.makedirs(dl, exist_ok=True)
    os.makedirs(STAGING, exist_ok=True)
    total, done = max(1, man["_bytes"]), 0
    for key in man["_changed"]:
        it = man["items"][key]
        target = os.path.join(dl, "part.bin")
        h = hashlib.sha256()
        with _open_range(man["zip"], man["version"], int(it["offset"]), int(it["size"])) as src, \
                open(target, "wb") as out:
            while True:
                if cancelled():
                    return False
                block = src.read(1 << 16)
                if not block:
                    break
                out.write(block)
                h.update(block)
                done += len(block)
                if progress:
                    progress(done, total)
        if h.hexdigest() != it["sha256"]:
            raise ValueError("%s is damaged (checksum differs)" % key)
        _unpack(target, STAGING)
        os.remove(target)
    with open(os.path.join(STAGING, "install.json"), "w", encoding="utf-8") as f:
        f.write(man["_raw"])
    with open(os.path.join(STAGING, "ready.json"), "w", encoding="utf-8") as f:
        json.dump({"version": man["version"], "from": current_version(), "changed": man["_changed"],
                   "removed": man["_removed"], "notes": man.get("notes", [])}, f, indent=1)
    _log("%s downloaded (%d parts, %.1f MB)" % (man["version"], len(man["_changed"]), done / 1048576.0))
    return True


def _helper(args):
    """Starts update_apply.py from a copy in the temp folder (the files of _portable are replaced) with the
    Python of this folder, detached from MSpektra."""
    src = os.path.join(ROOT, "_portable", "update_apply.py")
    dst = os.path.join(tempfile.gettempdir(), "mspektra_update_%d.py" % os.getpid())
    shutil.copyfile(src, dst)
    pyw = os.path.join(ROOT, "python", "pythonw.exe")
    if not os.path.isfile(pyw):
        pyw = sys.executable
    flags = 0x00000008 | 0x00000200  # DETACHED_PROCESS | CREATE_NEW_PROCESS_GROUP
    subprocess.Popen([pyw, "-I", dst, "--root", ROOT, "--pid", str(os.getpid())] + list(args),
                     cwd=tempfile.gettempdir(), creationflags=flags, close_fds=True)
    _log("installer started: %s" % " ".join(args or ["install"]))


def apply_and_restart():
    """The downloaded update is installed after MSpektra has closed, then MSpektra starts again."""
    _helper([])
    _close_all()


def previous():
    """The version kept in _update\\old (About > Back to ...), or None."""
    info = _read_json(os.path.join(OLD, "info.json"))
    if info and os.path.isfile(os.path.join(OLD, "journal.json")):
        return info.get("from")
    return None


def back_to_previous():
    _helper(["--rollback"])
    _close_all()


def _close_all():
    import wx
    for w in list(wx.GetTopLevelWindows()):
        try:
            w.Close()
        except Exception:
            pass


def mark_started():
    """The first window of a new version is up: the installer does not put the old version back."""
    try:
        if os.path.isfile(PENDING):
            os.remove(PENDING)
            _log("version %s started after the update" % current_version())
    except OSError:
        pass


# ------------------------------------------------------------------ interface
def on_update(listener):
    """listener(info): called (on the main thread) when an update is found or downloaded; info has
    "version" and "staged" (True: downloaded, restart to install). Called at once when it is known."""
    _STATE["listeners"].append(listener)
    if _STATE["result"] is not None:
        _call(listener, _STATE["result"])


def _call(listener, info):
    try:
        listener(info)
    except RuntimeError:  # (its window was closed)
        pass
    except Exception as ex:
        _log("update button: %s" % ex)


def _found(info):
    _STATE["result"] = info
    for li in list(_STATE["listeners"]):
        _call(li, info)


def bar(parent):
    """The information bar of an update (start screen, under the top bar of LCMS and HRMS Analysis): hidden
    until a newer version is found; Update opens the download window, Restart now installs a downloaded one,
    x hides the bar until MSpektra starts again. Returns the bar (a wx.Panel in the parent's sizer)."""
    b = _bar_class()(parent)
    on_update(b.set_info)
    return b


def _bar_class():
    import wx
    import unidec_theme as T
    from unidec_theme import C, ui_font

    class InfoBar(wx.Panel):
        """Light blue bar with a rounded border (Windows 11 information bar): (i), Update available, the
        version, its button and x."""

        def __init__(self, parent):
            wx.Panel.__init__(self, parent, style=wx.BORDER_NONE)
            self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
            self.info, self.closed = None, False
            self.btn = T._cls("FlatButton")(self, label="Update", kind="primary", height=26, padx=12)
            self.btn.Bind(wx.EVT_BUTTON, self.on_button)
            self.SetMinSize(wx.Size(-1, self.FromDIP(38)))
            self.Bind(wx.EVT_PAINT, self._paint)
            self.Bind(wx.EVT_SIZE, lambda e: (self._place(), self.Refresh(), e.Skip()))
            self.Bind(wx.EVT_LEFT_UP, self._click)
            self.Bind(wx.EVT_MOTION, self._motion)
            self.Hide()

        def _texts(self):
            if not self.info:
                return "", ""
            if self.info.get("staged"):
                return "Update ready", "MSpektra %s is ready to install." % self.info["version"]
            return "Update available", "MSpektra %s is ready to download." % self.info["version"]

        def _x_rect(self):
            w, h = self.GetClientSize()
            s = self.FromDIP(22)
            return wx.Rect(w - s - self.FromDIP(8), (h - s) // 2, s, s)

        def _place(self):
            """The button right after the text."""
            dc = wx.ClientDC(self)
            head, text = self._texts()
            dc.SetFont(ui_font(9.5, 600))
            w1 = dc.GetTextExtent(head)[0]
            dc.SetFont(ui_font(9.5, 400))
            w2 = dc.GetTextExtent(text)[0]
            x = self.FromDIP(12 + 16 + 10) + w1 + self.FromDIP(10) + w2 + self.FromDIP(12)
            bw, bh = self.btn.GetSize()
            self.btn.SetPosition(wx.Point(x, (self.GetClientSize()[1] - bh) // 2))

        def set_info(self, info):
            self.info = info
            self.btn.SetLabel("Restart now" if info.get("staged") else "Update")
            if not self.closed:
                self.Show()
                self._relayout()
            self._place()
            self.Refresh()

        def _relayout(self):
            top = self.GetTopLevelParent()
            grow = getattr(top, "info_bar_grows", False)  # (the start screen gets taller instead of its tiles smaller)
            if grow and self.IsShown() and not getattr(self, "_grown", False):
                self._grown = True
                w, h = top.GetSize()
                top.SetSize(wx.Size(w, h + self.FromDIP(38 + 8)))
            self.GetParent().Layout()
            top.Layout()

        def on_button(self, e=None):
            if not self.info:
                return
            if self.info.get("staged"):
                wx.CallAfter(apply_and_restart)
            else:
                wx.CallAfter(show_dialog, self.GetTopLevelParent(), _STATE["result"] or self.info)

        def _motion(self, e):
            self.SetCursor(wx.Cursor(wx.CURSOR_HAND if self._x_rect().Contains(e.GetPosition()) else wx.CURSOR_ARROW))

        def _click(self, e):
            if self._x_rect().Contains(e.GetPosition()):
                self.closed = True
                self.Hide()
                self.GetParent().Layout()

        def _paint(self, e):
            dc = wx.AutoBufferedPaintDC(self)
            dc.SetBackground(wx.Brush(self.GetParent().GetBackgroundColour()))
            dc.Clear()
            gc = T.crisp(dc)
            w, h = self.GetClientSize()
            k = self.FromDIP(10) / 10.0
            gc.SetBrush(wx.Brush(wx.Colour("#EAF1FB")))
            gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(wx.Colour("#C7D9F4")).Width(1)))
            gc.DrawRoundedRectangle(0.5, 0.5, w - 1, h - 1, 4 * k)
            cx, cy, r = self.FromDIP(12) + 8 * k, h / 2.0, 8 * k  # (i)
            gc.SetPen(wx.TRANSPARENT_PEN)
            gc.SetBrush(wx.Brush(wx.Colour(C["accent"])))
            gc.DrawEllipse(cx - r, cy - r, 2 * r, 2 * r)
            gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(wx.Colour("#FFFFFF")).Width(1.6 * k).Cap(wx.CAP_ROUND)))
            gc.StrokeLine(cx, cy - 1 * k, cx, cy + 3 * k)
            gc.StrokeLine(cx, cy - 3.6 * k, cx, cy - 3.4 * k)
            head, text = self._texts()
            x = self.FromDIP(12 + 16 + 10)
            gc.SetFont(ui_font(9.5, 600), wx.Colour(C["text"]))
            tw, th = gc.GetTextExtent(head)
            gc.DrawText(head, x, (h - th) / 2.0)
            gc.SetFont(ui_font(9.5, 400), wx.Colour("#444B55"))
            gc.DrawText(text, x + tw + self.FromDIP(10), (h - th) / 2.0)
            xr = self._x_rect()  # x
            gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(wx.Colour("#5B6475")).Width(1.4 * k).Cap(wx.CAP_ROUND)))
            q = 4.5 * k
            mx, my = xr.x + xr.width / 2.0, xr.y + xr.height / 2.0
            gc.StrokeLine(mx - q, my - q, mx + q, my + q)
            gc.StrokeLine(mx - q, my + q, mx + q, my - q)

    return InfoBar


def start(delay_ms=5000):
    """Looks for an update once per process, a few seconds after the first window (in the background)."""
    if _STATE["started"]:
        return
    _STATE["started"] = True
    mark_started()
    if not enabled():
        return
    import wx
    st = staged()
    if st:
        wx.CallAfter(_found, {"version": st["version"], "staged": True, "notes": st.get("notes", [])})
        return

    def work():
        man = check()
        if man is not None:
            wx.CallAfter(_found, {"version": man["version"], "staged": False, "man": man,
                                  "notes": man.get("notes", [])})

    wx.CallLater(delay_ms, lambda: threading.Thread(target=work, name="update check", daemon=True).start())


def show_dialog(parent, info):
    import wx
    import unidec_theme as T
    dlg = UpdateDialog(parent, info)
    dlg.ShowModal()
    dlg.Destroy()


def _dialog_class():
    import wx
    import unidec_theme as T
    from unidec_theme import C, ui_font

    class UpdateDialog(wx.Dialog):
        """Version, what is new, a progress bar; Update downloads, Restart now installs."""

        def __init__(self, parent, info):
            wx.Dialog.__init__(self, parent, title="Update")
            self.SetBackgroundColour(wx.Colour(C["panel"]))
            self.info = info
            self.man = info.get("man")
            self._cancel = False
            self._busy = False
            pad = self.FromDIP(16)
            vs = wx.BoxSizer(wx.VERTICAL)
            head = wx.StaticText(self, label="MSpektra %s" % info["version"])
            head.SetFont(ui_font(13, 600))
            vs.Add(head, 0, wx.LEFT | wx.RIGHT | wx.TOP, pad)
            notes = [n for n in info.get("notes", []) if n][:8]
            if notes:
                txt = wx.StaticText(self, label="\n".join("• " + n for n in notes))
                txt.SetFont(ui_font(9.5))
                txt.Wrap(self.FromDIP(440))
                vs.Add(txt, 0, wx.LEFT | wx.RIGHT | wx.TOP, pad)
            self.gauge = wx.Gauge(self, range=1000, size=self.FromDIP(wx.Size(440, 6)))
            vs.Add(self.gauge, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, pad)
            self.status = wx.StaticText(self, label=self._ready_text())
            self.status.SetForegroundColour(wx.Colour(C["muted"]))
            self.status.SetFont(ui_font(9))
            vs.Add(self.status, 0, wx.LEFT | wx.RIGHT | wx.TOP, self.FromDIP(6))
            flat = T._cls("FlatButton")
            self.later = flat(self, label="Later", kind="secondary", height=32)
            self.go = flat(self, label="", kind="primary", height=32, min_width=110)
            self._set_go()
            self.later.Bind(wx.EVT_BUTTON, self.on_later)
            self.go.Bind(wx.EVT_BUTTON, self.on_go)
            bs = wx.BoxSizer(wx.HORIZONTAL)
            bs.AddStretchSpacer(1)
            bs.Add(self.later, 0, wx.RIGHT, self.FromDIP(8))
            bs.Add(self.go, 0)
            vs.Add(bs, 0, wx.EXPAND | wx.ALL, pad)
            self.SetSizerAndFit(vs)
            self.CentreOnParent()
            self.Bind(wx.EVT_CLOSE, self.on_later)

        def _ready_text(self):
            if self.info.get("staged"):
                return "Downloaded"
            if self.man and self.man.get("_full"):
                return "This version needs the full download"
            if self.man:
                return "%.1f MB" % (self.man["_bytes"] / 1048576.0)
            return ""

        def _set_go(self):
            if self.info.get("staged"):
                self.go.SetLabel("Restart now")
                self.gauge.SetValue(1000)
            elif self.man and self.man.get("_full"):
                self.go.SetLabel("Download page")
            else:
                self.go.SetLabel("Update")
            self.Layout()

        def on_later(self, e=None):
            if self._busy:
                self._cancel = True
                return
            self.EndModal(wx.ID_CANCEL)

        def on_go(self, e=None):
            if self._busy:
                return
            if self.info.get("staged"):
                self.EndModal(wx.ID_OK)
                wx.CallAfter(apply_and_restart)
                return
            if self.man and self.man.get("_full"):
                import webbrowser
                webbrowser.open("%s/tag/v%s" % (RELEASES, self.info["version"]))
                self.EndModal(wx.ID_OK)
                return
            self._busy = True
            self.go.Enable(False)
            self.later.SetLabel("Cancel")
            self.status.SetLabel("Downloading")

            def progress(done, total):
                wx.CallAfter(self._progress, done, total)

            def work():
                try:
                    ok = download(self.man, progress, lambda: self._cancel)
                    err = None
                except Exception as ex:
                    ok, err = False, str(ex)
                wx.CallAfter(self._done, ok, err)

            threading.Thread(target=work, name="update download", daemon=True).start()

        def _progress(self, done, total):
            try:
                self.gauge.SetValue(int(1000 * done / max(total, 1)))
                self.status.SetLabel("%.1f of %.1f MB" % (done / 1048576.0, total / 1048576.0))
            except RuntimeError:
                pass

        def _done(self, ok, err):
            try:
                self._busy = False
                self.later.SetLabel("Later")
                self.go.Enable(True)
                if ok:
                    self.info = dict(self.info, staged=True)
                    _found(self.info)
                    self.status.SetLabel("Downloaded")
                    self._set_go()
                elif self._cancel:
                    self.EndModal(wx.ID_CANCEL)
                else:
                    self.gauge.SetValue(0)
                    self.status.SetLabel("Download failed: %s" % (err or "try again later"))
            except RuntimeError:
                pass

    return UpdateDialog


def UpdateDialog(parent, info):
    return _dialog_class()(parent, info)

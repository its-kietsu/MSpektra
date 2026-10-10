"""Updates of MSpektra from its GitHub releases.

Each release carries, next to the full zip, a signed list of the parts of the program (update.json and
update.json.sig: the version, a hash of every part, the notes) and one small zip per part
(update-<part>-<hash>.zip). At start MSpektra reads the list of the latest release; when it is newer, the
"Update" button appears. A click downloads only the parts that differ from those installed (install.json
in the program folder), checks the signature of the list and the SHA-256 of every zip, and unpacks them into
_update\\staging. "Restart now" starts update_apply.py (from a copy in the temp folder, with the Python of
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
import zipfile
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
_STATE = {"started": False, "result": None, "listeners": [], "popped": False}


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
    man["_full"] = man.get("runtime") != inst.get("runtime")
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
        target = os.path.join(dl, it["asset"])
        h = hashlib.sha256()
        with _open(it["asset"], man["version"], timeout=60) as src, open(target, "wb") as out:
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
            raise ValueError("%s is damaged (checksum differs)" % it["asset"])
        with zipfile.ZipFile(target) as z:
            for info in z.infolist():
                dest = _inside(STAGING, info.filename)
                if info.is_dir():
                    os.makedirs(dest, exist_ok=True)
                    continue
                os.makedirs(os.path.dirname(dest), exist_ok=True)
                with z.open(info) as zs, open(dest, "wb") as zd:
                    shutil.copyfileobj(zs, zd)
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
    if not _STATE["popped"]:  # once per start: a notice in the corner of the window, without a click
        _STATE["popped"] = True
        try:
            _notice(info)
        except Exception as ex:
            _log("update notice: %s" % ex)


def _notice(info):
    """A small notice at the top right of the MSpektra window: the new version, Update and Later. It stays
    until one of them is clicked (or its window closes)."""
    import wx
    import unidec_theme as T
    from unidec_theme import C, ui_font
    tops = [w for w in wx.GetTopLevelWindows() if w.IsShown() and not isinstance(w, wx.Dialog)]
    if not tops:
        return
    act = wx.GetActiveWindow()
    parent = act.GetTopLevelParent() if act is not None and act.GetTopLevelParent() in tops else tops[0]
    f = wx.Frame(parent, style=wx.FRAME_TOOL_WINDOW | wx.FRAME_FLOAT_ON_PARENT | wx.FRAME_NO_TASKBAR
                 | wx.BORDER_SIMPLE)
    p = wx.Panel(f)
    p.SetBackgroundColour(wx.Colour(C["panel"]))
    vs = wx.BoxSizer(wx.VERTICAL)
    head = wx.StaticText(p, label="MSpektra %s %s" % (info["version"], "is ready to install" if info.get("staged")
                                                       else "is available"))
    head.SetFont(ui_font(10.5, 600))
    vs.Add(head, 0, wx.LEFT | wx.RIGHT | wx.TOP, p.FromDIP(14))
    notes = [n for n in info.get("notes", []) if n][:2]
    if notes:
        t = wx.StaticText(p, label="\n".join(notes))
        t.SetFont(ui_font(9))
        t.SetForegroundColour(wx.Colour(C["muted"]))
        t.Wrap(p.FromDIP(300))
        vs.Add(t, 0, wx.LEFT | wx.RIGHT | wx.TOP, p.FromDIP(14))
    flat = T._cls("FlatButton")
    later = flat(p, label="Later", kind="secondary", height=30)
    go = flat(p, label="Restart now" if info.get("staged") else "Update", kind="primary", height=30, min_width=90)
    bs = wx.BoxSizer(wx.HORIZONTAL)
    bs.AddStretchSpacer(1)
    bs.Add(later, 0, wx.RIGHT, p.FromDIP(8))
    bs.Add(go, 0)
    vs.Add(bs, 0, wx.EXPAND | wx.ALL, p.FromDIP(14))
    p.SetSizer(vs)
    vs.Fit(f)
    f.SetClientSize(p.GetBestSize())

    def place():
        try:
            r = parent.GetScreenRect()
            w, h = f.GetSize()
            f.SetPosition(wx.Point(r.x + r.width - w - parent.FromDIP(24), r.y + parent.FromDIP(110)))
        except RuntimeError:
            pass

    def update(e):
        f.Destroy()
        if info.get("staged"):
            wx.CallAfter(apply_and_restart)
        else:
            wx.CallAfter(show_dialog, parent, _STATE["result"] or info)

    later.Bind(wx.EVT_BUTTON, lambda e: f.Destroy())
    go.Bind(wx.EVT_BUTTON, update)
    try:
        parent.Bind(wx.EVT_MOVE, lambda e: (place(), e.Skip()))
        parent.Bind(wx.EVT_SIZE, lambda e: (place(), e.Skip()))
    except Exception:
        pass
    place()
    f.Show()


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

"""
Modern interface add-on for the portable UniDec package.

* Sharp on scaled displays: UniDec tells Windows it handles display scaling
  itself (per-monitor DPI awareness), so text and plots are no longer
  stretched and blurred at 125 % / 150 % scaling. Fixed control sizes from
  UniDec's code are scaled to match.
* Main UniDec window: labelled toolbar (Open, Process, Deconvolve, Pick peaks,
  Run all, Replot), clean section headers with a thin colour mark, blue main
  buttons, renamed quick buttons (Processing / UniDec / Peaks), plots on white
  cards with a title (the titles are not part of the saved figures).
* Peak list: white rows with a small coloured marker instead of fully
  coloured rows.
* Windows reopen at their last size and position; file dialogs start in the
  last used data folder.
* Launcher: light frosted glass tiles with an icon for each tool.

Every button keeps its function and keyboard shortcut (the new buttons use the
same IDs as the original ones). UniDec's own code is not modified.
"""
import json
import os

_ST = {"settings_path": None, "system_scale": 1.0, "dpi_aware": False, "installed": False}

# light, cool palette (white surfaces on a pale blue-grey, deep blue accent)
C = {
    "bg": "#F1F4F9", "panel": "#FFFFFF", "line": "#E5EAF2", "line2": "#D2DAE6",
    "text": "#131A24", "muted": "#57616E", "faint": "#67707C",
    "accent": "#2A62C4", "accent_hover": "#3470D4", "accent_down": "#16407F", "accent_end": "#16407F",
    "accent_text": "#17418F", "accent_bg": "#E6EDFA",
    "hover": "#EFF2F7", "down": "#E6EDFA", "sel": "#E6EDFA", "tint": "#F5F8FC",
}
GROUP = {"blue": "#2F6FCB", "yellow": "#C99A1E", "red": "#D2456F", "green": "#0BA064"}

# Plus Jakarta Sans (SIL OFL 1.1, bundled in _portable\fonts) for headings,
# toolbar, launcher; native form fields keep the Windows system font
# Regular and Bold share one GDI family; Medium/SemiBold map onto them (their
# separate families rendered word spaces too narrow in testing)
_FACES = {400: "Plus Jakarta Sans", 500: "Plus Jakarta Sans", 600: "Plus Jakarta Sans",
          700: "Plus Jakarta Sans", 800: "Plus Jakarta Sans ExtraBold"}


def register_fonts(folder):
    """Make the bundled fonts available to this process only (nothing is
    installed on the computer)."""
    if os.name != "nt" or not folder or not os.path.isdir(folder):
        return 0
    import ctypes
    try:
        import wx
    except Exception:
        wx = None
    n = 0
    for name in sorted(os.listdir(folder)):
        if name.lower().endswith(".ttf"):
            path = os.path.join(folder, name)
            gdi_ok = False
            try:
                gdi_ok = ctypes.windll.gdi32.AddFontResourceExW(path, 0x10, 0) > 0
            except Exception:
                pass
            # GDI+ can reject a private font even when GDI accepts it. wx logs
            # that failure as a modal error; an optional heading font must
            # fall back to Segoe UI without interrupting the user.
            private_ok = False
            if wx is not None:
                quiet = wx.LogNull()
                try:
                    private_ok = bool(wx.Font.AddPrivateFont(path))
                except Exception:
                    pass
                finally:
                    del quiet  # restore logging immediately, including real app errors
            if gdi_ok and private_ok:
                n += 1
    _ST.pop("jakarta_ok", None)
    _ST["fonts_registered"] = n
    return n


def _jakarta_ok():
    if "jakarta_ok" not in _ST:
        try:
            import wx
            _ST["jakarta_ok"] = bool(_ST.get("fonts_registered")) and \
                wx.FontEnumerator.IsValidFacename("Plus Jakarta Sans")
        except Exception:
            _ST["jakarta_ok"] = False
    return _ST["jakarta_ok"]


def ui_font(points, weight=400):
    """Headings (13 pt and larger): Plus Jakarta Sans if available. Smaller
    text: Segoe UI, the Windows interface font, which is hinted for ClearType
    and stays sharp at small sizes (Plus Jakarta Sans looked soft there)."""
    import wx
    if points >= 13 and _jakarta_ok():
        f = wx.Font(wx.FontInfo(points).FaceName(_FACES.get(weight, _FACES[400])))
        if weight in (600, 700):
            f.SetWeight(wx.FONTWEIGHT_BOLD)
        return f
    f = wx.Font(wx.FontInfo(points).FaceName("Segoe UI"))
    f.SetWeight({500: wx.FONTWEIGHT_MEDIUM, 600: wx.FONTWEIGHT_SEMIBOLD, 700: wx.FONTWEIGHT_BOLD,
                 800: wx.FONTWEIGHT_EXTRABOLD}.get(weight, wx.FONTWEIGHT_NORMAL))
    return f

PLOT_TITLES = [("plot1", "MS data and fit"), ("plot2", "Mass distribution"), ("plot3", "m/z vs charge"),
               ("plot4", "Individual peaks"), ("plot5", "Mass vs charge"), ("plot6", "Peak intensities")]

ICONS = {
    "open": '<path d="M3 7a2 2 0 0 1 2-2h4l2 2h8a2 2 0 0 1 2 2v8a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2z"/>',
    "process": '<path d="M4 6h16M7 12h10M10 18h4"/>',
    "deconvolve": '<path d="M3 18c2 0 2-9 4-9s2 6 4 6 2-11 4-11 2 14 4 14h2"/>',
    "pick": '<path d="M3 20h18M7 20v-9M12 20V6M17 20v-6M9.5 3l2.5 2.5L14.5 3"/>',
    "run": '<path d="M8 5v14l11-7z" fill="{c}"/>',
    "replot": '<path d="M20 11a8 8 0 1 0-2.3 5.7M20 5v6h-6"/>',
}
ICONS.setdefault("report", '<path d="M6 3h8l4 4v14H6z"/><path d="M14 3v4h4"/><path d="M9 12h6M9 15.5h6M9 19h4"/>')


def _log(*args):
    try:
        print("[interface]", *args)
    except Exception:
        pass


# --------------------------------------------------------------------------
# settings (shared file with the plot style add-on)
# --------------------------------------------------------------------------
def read_settings(p, tries=20):
    """The settings file p as a dict ({} when there is none or it is empty). Read again a few times when
    another process holds it right now (Windows: it is being swapped in); raises when it still cannot be read
    or holds no valid settings."""
    import time
    bad = 0
    for attempt in range(tries):
        try:
            if not os.path.isfile(p):
                return {}
            with open(p, "r", encoding="utf-8") as fh:
                text = fh.read()
            if not text.strip():
                return {}  # (nothing to lose)
            d = json.loads(text)
            if not isinstance(d, dict):
                raise ValueError("the settings file holds no settings")
            return d
        except ValueError:  # (a damaged file stays damaged: read again twice only)
            bad += 1
            if bad >= 3 or attempt == tries - 1:
                raise
            time.sleep(0.01)
        except OSError:
            if attempt == tries - 1:
                raise
            time.sleep(0.01 + 0.005 * attempt)


_LOCK = {}


class _settings_lock(object):
    """Only one process at a time reads, changes and writes the settings file (a lock on the file p + ".lock",
    held while the file is read, changed and swapped in). After timeout s the save goes on without it."""

    def __init__(self, p, timeout=10.0):
        self.p, self.timeout, self.fd, self.mode = p, timeout, None, None

    def __enter__(self):
        import threading
        import time
        lk = _LOCK.setdefault("thread", threading.Lock())
        lk.acquire()
        try:
            self.fd = os.open(self.p + ".lock", os.O_RDWR | os.O_CREAT, 0o666)
        except OSError:
            self.fd = None
            return self
        t_end = time.time() + self.timeout
        while True:
            try:
                try:
                    import msvcrt
                    os.lseek(self.fd, 0, 0)
                    msvcrt.locking(self.fd, msvcrt.LK_NBLCK, 1)
                    self.mode = "msvcrt"
                except ImportError:
                    import fcntl
                    fcntl.flock(self.fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    self.mode = "fcntl"
                return self
            except OSError:
                if time.time() > t_end:
                    _log("settings lock not taken in %.0f s: saving without it" % self.timeout)
                    return self
                time.sleep(0.002)

    def __exit__(self, *exc):
        try:
            if self.fd is not None:
                try:
                    if self.mode == "msvcrt":
                        import msvcrt
                        os.lseek(self.fd, 0, 0)
                        msvcrt.locking(self.fd, msvcrt.LK_UNLCK, 1)
                    elif self.mode == "fcntl":
                        import fcntl
                        fcntl.flock(self.fd, fcntl.LOCK_UN)
                except OSError:
                    pass
                os.close(self.fd)
        finally:
            self.fd = None
            _LOCK["thread"].release()
        return False


def update_settings(p, update):
    """update (a dict) written into the settings file p, keeping every other setting. The file is read, changed
    and swapped in under a lock shared by every process of the program; when it exists but cannot be read, it
    is kept unchanged and this save is dropped (never a file with the other settings lost). True when saved."""
    if not p:
        return False
    # written whole, then swapped in: a crash cannot leave a truncated file. The
    # temporary name is per process: the start screen, the analysis windows and
    # the Deconvolute window (own processes) save to the same file
    tmp = "%s.%d.tmp" % (p, os.getpid())
    try:
        os.makedirs(os.path.dirname(p) or ".", exist_ok=True)
        with _settings_lock(p):
            data = read_settings(p)  # raises: the save is dropped, the file kept
            data.update(update)
            with open(tmp, "w", encoding="utf-8") as fh:
                json.dump(data, fh, indent=1)
            for attempt in range(60):
                try:
                    os.replace(tmp, p)
                    break
                except PermissionError:  # Windows: another process is reading the file right now
                    if attempt == 59:
                        raise
                    import time
                    time.sleep(0.005 + 0.002 * attempt)
        return True
    except Exception as e:
        _log("settings not saved:", e)
        try:
            os.remove(tmp)
        except OSError:
            pass
        return False


def _load():
    p = _ST["settings_path"]
    if not p:
        return {}
    try:
        return read_settings(p)
    except Exception as e:
        if _ST.get("load_error") != str(e):  # (once, not at every read)
            _ST["load_error"] = str(e)
            _log("settings not read:", e)
        return {}


def _save(update):
    p = _ST["settings_path"]
    if not p:
        return
    _ST["gen"] = _ST.get("gen", 0) + 1  # (readers that cache the file: two saves within one file time tick)
    revs = _ST.setdefault("revs", {})  # saves of each setting by this process (undo.py: changed here or elsewhere)
    for k in update:
        revs[k] = revs.get(k, 0) + 1
    update_settings(p, update)


# --------------------------------------------------------------------------
# 1. DPI awareness (call before any window exists)
# --------------------------------------------------------------------------
def enable_dpi_awareness():
    """Per-monitor DPI awareness; returns the system scale factor."""
    if os.name != "nt":
        return 1.0
    import ctypes
    aware = False
    try:
        aware = bool(ctypes.windll.user32.SetProcessDpiAwarenessContext(ctypes.c_void_p(-4)))  # PMv2
    except Exception:
        pass
    if not aware:
        try:
            aware = ctypes.windll.shcore.SetProcessDpiAwareness(2) == 0
        except Exception:
            pass
    if not aware:
        try:
            aware = bool(ctypes.windll.user32.SetProcessDPIAware())
        except Exception:
            pass
    scale = 1.0
    try:
        scale = ctypes.windll.user32.GetDpiForSystem() / 96.0
    except Exception:
        try:
            hdc = ctypes.windll.user32.GetDC(0)
            scale = ctypes.windll.gdi32.GetDeviceCaps(hdc, 88) / 96.0
            ctypes.windll.user32.ReleaseDC(0, hdc)
        except Exception:
            pass
    _ST["dpi_aware"] = aware
    _ST["system_scale"] = scale if scale > 0 else 1.0
    return _ST["system_scale"]


def _patch_display_size():
    """UniDec picks plot sizes from the display size in 96-dpi pixels; keep
    that behaviour when the process is DPI aware."""
    import wx
    s = _ST["system_scale"]
    if not _ST["dpi_aware"] or s <= 1.01:
        return
    orig = wx.GetDisplaySize

    def GetDisplaySize():
        sz = orig()
        return wx.Size(int(round(sz[0] / s)), int(round(sz[1] / s)))

    wx.GetDisplaySize = GetDisplaySize


def _scale(win):
    try:
        return float(win.GetDPIScaleFactor()) if _ST["dpi_aware"] else 1.0
    except Exception:
        return _ST["system_scale"] if _ST["dpi_aware"] else 1.0


def _skip_scaling(w):
    name = w.__class__.__name__
    if name in ("CaptionBar", "FigureCanvasWxAgg", "FigureCanvasWx", "VisibleCanvas", "HtmlWindow", "FlatButton",
                "ModernToolbar", "Chip", "GlassLauncher"):
        return True
    try:
        from unidec.modules.plotting.PlotBase import PlotBase
        if isinstance(w, PlotBase):
            return True
    except Exception:
        pass
    return False


def scale_fixed_sizes(root):
    """Scale explicit (96-dpi) minimum sizes set by UniDec's code."""
    import wx
    f = _scale(root)
    if f <= 1.01 or getattr(root, "_udp_scaled", False):
        return
    root._udp_scaled = True

    def walk(w):
        for c in w.GetChildren():
            if isinstance(c, wx.TopLevelWindow) or _skip_scaling(c):
                continue
            try:
                ms = c.GetMinSize()
                if ms.width > 0 or ms.height > 0:
                    c.SetMinSize(wx.Size(int(ms.width * f) if ms.width > 0 else -1,
                                         int(ms.height * f) if ms.height > 0 else -1))
            except Exception:
                pass
            walk(c)

    walk(root)


def _from_unidec(win):
    """A window of UniDec's own code (its sizes are in 96-dpi pixels); the
    windows of MSpektra give theirs with FromDIP."""
    mod = type(win).__module__ or ""
    return mod == "unidec" or mod.startswith("unidec.")


# --------------------------------------------------------------------------
# 1b. windows inside the screen
# --------------------------------------------------------------------------
def display_area(win=None, own_first=False):
    """Work area (screen pixels, without the task bar) of the display a top
    level window belongs to: the display of its parent window, else the one
    it is on, else the one under the mouse, else the main display.
    own_first: the display it is on first (a window the user may have moved)."""
    import wx

    def index(w):
        try:
            return wx.Display.GetFromWindow(w)
        except Exception:
            return wx.NOT_FOUND

    n = wx.NOT_FOUND
    if win is not None:
        if own_first:
            n = index(win)
        par = win.GetParent()
        top = wx.GetTopLevelParent(par) if par is not None else None
        if n == wx.NOT_FOUND and top is not None:
            try:
                if top.IsShown() and not top.IsIconized():
                    n = index(top)
            except RuntimeError:
                pass
        if n == wx.NOT_FOUND and par is None:
            n = index(win)
    if n == wx.NOT_FOUND:
        try:
            n = wx.Display.GetFromPoint(wx.GetMousePosition())
        except Exception:
            n = wx.NOT_FOUND
    if n == wx.NOT_FOUND and win is not None:
        n = index(win)
    try:
        return wx.Display(n if n != wx.NOT_FOUND and 0 <= n < wx.Display.GetCount() else 0).GetClientArea()
    except Exception:
        return wx.Rect(0, 0, 1280, 720)


def fit_rect(rect, area, best=None, gap=0):
    """Pure geometry of fit_to_screen: (x, y, w, h) of a window with the
    rect `rect` (x, y, w, h) inside the work area `area`: at least `best`
    (its contents) when that fits, at most the area less `gap` on each side,
    moved only as far as needed to lie inside the area. A window made
    smaller or larger keeps its centre."""
    x, y, w, h = [int(v) for v in rect]
    ax, ay, aw, ah = [int(v) for v in area]
    nw, nh = w, h
    if best is not None:
        nw, nh = max(nw, int(best[0])), max(nh, int(best[1]))
    nw = max(1, min(nw, aw - 2 * gap))
    nh = max(1, min(nh, ah - 2 * gap))
    x += (w - nw) // 2
    y += (h - nh) // 2
    x = max(ax, min(x, ax + aw - nw))
    y = max(ay, min(y, ay + ah - nh))
    return x, y, nw, nh


def fit_to_screen(win, grow=True, area=None):
    """Keeps a top level window (frame or dialog) inside the work area of its
    display: never larger than it, as large as its contents when they fit
    (grow), moved onto it as far as needed. Maximised, minimised and full
    screen windows are left alone. Applied to every window when it is shown
    (_install_generic_show); windows that change their size while open call
    it again. Returns the work area used (None: nothing done)."""
    import wx
    try:
        if win.IsMaximized() or win.IsIconized() or win.IsFullScreen():
            return None
    except Exception:
        return None
    if area is None:
        area = display_area(win, own_first=win.IsShown() or getattr(win, "_udp_fitted", False))
    gap = win.FromDIP(4)
    room = wx.Size(max(1, area.width - 2 * gap), max(1, area.height - 2 * gap))
    # a minimum size larger than the screen would keep the window too large
    mn = win.GetMinSize()
    if mn.width > room.width or mn.height > room.height:
        win.SetMinSize(wx.Size(min(mn.width, room.width) if mn.width > 0 else -1,
                               min(mn.height, room.height) if mn.height > 0 else -1))
    best = None
    if grow and win.GetSizer() is not None:
        b = win.GetBestSize()
        best = (b.width, b.height)
    r = win.GetRect()
    new = fit_rect((r.x, r.y, r.width, r.height), (area.x, area.y, area.width, area.height), best, gap)
    if new != (r.x, r.y, r.width, r.height):
        win.SetSize(*new)
    return area


# --------------------------------------------------------------------------
# small custom controls
# --------------------------------------------------------------------------
def _bitmap_from_svg(body, size, colour):
    import wx
    try:
        from wx.svg import SVGimage
        svg = ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none" stroke="%s" '
               'stroke-width="1.8" stroke-linecap="round" stroke-linejoin="round">%s</svg>'
               % (colour, body.replace("{c}", colour)))
        img = SVGimage.CreateFromBytes(svg.encode("utf-8"))
        return img.ConvertToScaledBitmap(wx.Size(size, size))
    except Exception as e:
        _log("icon:", e)
        return None


def _controls():
    import wx

    class FlatButton(wx.Control):
        """Owner-drawn button; posts a normal EVT_BUTTON with its own ID."""

        def __init__(self, parent, id=wx.ID_ANY, label="", icon=None, kind="ghost", tooltip=None,
                     height=30, padx=12, min_width=0):
            wx.Control.__init__(self, parent, id, style=wx.BORDER_NONE)
            self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
            self._label, self._icon, self._kind = label, icon, kind
            self._h, self._padx, self._minw = height, padx, min_width
            self._hover = self._down = False
            self._bmps = {}
            for ev, fn in ((wx.EVT_PAINT, self._on_paint), (wx.EVT_ENTER_WINDOW, self._on_enter),
                           (wx.EVT_LEAVE_WINDOW, self._on_leave), (wx.EVT_LEFT_DOWN, self._on_down),
                           (wx.EVT_LEFT_DCLICK, self._on_down), (wx.EVT_LEFT_UP, self._on_up),
                           (wx.EVT_MOUSE_CAPTURE_LOST, self._on_lost), (wx.EVT_SIZE, lambda e: self.Refresh())):
                self.Bind(ev, fn)
            if tooltip:
                self.SetToolTip(tooltip)
            pts = 9 if height < 28 else 9.5
            self.SetFont(ui_font(pts, 700 if kind == "primary" else 600 if icon else 500))
            self._fit()

        def _fit(self):
            dc = wx.ClientDC(self)
            dc.SetFont(self.GetFont())
            tw, th = dc.GetTextExtent(self._label or "Xg")
            w = 2 * self.FromDIP(self._padx) + (tw if self._label else 0)
            if self._icon:
                w += self.FromDIP(16) + (self.FromDIP(7) if self._label else 0)
            w = max(w, self.FromDIP(self._minw))
            hh = self.FromDIP(int(round(self._h)))
            self.SetMinSize(wx.Size(w, hh))
            self.SetInitialSize(wx.Size(w, hh))

        def SetLabel(self, label):
            self._label = label
            self._fit()
            self.Refresh()

        def GetLabel(self):
            return self._label

        def Enable(self, enable=True):
            r = wx.Control.Enable(self, enable)
            self.Refresh()
            return r

        def _on_enter(self, e):
            self._hover = True
            self.Refresh()

        def _on_leave(self, e):
            self._hover = False
            self.Refresh()

        def _on_down(self, e):
            if not self.IsEnabled():
                return
            self._down = True
            if not self.HasCapture():
                self.CaptureMouse()
            self.Refresh()

        def _on_up(self, e):
            was = self._down
            self._down = False
            if self.HasCapture():
                self.ReleaseMouse()
            self.Refresh()
            if was and self.IsEnabled() and self.GetClientRect().Contains(e.GetPosition()):
                evt = wx.CommandEvent(wx.EVT_BUTTON.typeId, self.GetId())
                evt.SetEventObject(self)
                wx.CallAfter(self._fire, evt)

        def _fire(self, evt):
            try:
                self.GetEventHandler().ProcessEvent(evt)
            except RuntimeError:
                pass  # window closed meanwhile

        def _on_lost(self, e):
            self._down = False
            self.Refresh()

        def _icon_bmp(self, colour):
            if not self._icon:
                return None
            key = colour
            if key not in self._bmps:
                body = ICONS.get(self._icon)
                if body is None:  # unknown icon: draw the label only (never break the paint)
                    self._bmps[key] = None
                    return None
                try:
                    self._bmps[key] = _bitmap_from_svg(body, self.FromDIP(16), colour)
                except Exception:
                    self._bmps[key] = None
            return self._bmps[key]

        def _on_paint(self, e):
            dc = wx.AutoBufferedPaintDC(self)
            bg = self.GetParent().GetBackgroundColour()
            dc.SetBackground(wx.Brush(bg))
            dc.Clear()
            gc = crisp(dc)
            w, h = self.GetClientSize()
            k, en = self._kind, self.IsEnabled()
            big = h >= self.FromDIP(30)
            # pill shape; primary buttons leave room for a soft shadow below
            inset = self.FromDIP(2) if (k == "primary" and big) else 0
            x0, y0, bw, bh = 1.0, 0.5 + (inset * 0.3), w - 2.0, h - 1.0 - inset
            r = bh / 2.0 if big or k == "ghost" else self.FromDIP(7)
            fg = C["text"] if en else "#9AA3AD"
            if k == "primary":
                if not en:
                    c1, c2 = "#9DB6E0", "#8AA3CF"
                elif self._down:
                    c1, c2 = "#1E4F9E", "#123565"
                elif self._hover:
                    c1, c2 = C["accent_hover"], "#1B4A92"
                else:
                    c1, c2 = C["accent"], C["accent_end"]
                if big and en:
                    gc.SetPen(wx.TRANSPARENT_PEN)
                    for i, a in ((3, 22), (1.5, 34)):
                        gc.SetBrush(wx.Brush(wx.Colour(23, 65, 143, a)))
                        gc.DrawRoundedRectangle(x0 + i * 0.5, y0 + i, bw - i, bh, r)
                gc.SetBrush(gc.CreateLinearGradientBrush(x0, y0, x0 + bw * 0.35, y0 + bh, wx.Colour(c1),
                                                         wx.Colour(c2)))
                gc.SetPen(wx.TRANSPARENT_PEN)
                gc.DrawRoundedRectangle(x0, y0, bw, bh, r)
                # top sheen
                gc.SetBrush(gc.CreateLinearGradientBrush(x0, y0, x0, y0 + bh * 0.55, wx.Colour(255, 255, 255, 46),
                                                         wx.Colour(255, 255, 255, 0)))
                gc.DrawRoundedRectangle(x0 + 1, y0 + 1, bw - 2, bh * 0.55, max(1, r - 1))
                fg = "#FFFFFF"
            elif k == "secondary":
                fill = C["down"] if self._down else C["hover"] if self._hover else C["panel"]
                gc.SetBrush(wx.Brush(wx.Colour(fill)))
                gc.SetPen(wx.Pen(wx.Colour(C["line2"]), 1))
                gc.DrawRoundedRectangle(x0, y0, bw, bh, r)
            else:  # ghost
                fill = C["down"] if self._down else C["hover"] if self._hover else None
                if fill:
                    gc.SetBrush(wx.Brush(wx.Colour(fill)))
                    gc.SetPen(wx.TRANSPARENT_PEN)
                    gc.DrawRoundedRectangle(x0, y0, bw, bh, r)
                if self._down or self._hover:
                    fg = C["accent_text"] if en else fg
            gc.SetFont(self.GetFont(), wx.Colour(fg))
            tw, th = gc.GetTextExtent(self._label) if self._label else (0, 0)
            icol = "#FFFFFF" if k == "primary" else (C["accent_text"] if (self._hover or self._down) else
                                                     "#3B4450") if en else "#9AA3AD"
            bmp = self._icon_bmp(icol)
            iw = self.FromDIP(16) if bmp else 0
            gap = self.FromDIP(7) if (bmp and self._label) else 0
            x = x0 + (bw - (iw + gap + tw)) / 2.0
            cy = y0 + bh / 2.0
            if bmp:
                gc.DrawBitmap(bmp, x, cy - iw / 2.0, iw, iw)
            if self._label:
                gc.DrawText(self._label, x + iw + gap, cy - th / 2.0)

    class Chip(wx.Control):
        def __init__(self, parent, label):
            wx.Control.__init__(self, parent, style=wx.BORDER_NONE)
            self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
            self._label, self._value = label, ""
            self._value_tooltip = ""
            self.SetFont(ui_font(9, 500))
            self._bold = ui_font(9, 700)
            self.Bind(wx.EVT_PAINT, self._on_paint)
            self._fit()
            self.Hide()

        def _fit(self):
            dc = wx.ClientDC(self)
            f = self.GetFont()
            dc.SetFont(f)
            w1, _ = dc.GetTextExtent(self._label + "  ")
            dc.SetFont(self._bold)
            w2, _ = dc.GetTextExtent(self._value or "0")
            s = wx.Size(min(w1 + w2 + self.FromDIP(22), self.FromDIP(230)), self.FromDIP(26))
            self.SetMinSize(s)
            self.SetInitialSize(s)

        def SetValue(self, v):
            v = (v or "").strip()
            self._value = v
            previous = self.GetToolTipText()
            if previous != self._value_tooltip:
                self._hint_tooltip = previous
            self._value_tooltip = ((getattr(self, "_hint_tooltip", "") + "\n") if getattr(self, "_hint_tooltip", "") else "") + self._label + ": " + v
            self.SetToolTip(self._value_tooltip)
            self._fit()
            self.Show(bool(v))
            self.GetParent().Layout()
            self.Refresh()

        def _on_paint(self, e):
            dc = wx.AutoBufferedPaintDC(self)
            dc.SetBackground(wx.Brush(self.GetParent().GetBackgroundColour()))
            dc.Clear()
            gc = crisp(dc)
            w, h = self.GetClientSize()
            gc.SetBrush(wx.Brush(wx.Colour(C["tint"])))
            gc.SetPen(wx.Pen(wx.Colour(C["line"]), 1))
            gc.DrawRoundedRectangle(0.5, 0.5, w - 1, h - 1, (h - 1) / 2.0)
            f = self.GetFont()
            gc.SetFont(f, wx.Colour(C["muted"]))
            t1 = self._label + "  "
            w1, th = gc.GetTextExtent(t1)
            x = self.FromDIP(11)
            gc.DrawText(t1, x, (h - th) / 2.0)
            gc.SetFont(self._bold, wx.Colour(C["text"]))
            value = self._value
            room = max(0, w - x - w1 - self.FromDIP(10))
            if gc.GetTextExtent(value)[0] > room:
                lo, hi = 0, len(value)
                while lo < hi:
                    mid = (lo + hi + 1) // 2
                    if gc.GetTextExtent(value[:mid] + "…")[0] <= room:
                        lo = mid
                    else:
                        hi = mid - 1
                value = value[:lo] + "…" if room > gc.GetTextExtent("…")[0] else ""
            gc.DrawText(value, x + w1, (h - th) / 2.0)

    class ModernToolbar(wx.Panel):
        def __init__(self, parent, responsive=False):
            wx.Panel.__init__(self, parent, style=wx.BORDER_NONE)
            self.SetBackgroundColour(C["panel"])
            self.responsive = responsive
            self._fit_pending = False
            self.sizer = wx.WrapSizer(wx.HORIZONTAL) if responsive else wx.BoxSizer(wx.HORIZONTAL)
            if responsive:
                # the rows of buttons with the same margin above and below: the WrapSizer put them at the top
                # edge (the view switch touched it) and left all the room below
                self._pad = self.FromDIP(8)
                outer = wx.BoxSizer(wx.VERTICAL)
                outer.Add(self.sizer, 0, wx.EXPAND | wx.TOP, self._pad)
                self.SetSizer(outer)
            else:
                self.SetSizer(self.sizer)
            self.Bind(wx.EVT_PAINT, self._on_paint)
            self.SetMinSize(wx.Size(-1, self.FromDIP(50)))
            if responsive:
                self.Bind(wx.EVT_SIZE, self._resized)

        def stretch(self):
            if not self.responsive:
                self.sizer.AddStretchSpacer(1)

        def Layout(self):
            result = wx.Panel.Layout(self)
            self._queue_fit()
            return result

        def _resized(self, event):
            self._queue_fit()
            event.Skip()

        def _queue_fit(self):
            if self.responsive and not self._fit_pending:
                self._fit_pending = True
                wx.CallAfter(self._fit_rows)

        def _fit_rows(self):
            self._fit_pending = False
            if not self or self.GetClientSize().width <= 0:
                return
            wx.Panel.Layout(self)
            bottom = max((it.GetRect().bottom for it in self.sizer.GetChildren() if it.IsShown()), default=0)
            height = max(self.FromDIP(50), bottom + self._pad)
            if self.GetMinSize().height != height:
                self.SetMinSize(wx.Size(-1, height))
                self.GetParent().Layout()

        def add(self, win, border=2):
            self.sizer.Add(win, 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT | wx.RIGHT, self.FromDIP(border))

        def separator(self):
            s = wx.Panel(self, size=wx.Size(1, self.FromDIP(24)))
            s.SetBackgroundColour(C["line"])
            self.sizer.Add(s, 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT | wx.RIGHT, self.FromDIP(6))

        def _on_paint(self, e):
            dc = wx.PaintDC(self)
            w, h = self.GetClientSize()
            dc.SetPen(wx.Pen(wx.Colour(C["line"]), 1))
            dc.DrawLine(0, h - 1, w, h - 1)

    return FlatButton, Chip, ModernToolbar


_CLS = {}


def _cls(name):
    if not _CLS:
        FlatButton, Chip, ModernToolbar = _controls()
        _CLS.update(FlatButton=FlatButton, Chip=Chip, ModernToolbar=ModernToolbar)
    return _CLS[name]


# --------------------------------------------------------------------------
# 2. section headers (fold panel captions)
# --------------------------------------------------------------------------
def _group_colour(style):
    try:
        c = style.GetFirstColour()
        r, g, b = c.Red(), c.Green(), c.Blue()
    except Exception:
        return None
    if b > r + 8 and b > g + 8:
        return GROUP["blue"]
    if r > b + 8 and g > b + 8:
        return GROUP["yellow"]
    if r > g + 8 and r > b + 8:
        return GROUP["red"]
    if g > r + 8 and g > b + 8:
        return GROUP["green"]
    return None


def _install_captions():
    import wx
    import wx.lib.agw.foldpanelbar as fpb
    cb = fpb.CaptionBar
    orig_best = cb.DoGetBestSize
    orig_mouse = cb.OnMouseEvent

    def DoGetBestSize(self):
        s = orig_best(self)
        try:
            if self.IsVertical():
                return wx.Size(s.width, max(s.height, self.FromDIP(30)))
        except Exception:
            pass
        return s

    def OnPaint(self, event):
        if not getattr(self, "_controlCreated", False):
            event.Skip()
            return
        dc = wx.AutoBufferedPaintDC(self)
        dc.SetBackground(wx.Brush(wx.Colour(C["panel"])))
        dc.Clear()
        if not self.IsVertical():
            dc.DrawRotatedText(self._caption, 4, self.GetClientSize()[1] - 4, 90)
            return
        gc = crisp(dc)
        w, h = self.GetClientSize()
        caption = RENAME.get(self._caption or "", self._caption or "")
        main = not caption.strip().lower().startswith("additional")
        colour = _group_colour(self._style) if main else None
        gc.SetPen(wx.Pen(wx.Colour(C["line"]), 1))
        gc.StrokeLine(0, h - 0.5, w, h - 0.5)
        if colour:
            gc.SetPen(wx.TRANSPARENT_PEN)
            gc.SetBrush(wx.Brush(wx.Colour(colour)))
            bh = min(h - 8, self.FromDIP(16))
            gc.DrawRoundedRectangle(0, (h - bh) / 2.0, self.FromDIP(3), bh, 1)
        f = ui_font(9.5, 700) if main else ui_font(9, 500)
        gc.SetFont(f, wx.Colour(C["text"] if main else C["muted"]))
        tw, th = gc.GetTextExtent(caption)
        gc.DrawText(caption, self.FromDIP(9), (h - th) / 2.0)
        # chevron
        cx, cy, s = w - self.FromDIP(16), h / 2.0, self.FromDIP(4)
        path = gc.CreatePath()
        if self._collapsed:
            path.MoveToPoint(cx - s / 2.0, cy - s)
            path.AddLineToPoint(cx + s / 2.0, cy)
            path.AddLineToPoint(cx - s / 2.0, cy + s)
        else:
            path.MoveToPoint(cx - s, cy - s / 2.0)
            path.AddLineToPoint(cx, cy + s / 2.0)
            path.AddLineToPoint(cx + s, cy - s / 2.0)
        gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(wx.Colour(C["muted"])).Width(1.4 * _scale(self))))
        gc.StrokePath(path)

    def OnMouseEvent(self, event):
        # a click anywhere on the header opens or closes the section
        if event.LeftDown():
            evt = fpb.CaptionBarEvent(fpb.wxEVT_CAPTIONBAR)
            evt.SetId(self.GetId())
            evt.SetEventObject(self)
            evt.SetBar(self)
            self.GetEventHandler().ProcessEvent(evt)
            return
        if event.Entering() or event.Moving():
            self.SetCursor(wx.Cursor(wx.CURSOR_HAND))
            return
        if event.LeftDClick():
            return  # the first click already toggled
        return orig_mouse(self, event)

    def init_wrap(orig_init):
        def __init__(self, *args, **kwargs):
            orig_init(self, *args, **kwargs)
            try:
                self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
            except Exception:
                pass
        return __init__

    cb.__init__ = init_wrap(cb.__init__)
    cb.DoGetBestSize = DoGetBestSize
    cb.OnPaint = OnPaint
    cb.OnMouseEvent = OnMouseEvent


# --------------------------------------------------------------------------
# 3. generic window theming
# --------------------------------------------------------------------------
_NATIVE_KEEP = ("Button", "BitmapButton", "ToggleButton", "TextCtrl", "ComboBox", "Choice", "SpinCtrl",
                "SpinCtrlDouble", "ListCtrl", "ListBox", "CheckListBox", "TreeCtrl", "Notebook",
                "FlatButton", "Chip", "GlassLauncher", "ModernToolbar")


def _whiten(w):
    import wx
    for c in w.GetChildren():
        if isinstance(c, wx.TopLevelWindow) or _skip_scaling(c):
            continue
        name = c.__class__.__name__
        if name not in _NATIVE_KEEP and not isinstance(c, (wx.Button, wx.TextCtrl, wx.ComboBox, wx.Choice)):
            try:
                c.SetBackgroundColour(C["panel"])
            except Exception:
                pass
        _whiten(c)


def theme_tree(root):
    import wx
    import wx.lib.agw.foldpanelbar as fpb

    def walk(w):
        for c in w.GetChildren():
            if isinstance(c, wx.TopLevelWindow) or _skip_scaling(c):
                continue
            try:
                if isinstance(c, fpb.FoldPanelBar):
                    c.SetBackgroundColour(C["panel"])
                    _whiten(c)
                    try:
                        c.GetParent().SetBackgroundColour(C["panel"])
                    except Exception:
                        pass
                    continue
                if isinstance(c, wx.SplitterWindow):
                    st = c.GetWindowStyleFlag() & ~(wx.SP_3D | wx.SP_3DSASH | wx.SP_3DBORDER | wx.SP_BORDER)
                    c.SetWindowStyleFlag(st | wx.SP_LIVE_UPDATE | wx.SP_NOBORDER | wx.SP_THIN_SASH)
                    c.SetBackgroundColour(C["line2"])
                    c.Refresh()
            except Exception as e:
                _log("theme:", e)
            walk(c)

    walk(root)


def _find(root, cls):
    import wx
    out = []

    def walk(w):
        for c in w.GetChildren():
            if isinstance(c, wx.TopLevelWindow):
                continue
            if isinstance(c, cls):
                out.append(c)
            walk(c)

    walk(root)
    return out


def relayout_foldpanels(root):
    """Re-stack fold panel contents after sizes changed (taller headers,
    scaled controls). The fold panel bar positions its children once, when
    they are added, so they have to be moved here."""
    import wx.lib.agw.foldpanelbar as fpb
    for bar in _find(root, fpb.FoldPanelBar):
        try:
            for i in range(bar.GetCount()):
                item = bar.GetFoldPanel(i)
                if not item.IsVertical():
                    continue
                cb = item._captionBar
                cb.SetSize(cb.GetSize().width, cb.FromDIP(28))
                y = item.GetCaptionLength()
                for wi in item._items:
                    w = wi._wnd
                    if w is not None:
                        wi._leftSpacing = item.FromDIP(6)
                        wi._rightSpacing = item.FromDIP(6)
                        sz = w.GetSizer()
                        if sz is not None:
                            best = sz.GetMinSize()
                            w.SetSize(max(w.GetSize().width, best.width), best.height)
                        w.SetPosition((wi._leftSpacing, y + wi._spacing))
                    y += wi.GetWindowLength(True)
                item._LastInsertPos = y
                item.ResizePanel()
            bar.RedisplayFoldPanelItems()
            if bar.GetCount():
                bar.RefreshPanelsFrom(bar.GetFoldPanel(0))
        except Exception as e:
            _log("fold panel layout:", e)


def _after_show_splitters(frame):
    """Keep the right-hand panes (peak list, controls) at their intended
    96-dpi widths on scaled displays."""
    import wx
    f = _scale(frame)
    if f <= 1.01:
        return
    for name in ("splitterwindow", "splitterwindow2"):
        sp = getattr(frame, name, None)
        if not isinstance(sp, wx.SplitterWindow) or not sp.IsSplit():
            continue
        try:
            w = sp.GetClientSize()[0]
            right = w - sp.GetSashPosition()
            if 0 < right < w:
                sp.SetSashPosition(max(50, w - int(right * f)))
        except Exception:
            pass


# --------------------------------------------------------------------------
# 4. main UniDec window: toolbar, buttons, cards
# --------------------------------------------------------------------------
def _replace_button(btn, kind, label=None, height=None, padx=12, min_width=0):
    """Swap a native button for a flat one with the same ID and place.
    Without `height` the original pixel height is kept (fold panels do not
    re-layout when a child grows)."""
    FlatButton = _cls("FlatButton")
    parent = btn.GetParent()
    tip = btn.GetToolTipText() if btn.GetToolTip() else None
    if height is None:
        hpx = max(btn.GetSize().height, btn.GetBestSize().height)
        height = hpx / max(_scale(btn), 1.0)
    new = FlatButton(parent, btn.GetId(), label or btn.GetLabel(), kind=kind, tooltip=tip, height=height,
                     padx=padx, min_width=min_width)
    sizer = btn.GetContainingSizer()
    if sizer is not None and sizer.Replace(btn, new, True):
        btn.Hide()
        return new
    new.Destroy()
    return None


def _build_toolbar(frame):
    import wx
    ctl = frame.controls
    FlatButton, Chip, ModernToolbar = _cls("FlatButton"), _cls("Chip"), _cls("ModernToolbar")
    tb = ModernToolbar(frame)
    tb.sizer.AddSpacer(tb.FromDIP(8))

    def mk(attr, label, icon, kind="ghost", tip=None):
        orig = getattr(ctl, attr, None)
        if orig is None:
            return None
        t = (orig.GetToolTipText() if orig.GetToolTip() else None) or tip
        b = FlatButton(tb, orig.GetId(), label, icon=icon, kind=kind, tooltip=t, height=34)
        tb.add(b)
        return b

    mk("openbutton", "Open", "open", tip="Open a data file (Ctrl+O)")
    tb.separator()
    mk("procbutton", "Process", "process", tip="Process data (Ctrl+D)")
    mk("udbutton", "Deconvolute", "deconvolve", tip="Run the deconvolution (Ctrl+R)")
    mk("ppbutton", "Pick peaks", "pick", tip="Peak detection (Ctrl+P)")
    tb.separator()
    mk("autobutton", "Run all", "run", kind="primary", tip="Process, deconvolute and pick peaks (Ctrl+E)")
    mk("replotbutton", "Replot", "replot", tip="Replot (Ctrl+N)")
    tb.sizer.AddStretchSpacer(1)
    chips = {}
    for key, label in (("file", "File"), ("points", "Points"), ("r2", "R\u00b2"), ("score", "UniScore")):
        ch = Chip(tb, label)
        tb.add(ch, 3)
        chips[key] = ch
    tb.sizer.AddSpacer(tb.FromDIP(8))
    frame._udp_chips = chips

    # hide the old icon row at the top of the control panel
    for attr in ("openbutton", "procbutton", "udbutton", "ppbutton", "autobutton"):
        b = getattr(ctl, attr, None)
        if b is not None:
            b.Hide()
    ctl.Layout()

    # toolbar above the window content
    old = frame.GetSizer()
    v = wx.BoxSizer(wx.VERTICAL)
    v.Add(tb, 0, wx.EXPAND)
    if old is not None:
        items = list(old.GetChildren())
        if len(items) == 1 and items[0].IsWindow():
            win = items[0].GetWindow()
            old.Detach(win)
            v.Add(win, 1, wx.EXPAND)
            frame.SetSizer(v, True)
        else:
            v.Add(old, 1, wx.EXPAND)
            frame.SetSizer(v, False)
    else:
        # a frame without sizer lays out its single child itself; with the
        # toolbar added it needs a sizer
        for c in frame.GetChildren():
            if c is tb or isinstance(c, (wx.TopLevelWindow, wx.StatusBar, wx.ToolBar)):
                continue
            v.Add(c, 1, wx.EXPAND)
        frame.SetSizer(v)
    frame._udp_toolbar = tb

    # mirror status bar values in the chips
    orig_set = frame.SetStatusText

    def SetStatusText(text, number=0, *a, **k):
        r = orig_set(text, number, *a, **k)
        try:
            _update_chip(frame, str(text))
        except Exception:
            pass
        return r

    frame.SetStatusText = SetStatusText


def _update_chip(frame, text):
    chips = getattr(frame, "_udp_chips", None)
    if not chips:
        return
    t = text.strip()
    if t.startswith("File:"):
        chips["file"].SetValue(t[5:].strip())
    elif t.startswith("Data Length:"):
        chips["points"].SetValue(t[12:].strip())
    elif t.startswith("R\u00b2"):
        v = t[2:].lstrip(":").strip()
        try:
            v = "%.3f" % float(v)
        except Exception:
            pass
        chips["r2"].SetValue(v)
    elif t.lower().startswith("uniscore"):
        chips["score"].SetValue(t.split(":", 1)[-1].strip())


def _restyle_controls(frame):
    ctl = frame.controls
    for attr in ("dataprepbutton", "rununidec", "plotbutton"):
        b = getattr(ctl, attr, None)
        if b is not None:
            _replace_button(b, "primary", label=RENAME.get(b.GetLabel(), b.GetLabel()))
    b = getattr(ctl, "plotbutton2", None)
    if b is not None:
        _replace_button(b, "secondary")
    names = {"mainbutton": "Main", "expandbutton": "All", "collapsebutton": "None",
             "bluebutton": "Processing", "yellowbutton": "Deconvolution", "redbutton": "Peaks"}
    for attr, label in names.items():
        b = getattr(ctl, attr, None)
        if b is not None:
            _replace_button(b, "ghost", label=label, height=26, padx=6)
    try:
        ctl.SetBackgroundColour(C["panel"])
        ctl.Layout()
    except Exception:
        pass


def _make_cards(frame):
    import wx
    if getattr(frame, "tabbed", 0) == 1 or getattr(frame.config, "imflag", 0) != 0:
        return
    gbs = getattr(frame, "sizerplot", None)
    panel = getattr(frame, "plotpanel", None)
    if gbs is None or panel is None:
        return
    header, pad, gap = panel.FromDIP(30), panel.FromDIP(6), panel.FromDIP(12)
    cards = []
    for name, title in PLOT_TITLES:
        plot = getattr(frame, name, None)
        if plot is None:
            continue
        item = gbs.FindItem(plot)
        if item is None:
            continue
        pos, span = item.GetPos(), item.GetSpan()
        gbs.Detach(plot)
        box = wx.BoxSizer(wx.VERTICAL)
        box.AddSpacer(header)
        box.Add(plot, 1, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, pad)
        gbs.Add(box, pos, span, flag=wx.EXPAND | wx.ALL, border=gap // 2)
        cards.append((plot, title))
    if not cards:
        return
    outer = wx.BoxSizer(wx.VERTICAL)
    outer.Add(gbs, 1, wx.EXPAND | wx.ALL, gap // 2)
    panel.SetSizer(outer, False)
    panel.SetBackgroundColour(C["bg"])
    frame._udp_cards = {"header": header, "pad": pad, "gap": gap, "cards": cards}
    panel.Bind(wx.EVT_PAINT, lambda e, f=frame: _paint_cards(f))
    try:
        panel.SetupScrolling(scroll_x=False, scrollToTop=False)
    except Exception:
        pass
    frame.resize_plots()


def _paint_cards(frame):
    import wx
    panel = frame.plotpanel
    dc = wx.PaintDC(panel)
    c = getattr(frame, "_udp_cards", None)
    if not c:
        return
    gc = crisp(dc)
    font = ui_font(10, 700)
    for plot, title in c["cards"]:
        try:
            if not plot.IsShown():
                continue
            r = plot.GetRect()
        except RuntimeError:
            continue
        x, y = r.x - c["pad"], r.y - c["header"]
        w, h = r.width + 2 * c["pad"], r.height + c["header"] + c["pad"]
        gc.SetPen(wx.Pen(wx.Colour(C["line"]), 1))
        gc.SetBrush(wx.Brush(wx.Colour(C["panel"])))
        gc.DrawRoundedRectangle(x + 0.5, y + 0.5, w - 1, h - 1, panel.FromDIP(8))
        gc.SetFont(font, wx.Colour(C["text"]))
        tw, th = gc.GetTextExtent(title)
        gc.DrawText(title, x + panel.FromDIP(12), y + (c["header"] - th) / 2.0 + panel.FromDIP(2))


def _install_resize_plots():
    from unidec.modules.gui_elements.mainwindow_base import MainwindowBase
    import wx
    orig = MainwindowBase.resize_plots

    def resize_plots(self, e=None):
        c = getattr(self, "_udp_cards", None)
        if not c:
            return orig(self, e)
        try:
            panel = self.plotpanel
            cw = panel.GetClientSize()[0]
            w = int((cw - 3 * c["gap"]) / 2 - 2 * c["pad"])
            if w < 50:
                return
            try:
                fs = self.config.figsize
                aspect = float(fs[1]) / float(fs[0])
            except Exception:
                aspect = 0.7
            h = max(int(w * aspect), panel.FromDIP(240))
            if c.get("last") == (w, h):
                return
            c["last"] = (w, h)
            for plot, _t in c["cards"]:
                plot.SetMinSize(wx.Size(w, h))
            panel.Layout()
            panel.FitInside()
            panel.Refresh()
        except RuntimeError:
            pass

    MainwindowBase.resize_plots = resize_plots


# --------------------------------------------------------------------------
# 5. peak list
# --------------------------------------------------------------------------
def _marker_bitmap(colour, marker, n):
    import wx
    bmp = wx.Bitmap.FromRGBA(n, n, 0, 0, 0, 0)
    dc = wx.MemoryDC(bmp)
    gc = wx.GraphicsContext.Create(dc)
    gc.SetBrush(wx.Brush(colour))
    gc.SetPen(wx.Pen(wx.Colour(255, 255, 255), 1))
    m = 1.0
    s = n - 2 * m
    path = gc.CreatePath()
    if marker == "o":
        path.AddCircle(n / 2.0, n / 2.0, s / 2.0)
    elif marker == "s":
        path.AddRectangle(m + 1, m + 1, s - 2, s - 2)
    elif marker == "d":
        path.MoveToPoint(n / 2.0, m)
        path.AddLineToPoint(n - m, n / 2.0)
        path.AddLineToPoint(n / 2.0, n - m)
        path.AddLineToPoint(m, n / 2.0)
    elif marker == "^":
        path.MoveToPoint(n / 2.0, m)
        path.AddLineToPoint(n - m, n - m)
        path.AddLineToPoint(m, n - m)
    elif marker == ">":
        path.MoveToPoint(m, m)
        path.AddLineToPoint(n - m, n / 2.0)
        path.AddLineToPoint(m, n - m)
    elif marker == "*":
        import math
        for i in range(10):
            r = s / 2.0 if i % 2 == 0 else s / 4.5
            a = -math.pi / 2 + i * math.pi / 5
            px, py = n / 2.0 + r * math.cos(a), n / 2.0 + r * math.sin(a)
            (path.MoveToPoint if i == 0 else path.AddLineToPoint)(px, py)
    else:  # "v" and anything else
        path.MoveToPoint(m, m)
        path.AddLineToPoint(n - m, m)
        path.AddLineToPoint(n / 2.0, n - m)
    path.CloseSubpath()
    gc.DrawPath(path)
    del gc
    dc.SelectObject(wx.NullBitmap)
    return bmp


def _install_peaklist():
    import wx
    from unidec.modules.gui_elements import peaklistsort as pls
    cls = pls.PeakListCtrlPanel
    orig_init = cls.__init__
    orig_add = cls.add_data

    def __init__(self, *args, **kwargs):
        orig_init(self, *args, **kwargs)
        try:
            lc = self.list_ctrl
            st = lc.GetWindowStyleFlag() & ~wx.BORDER_MASK
            lc.SetWindowStyleFlag(st | wx.BORDER_NONE)
            f = _scale(self)
            if f > 1.01:
                for col in (0, 2, 4):
                    lc.SetColumnWidth(col, int(lc.GetColumnWidth(col) * f))
            self.SetBackgroundColour(C["panel"])
            head = wx.Panel(self)
            head.SetBackgroundColour(C["panel"])
            hs = wx.BoxSizer(wx.HORIZONTAL)
            t = wx.StaticText(head, label="Peaks")
            t.SetFont(ui_font(10, 700))
            t.SetForegroundColour(C["text"])
            cnt = wx.StaticText(head, label="")
            cnt.SetFont(ui_font(9, 500))
            cnt.SetForegroundColour(C["muted"])
            hs.Add(t, 0, wx.ALIGN_CENTER_VERTICAL)
            hs.AddStretchSpacer(1)
            hs.Add(cnt, 0, wx.ALIGN_CENTER_VERTICAL)
            head.SetSizer(hs)
            self.GetSizer().Insert(0, head, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP | wx.BOTTOM,
                                   self.FromDIP(8))
            self._udp_count = cnt
            self.Layout()
        except Exception as e:
            _log("peak list:", e)

    def add_data(self, pks, *args, **kwargs):
        result = orig_add(self, pks, *args, **kwargs)
        try:
            lc = self.list_ctrl
            n = self.FromDIP(12)
            il = wx.ImageList(n, n)
            row = 0
            for p in pks.peaks:
                if getattr(p, "ignore", 0):
                    continue
                if row >= lc.GetItemCount():
                    break
                col = wx.Colour(int(round(p.color[0] * 255)), int(round(p.color[1] * 255)),
                                int(round(p.color[2] * 255)))
                idx = il.Add(_marker_bitmap(col, getattr(p, "marker", "o"), n))
                lc.SetItemBackgroundColour(row, wx.Colour(C["panel"]))
                lc.SetItemTextColour(row, wx.Colour(C["text"]))
                lc.SetItemText(row, "")
                lc.SetItemImage(row, idx)
                row += 1
            lc.AssignImageList(il, wx.IMAGE_LIST_SMALL)
            f = _scale(self)
            if f > 1.01:
                lc.SetColumnWidth(3, int(65 * f))
            cnt = getattr(self, "_udp_count", None)
            if cnt is not None:
                cnt.SetLabel("%d found" % row if row != 1 else "1 found")
                cnt.GetParent().Layout()
        except Exception as e:
            _log("peak list rows:", e)
        return result

    cls.__init__ = __init__
    cls.add_data = add_data


# --------------------------------------------------------------------------
# 6. window geometry and last folder
# --------------------------------------------------------------------------
def _restore_geometry(frame):
    import wx
    geo = (_load().get("windows") or {}).get(frame.__class__.__name__)
    if not geo:
        return False
    try:
        x, y, w, h = [int(v) for v in geo["rect"]]
        if w < 300 or h < 200:
            return False
        # saved on a display that is gone or smaller now (a large monitor, then
        # only the laptop): onto the display nearest to it, not larger than it
        n = wx.Display.GetFromPoint(wx.Point(x + w // 2, y + 20))  # the middle of its title bar
        if n == wx.NOT_FOUND:
            n = wx.Display.GetFromPoint(wx.GetMousePosition())
        a = wx.Display(n if n != wx.NOT_FOUND else 0).GetClientArea()
        x, y, w, h = fit_rect((x, y, w, h), (a.x, a.y, a.width, a.height), gap=frame.FromDIP(4))
        frame.SetSize(x, y, w, h)
        if geo.get("maximized"):
            frame.Maximize(True)
        return True
    except Exception:
        return False


def _track_geometry(frame):
    import wx
    state = {"pending": None}

    def save():
        state["pending"] = None
        try:
            if not frame or frame.IsIconized():
                return
            data = _load()
            wins = data.get("windows") or {}
            entry = wins.get(frame.__class__.__name__, {})
            entry["maximized"] = bool(frame.IsMaximized())
            if not frame.IsMaximized():
                r = frame.GetRect()
                entry["rect"] = [r.x, r.y, r.width, r.height]
            if "rect" in entry:
                wins[frame.__class__.__name__] = entry
                _save({"windows": wins})
        except RuntimeError:
            pass
        except Exception as e:
            _log("geometry:", e)

    def later(evt):
        evt.Skip()
        if state["pending"] is None:
            state["pending"] = wx.CallLater(800, save)
        else:
            state["pending"].Restart(800)

    frame.Bind(wx.EVT_SIZE, later)
    frame.Bind(wx.EVT_MOVE, later)

    def on_close(evt):
        save()
        evt.Skip()

    frame.Bind(wx.EVT_CLOSE, on_close)


def _install_launch():
    import wx
    from unidec.modules.gui_elements.mainwindow_base import MainwindowBase

    def launch(self):
        try:
            apply_main_theme(self)
        except Exception as e:
            _log("theme failed:", e)
        if not _restore_geometry(self):
            self.Centre()
        if _load().get("start_maximized", True):
            self.Maximize(True)  # opens maximised; the stored size is used when it is restored
        self.Show(True)
        try:
            _track_geometry(self)
        except Exception as e:
            _log("geometry tracking:", e)
        wx.CallAfter(_after_show_splitters, self)

    MainwindowBase.launch = launch


def _install_last_folder():
    import wx
    last = _load().get("last_folder")
    try:
        from unidec.modules.isolated_packages import FileDialogs
        if last and os.path.isdir(last):
            FileDialogs.default_dir = last
    except Exception:
        FileDialogs = None
    import importlib
    try:
        g = importlib.import_module("unidec.GUniDec")
    except Exception as e:
        _log("last folder:", e)
        return
    app_cls = g.UniDecApp
    orig_open_file = app_cls.on_open_file

    def on_open_file(self, filename, directory=None, *args, **kwargs):
        try:
            if directory and os.path.isdir(directory):
                _save({"last_folder": directory})
                if FileDialogs is not None:
                    FileDialogs.default_dir = directory
        except Exception:
            pass
        # LabSolutions .lcd files are read by LCMS Analysis (UniDec cannot):
        # open them there instead of failing with "Unsupported file type"
        try:
            full = os.path.join(directory, filename) if directory else filename
            if str(filename).lower().endswith(".lcd") and os.path.isfile(full):
                import unilcms
                unilcms.open_window(full)
                print("Opened in LCMS Analysis:", full)
                return None
        except Exception as e:
            _log("open .lcd in LCMS Analysis:", e)
        return orig_open_file(self, filename, directory, *args, **kwargs)

    def on_open(self, e=None):
        start = _load().get("last_folder") or ""
        if start and not os.path.isdir(start):
            start = ""
        dlg = wx.FileDialog(self.view, "Open a data file", start, "",
                            "*.*")
        try:
            if dlg.ShowModal() == wx.ID_OK:
                self.view.SetStatusText("Opening", number=5)
                filename = dlg.GetFilename()
                print("Opening: ", filename)
                if os.path.splitext(filename)[1] == ".zip":
                    print("Can't open zip, try Load State.")
                    return
                self.on_open_file(filename, dlg.GetDirectory())
        finally:
            dlg.Destroy()

    app_cls.on_open_file = on_open_file
    app_cls.on_open = on_open


# --------------------------------------------------------------------------
# 7. launcher (light frosted glass)
# --------------------------------------------------------------------------
APP_NAME = "MSpektra"
APP_VERSION = "4.23"  # +0.01 small change, +0.1 large change, +1.0 big change
WORKSPACES = [
    ("LCMS Analysis", "LC-MS and HPLC files", "lcms"),
    ("HRMS Analysis", "Bruker .d and mzML files", "hrms"),
    ("Deconvolute", "Text, JCAMP-DX, mzML, Thermo and Waters spectra", "deconv"),
]
# visible texts of the deconvolution window, renamed (UniDec's code is unchanged)
RENAME = {"UniDec Parameters": "Deconvolution Parameters", "Run UniDec": "Run Deconvolution", "UniDec": "Bayesian deconvolution"}
_LAUNCH_ORDER = ["UniDec", "UniLCMS", "MetaUniDec", "UniDecCD", "UniChromCD", "IsoDec",
                 "UniChrom", "Data Collector", "UltraMeta Data Collector", "Import Wizard", "HDF5 Import Wizard",
                 "UniDec Processing HEKPipeline", "UniDec API Shell"]
_DECONV = {"UniDec", "UniLCMS", "MetaUniDec", "UniDecCD", "UniChromCD", "IsoDec"}

# glyphs in a 24 x 24 box: ("line", [points]), ("curve", [start, (c1, c2, end), ...]),
# ("dots", [centres]), ("bars", [(x, top)]), ("rect", (x, y, w, h)), ("ellipse", (x, y, w, h))
_GLYPHS = {
    "UniDec": [("line", [(3, 20), (5, 20), (6, 14), (7, 20), (10, 20), (12, 4), (14, 20), (17, 20), (18, 12),
                         (19, 20), (21, 20)])],
    "MetaUniDec": [("line", [(3, 8), (7, 8), (8.5, 3.5), (10, 8), (21, 8)]),
                   ("line", [(3, 13.5), (9, 13.5), (10.5, 9), (12, 13.5), (21, 13.5)]),
                   ("line", [(3, 19), (11, 19), (12.5, 14.5), (14, 19), (21, 19)])],
    "UniChrom": [("line", [(3, 20), (21, 20)]),
                 ("curve", [(3, 19), ((7, 19), (7, 5), (10.5, 5)), ((14, 5), (13.5, 19), (16.5, 17)),
                            ((18, 16), (18.5, 19), (21, 19))])],
    "UniLCMS": [("line", [(3, 20), (21, 20)]),
                ("curve", [(3, 19.5), ((6.5, 19.5), (6.5, 5), (9, 5)), ((11.5, 5), (11.5, 19.5), (14, 19.5))]),
                ("bars", [(17, 12), (20, 8)])],
    "LCMS Analysis": [("line", [(3, 20), (21, 20)]),
                     ("curve", [(3, 19.5), ((6.5, 19.5), (6.5, 5), (9, 5)), ((11.5, 5), (11.5, 19.5), (14, 19.5))]),
                     ("bars", [(17, 12), (20, 8)])],
    "HRMS Analysis": [("line", [(3, 20), (21, 20)]), ("bars", [(6, 6), (9.5, 9), (13, 14), (16.5, 17)]),
                     ("line", [(18, 5), (21, 5)]), ("line", [(19.5, 3.5), (19.5, 6.5)])],
    "Deconvolute": [("line", [(3, 20), (21, 20)]), ("bars", [(5, 15), (8, 10.5), (11, 9), (14, 13)]),
                    ("line", [(16.5, 8), (19.5, 8)]), ("line", [(18, 6.5), (19.5, 8), (18, 9.5)]),
                    ("line", [(20.5, 4), (20.5, 12)])],
    "app": [("line", [(3.5, 19.5), (20.5, 19.5)]), ("bars", [(6, 16), (9.5, 10), (13, 5), (16.5, 9), (20, 15)])],
    "UniDecCD": [("line", [(3, 3), (3, 21), (21, 21)]),
                 ("dots", [(7, 16), (9.5, 11), (12, 15), (14, 8), (17, 12), (19, 6), (8, 6)])],
    "UniChromCD": [("line", [(3, 21), (21, 21)]),
                   ("curve", [(3, 20), ((8, 20), (8, 9), (12, 9)), ((16, 9), (16, 20), (21, 20))]),
                   ("dots", [(8, 5), (12, 3.5), (16, 5.5), (11, 13), (14, 14)])],
    "IsoDec": [("line", [(3, 21), (21, 21)]), ("bars", [(5, 14), (8.5, 7), (12, 4), (15.5, 9), (19, 15)])],
    "Data Collector": [("line", [(3, 3), (3, 21), (21, 21)]), ("line", [(6, 17), (10, 11), (14, 14), (19, 6)]),
                       ("dots", [(6, 17), (10, 11), (14, 14), (19, 6)])],
    "UltraMeta Data Collector": [("rect", (3, 3, 8, 8)), ("rect", (13, 3, 8, 8)), ("rect", (3, 13, 8, 8)),
                                 ("rect", (13, 13, 8, 8)), ("line", [(5, 9), (7, 6), (9, 7.5)]),
                                 ("line", [(15, 9), (17, 5), (19, 8)]), ("line", [(5, 19), (7, 16), (9, 17)]),
                                 ("line", [(15, 18), (17, 15.5), (19, 15)])],
    "Import Wizard": [("line", [(6, 3), (14, 3), (18, 7), (18, 21), (6, 21), (6, 3)]),
                      ("line", [(12, 9), (12, 17)]), ("line", [(9, 14), (12, 17), (15, 14)])],
    "HDF5 Import Wizard": [("ellipse", (5, 3, 14, 5)), ("line", [(5, 5.5), (5, 18.5)]),
                           ("line", [(19, 5.5), (19, 18.5)]), ("ellipse", (5, 16, 14, 5)),
                           ("ellipse", (5, 9.5, 14, 5))],
    "UniDec Processing HEKPipeline": [("rect", (2, 9, 5, 6)), ("rect", (9.5, 9, 5, 6)), ("rect", (17, 9, 5, 6)),
                                      ("line", [(7, 12), (9.5, 12)]), ("line", [(14.5, 12), (17, 12)])],
    "UniDec API Shell": [("rect", (3, 4, 18, 16)), ("line", [(7, 9.5), (10, 12), (7, 14.5)]),
                         ("line", [(12, 14.5), (16, 14.5)])],
}


class _Crisp(object):
    """A GraphicsContext for shapes, gradients and shadows whose text is
    drawn with the plain DC instead: GDI text with ClearType on whole pixels
    is far sharper than the anti-aliased text of the graphics renderer. The
    shapes drawn so far are flushed first, so the drawing order is kept."""

    def __init__(self, gc, dc):
        import wx
        self.__dict__.update(gc=gc, dc=dc, _font=None, _col=wx.Colour(0, 0, 0))

    def __getattr__(self, name):
        return getattr(self.__dict__["gc"], name)

    def SetFont(self, font, colour=None):
        import wx
        col = colour if colour is not None else wx.Colour(0, 0, 0)
        self.__dict__.update(_font=font, _col=col)
        self.dc.SetFont(font)
        try:
            self.gc.SetFont(font, col)  # for any measurement passed on to the context
        except Exception:
            pass

    def GetTextExtent(self, text):
        w, h = self.dc.GetTextExtent(text)
        return float(w), float(h)

    def DrawText(self, text, x, y, *args):
        import wx
        try:
            self.gc.Flush()
        except Exception:
            pass
        dc = self.dc
        if self._font is not None:
            dc.SetFont(self._font)
        dc.SetBackgroundMode(wx.TRANSPARENT)
        dc.SetTextForeground(self._col)
        dc.DrawText(text, int(round(x)), int(round(y)))

    def finish(self):
        try:
            self.gc.Flush()
        except Exception:
            pass


def crisp(dc):
    """GraphicsContext on dc with sharp (GDI) text; see _Crisp."""
    import wx
    return _Crisp(wx.GraphicsContext.Create(dc), dc)


def _glass_class():
    import wx

    NAVY = (14, 28, 48)

    class GlassLauncher(wx.Panel):
        """Custom painted launcher: pale blue-grey backdrop with a faint
        charge state series, white frosted glass tiles, Plus Jakarta Sans."""

        def __init__(self, parent, tools, logo_path, version, more=None):
            wx.Panel.__init__(self, parent, style=wx.BORDER_NONE | wx.WANTS_CHARS)
            self.more = more or []
            self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
            self.SetBackgroundColour(wx.Colour(C["bg"]))
            self.tools = tools
            self.version = version
            self._logo = None
            if logo_path and os.path.isfile(logo_path):
                try:
                    self._logo = wx.Image(logo_path)
                except Exception:
                    self._logo = None
            self._logo_bmp = {}
            self._rects = []
            self._hover = self._down = self._focus = -1
            self.k = self.FromDIP(1000) / 1000.0
            for ev, fn in ((wx.EVT_PAINT, self._on_paint), (wx.EVT_MOTION, self._on_motion),
                           (wx.EVT_LEAVE_WINDOW, self._on_leave), (wx.EVT_LEFT_DOWN, self._on_down),
                           (wx.EVT_LEFT_DCLICK, self._on_down), (wx.EVT_LEFT_UP, self._on_up),
                           (wx.EVT_SIZE, lambda e: self.Refresh()), (wx.EVT_KEY_DOWN, self._on_key),
                           (wx.EVT_RIGHT_UP, lambda e: self._more_menu(e.GetPosition()))):
                self.Bind(ev, fn)
            self.show_more = False  # classic tools: right click on the background only
            self.SetMinSize(self.FromDIP(wx.Size(980, 500)))

        # ---------------------------------------------------------- input
        def _hit(self, pos):
            for (x, y, w, h, i) in self._rects:
                if x <= pos.x <= x + w and y <= pos.y <= y + h:
                    return i
            return -1

        def _on_motion(self, e):
            i = self._hit(e.GetPosition())
            if i != self._hover:
                self._hover = i
                self.SetCursor(wx.Cursor(wx.CURSOR_HAND if i >= 0 else wx.CURSOR_ARROW))
                self.Refresh()

        def _on_leave(self, e):
            if self._hover != -1 or self._down != -1:
                self._hover = self._down = -1
                self.Refresh()

        def _on_down(self, e):
            self._down = self._hit(e.GetPosition())
            self.SetFocus()
            self.Refresh()

        def _on_up(self, e):
            i = self._hit(e.GetPosition())
            was, self._down = self._down, -1
            self.Refresh()
            if i >= 0 and i == was:
                self._activate(i)

        def _more_menu(self, at=None):
            if not self.more:
                return
            m = wx.Menu()
            head = m.Append(wx.ID_ANY, "More tools")
            head.Enable(False)
            m.AppendSeparator()
            for title, desc, tid in self.more:
                it = m.Append(wx.ID_ANY, title, desc)

                def fire(ev, tid=tid):
                    evt = wx.CommandEvent(wx.EVT_BUTTON.typeId, tid)
                    evt.SetEventObject(self)
                    wx.CallAfter(self.GetEventHandler().ProcessEvent, evt)
                self.Bind(wx.EVT_MENU, fire, it)
            r = [rc for rc in self._rects if rc[4] == len(self.tools)]
            if at is not None:
                pos = at
            else:
                pos = wx.Point(int(r[0][0]), int(r[0][1] + r[0][3])) if r else wx.DefaultPosition
            try:
                wx.ToolTip.Enable(False)
                self.PopupMenu(m, pos)
            finally:
                wx.ToolTip.Enable(True)
            m.Destroy()

        def _on_key(self, e):
            n = len(self.tools) + (1 if (self.more and self.show_more) else 0)
            key = e.GetKeyCode()
            f = self._focus if self._focus >= 0 else 0
            if key in (wx.WXK_RIGHT, wx.WXK_TAB):
                f = (f + 1) % n
            elif key == wx.WXK_LEFT:
                f = (f - 1) % n
            elif key == wx.WXK_DOWN:
                f = min(n - 1, f + 3)
            elif key == wx.WXK_UP:
                f = max(0, f - 3)
            elif key in (wx.WXK_RETURN, wx.WXK_NUMPAD_ENTER, wx.WXK_SPACE) and self._focus >= 0:
                self._activate(self._focus)
                return
            else:
                e.Skip()
                return
            self._focus = f
            self.Refresh()

        def _activate(self, i):
            if i == len(self.tools):
                self._more_menu()
                return
            evt = wx.CommandEvent(wx.EVT_BUTTON.typeId, self.tools[i][2])
            evt.SetEventObject(self)
            wx.CallAfter(self.GetEventHandler().ProcessEvent, evt)

        def _grow(self, h):
            self.SetMinSize(wx.Size(self.GetMinSize().width, h))
            top = self.GetTopLevelParent()
            top.Fit()
            top.Centre()
            fit_to_screen(top, grow=False)
            self.Refresh()

        # ---------------------------------------------------------- helpers
        def _rr(self, gc, x, y, w, h, r):
            gc.DrawRoundedRectangle(x, y, w, h, r)

        def _shadow(self, gc, x, y, w, h, r, lift=1.0):
            """Layered soft shadow in navy (approximates a CSS box-shadow)."""
            k = self.k
            gc.SetPen(wx.TRANSPARENT_PEN)
            for dy, spread, a in ((30 * lift, 20, 9), (16 * lift, 11, 13), (7 * lift, 5, 16), (2, 1.5, 18)):
                gc.SetBrush(wx.Brush(wx.Colour(NAVY[0], NAVY[1], NAVY[2], a)))
                self._rr(gc, x - spread * 0.35 * k, y + dy * 0.45 * k, w + spread * 0.7 * k,
                         h + spread * 0.3 * k, r + spread * 0.4 * k)

        def _glass(self, gc, x, y, w, h, r, state, focus=False):
            k = self.k
            x, y, w, h = round(x), round(y), round(w), round(h)
            if state == 1:
                y -= round(2 * k)
            self._shadow(gc, x, y, w, h, r, 1.35 if state == 1 else 0.6 if state == 2 else 1.0)
            # frosted white body, slightly translucent so the backdrop tints it
            fill = {0: 205, 1: 238, 2: 250}[state]
            gc.SetBrush(wx.Brush(wx.Colour(255, 255, 255, fill)))
            edge = wx.Colour(59, 130, 217, 150) if (state or focus) else wx.Colour(210, 218, 230, 200)
            gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(edge).Width(1)))
            self._rr(gc, x, y, w, h, r)
            # sheen: bright top-left fading out (158 degrees)
            gc.SetPen(wx.TRANSPARENT_PEN)
            gc.SetBrush(gc.CreateLinearGradientBrush(x, y, x + w * 0.38, y + h * 0.95,
                                                     wx.Colour(255, 255, 255, 215), wx.Colour(255, 255, 255, 0)))
            self._rr(gc, x + 1, y + 1, w - 2, h - 2, max(1, r - 1))
            # refraction: bright inner top edge, faint inner sides
            gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(wx.Colour(255, 255, 255, 240)).Width(1)))
            gc.StrokeLine(x + r, y + 1, x + w - r, y + 1)
            return y

        def _badge(self, gc, x, y, s, style):
            k = self.k
            x, y, s = round(x), round(y), round(s)
            r = s * 0.3
            gc.SetPen(wx.TRANSPARENT_PEN)
            if style == "brand":
                c1, c2, pen = wx.Colour("#FBE07A"), wx.Colour("#E8B425"), wx.Colour(255, 255, 255, 150)
                glow = (201, 154, 30)
            elif style == "accent":
                c1, c2, pen = wx.Colour(C["accent"]), wx.Colour(C["accent_end"]), wx.Colour(255, 255, 255, 70)
                glow = (23, 65, 143)
            else:
                c1, c2, pen = wx.Colour("#EEF3FC"), wx.Colour(C["accent_bg"]), wx.Colour(C["line2"])
                glow = None
            if glow:
                for dy, a in ((5, 26), (2.5, 40)):
                    gc.SetBrush(wx.Brush(wx.Colour(glow[0], glow[1], glow[2], a)))
                    self._rr(gc, x + 1 * k, y + dy * k, s - 2 * k, s, r)
            gc.SetBrush(gc.CreateLinearGradientBrush(x, y, x + s * 0.4, y + s, c1, c2))
            gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(pen).Width(1)))
            self._rr(gc, x, y, s, s, r)
            if style != "soft":
                gc.SetPen(wx.TRANSPARENT_PEN)
                gc.SetBrush(gc.CreateLinearGradientBrush(x, y, x, y + s * 0.55, wx.Colour(255, 255, 255, 80),
                                                         wx.Colour(255, 255, 255, 0)))
                self._rr(gc, x + 1.5 * k, y + 1.2 * k, s - 3 * k, s * 0.5, r * 0.8)

        def _glyph(self, gc, name, x, y, s, colour):
            """Icon strokes of a whole number of pixels, centred on the pixel
            grid, so the edges stay sharp (fractional widths and positions
            were smeared over two pixels by the anti-aliasing)."""
            import math
            parts = _GLYPHS.get(name)
            if not parts:
                return
            x, y = round(x), round(y)
            f = s / 24.0
            lw = max(2, int(round(1.7 * f)))

            def snap(v, w):
                return math.floor(v) + 0.5 if w % 2 else float(round(v))
            pen = gc.CreatePen(wx.GraphicsPenInfo(colour).Width(lw).Cap(wx.CAP_ROUND).Join(wx.JOIN_ROUND))
            gc.SetPen(pen)
            gc.SetBrush(wx.TRANSPARENT_BRUSH)
            P = lambda p, w=lw: (snap(x + p[0] * f, w), snap(y + p[1] * f, w))
            for kind, data in parts:
                if kind == "line":
                    gc.StrokeLines([wx.Point2D(*P(p)) for p in data])
                elif kind == "curve":
                    path = gc.CreatePath()
                    path.MoveToPoint(*P(data[0]))
                    for c1, c2, end in data[1:]:
                        path.AddCurveToPoint(*(P(c1) + P(c2) + P(end)))
                    gc.StrokePath(path)
                elif kind == "dots":
                    gc.SetBrush(wx.Brush(colour))
                    gc.SetPen(wx.TRANSPARENT_PEN)
                    for p in data:
                        cx, cy = P(p)
                        gc.DrawEllipse(cx - 1.4 * f, cy - 1.4 * f, 2.8 * f, 2.8 * f)
                    gc.SetPen(pen)
                    gc.SetBrush(wx.TRANSPARENT_BRUSH)
                elif kind == "bars":
                    bw = max(2, int(round(2.0 * f)))  # square ends, a clear gap between bars
                    bp = gc.CreatePen(wx.GraphicsPenInfo(colour).Width(bw).Cap(wx.CAP_BUTT))
                    gc.SetPen(bp)
                    for bx, top in data:
                        gc.StrokeLine(*(P((bx, 20), bw) + P((bx, top), bw)))
                    gc.SetPen(pen)
                elif kind == "rect":
                    rx, ry, rw, rh = data
                    gc.DrawRoundedRectangle(x + rx * f, y + ry * f, rw * f, rh * f, 2 * f)
                elif kind == "ellipse":
                    ex, ey, ew, eh = data
                    gc.DrawEllipse(x + ex * f, y + ey * f, ew * f, eh * f)

        def _logo_bitmap(self, h):
            if self._logo is None:
                return None
            if h not in self._logo_bmp:
                img = self._logo.Copy()
                w = int(img.GetWidth() * h / float(img.GetHeight()))
                self._logo_bmp[h] = wx.Bitmap(img.Scale(w, h, wx.IMAGE_QUALITY_HIGH))
            return self._logo_bmp[h]

        def _spaced(self, gc, text, x, y, spacing):
            for ch in text:
                gc.DrawText(ch, x, y)
                x += gc.GetTextExtent(ch)[0] + spacing
            return x

        def _wrap(self, gc, text, maxw, maxlines):
            words, lines, cur = text.split(), [], ""
            for wd in words:
                t = (cur + " " + wd).strip()
                if gc.GetTextExtent(t)[0] <= maxw or not cur:
                    cur = t
                else:
                    lines.append(cur)
                    cur = wd
            if cur:
                lines.append(cur)
            if len(lines) > maxlines:
                lines = lines[:maxlines]
                lines[-1] = lines[-1] + " \u2026"
            return lines

        # ---------------------------------------------------------- painting
        def _backdrop(self, gc, W, H):
            """Plain pale gradient in the palette, with a faint charge state
            series (UniDec's own subject) along the bottom edge."""
            k = self.k
            gc.SetPen(wx.TRANSPARENT_PEN)
            gc.SetBrush(gc.CreateLinearGradientBrush(0, 0, 0, H, wx.Colour("#F5F7FA"), wx.Colour("#EDF1F6")))
            gc.DrawRectangle(0, 0, W, H)
            import math
            base = H - 14 * k
            x0, x1 = W * 0.66, W * 0.97
            zmin, zmax, z0 = 9, 26, 15.5
            gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(wx.Colour(23, 65, 143, 20)).Width(1)))
            gc.StrokeLine(x0 - 30 * k, base, x1 + 10 * k, base)
            pen = gc.CreatePen(wx.GraphicsPenInfo(wx.Colour(23, 65, 143, 30)).Width(2.2 * k).Cap(wx.CAP_ROUND))
            gc.SetPen(pen)
            for z in range(zmin, zmax + 1):
                f = (1.0 / z - 1.0 / zmax) / (1.0 / zmin - 1.0 / zmax)  # m/z = M/z spacing
                x = x0 + (x1 - x0) * f
                h = H * 0.13 * math.exp(-((z - z0) ** 2) / (2 * 3.2 ** 2))
                if h > 2 * k:
                    gc.StrokeLine(x, base, x, base - h)

        def _on_paint(self, e):
            dc = wx.AutoBufferedPaintDC(self)
            gc = _Crisp(wx.GraphicsContext.Create(dc), dc)
            try:
                self._paint(gc)
            finally:
                gc.finish()

        def _paint(self, gc):
            k = self.k
            W, H = self.GetClientSize()
            self._backdrop(gc, W, H)

            M = 40 * k
            # header: app badge, name, version
            bs = round(64 * k)
            bx, by = round(M), round(M)
            self._badge(gc, bx, by, bs, "accent")
            m = round(bs * 0.18)
            self._glyph(gc, "app", bx + m, by + m, bs - 2 * m, wx.Colour(255, 255, 255))
            tx = bx + bs + 20 * k
            gc.SetFont(ui_font(24, 800), wx.Colour(C["text"]))
            tw, th = gc.GetTextExtent(APP_NAME)
            gc.DrawText(APP_NAME, tx, by - 2 * k)
            gc.SetFont(ui_font(11, 600), wx.Colour(C["faint"]))
            ver = APP_VERSION
            vw, vh = gc.GetTextExtent(ver)
            px, py = tx + tw + 12 * k, by - 2 * k + th * 0.5 - vh / 2.0 - 3 * k
            gc.SetBrush(wx.Brush(wx.Colour(23, 65, 143, 16)))
            gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(wx.Colour(23, 65, 143, 30)).Width(1)))
            self._rr(gc, px, py, vw + 16 * k, vh + 6 * k, (vh + 6 * k) / 2.0)
            gc.SetFont(ui_font(11, 600), wx.Colour(C["accent_text"]))
            gc.DrawText(ver, px + 8 * k, py + 3 * k)

            # three workspace tiles
            self._rects = []
            y = by + bs + 40 * k
            gap = 20 * k
            n = len(self.tools)
            cols = max(1, min(3, n))
            tw_ = (W - 2 * M - (cols - 1) * gap) / cols
            # tile height from its content (title, up to 3 lines of text)
            gc.SetFont(ui_font(14, 700), wx.Colour(C["text"]))
            title_h = max(gc.GetTextExtent(t[0])[1] for t in self.tools) if self.tools else 20 * k
            gc.SetFont(ui_font(9.5, 400), wx.Colour(C["muted"]))
            lh = gc.GetTextExtent("Xg")[1]
            nl = max([len(self._wrap(gc, t[1].replace("\n", " "), tw_ - 44 * k, 3)) for t in self.tools] or [1])
            th_ = max(170 * k, 22 * k + 52 * k + 16 * k + title_h + 6 * k + nl * lh * 1.18 + 20 * k)
            need = y + ((n + cols - 1) // cols) * (th_ + gap) + 60 * k
            if need > H + 1 and not getattr(self, "_grown", False):
                self._grown = True
                wx.CallAfter(self._grow, int(need))
            for i in range(n):
                cx = M + (i % cols) * (tw_ + gap)
                cy = y + (i // cols) * (th_ + gap)
                self._tile(gc, i, cx, cy, tw_, th_)
            y += ((n + cols - 1) // cols) * (th_ + gap)

            # "More tools" link (hidden: the classic tools are in the right click menu)
            if self.more and self.show_more:
                idx = len(self.tools)
                hov = self._hover == idx or self._focus == idx
                gc.SetFont(ui_font(9.5, 600), wx.Colour(C["accent_text"] if hov else C["muted"]))
                label = "More tools"
                lw, lh = gc.GetTextExtent(label)
                lx, ly = M + 2 * k, y + 2 * k
                gc.DrawText(label, lx, ly)
                cxv, cyv, sv = lx + lw + 9 * k, ly + lh / 2.0, 3.5 * k
                gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(wx.Colour(C["accent_text"] if hov else C["muted"]))
                                       .Width(1.4 * k).Cap(wx.CAP_ROUND)))
                gc.StrokeLines([wx.Point2D(cxv - sv, cyv - sv / 2.0), wx.Point2D(cxv, cyv + sv / 2.0),
                                wx.Point2D(cxv + sv, cyv - sv / 2.0)])
                if hov:
                    gc.StrokeLine(lx, ly + lh + 1 * k, lx + lw, ly + lh + 1 * k)
                self._rects.append((lx - 4 * k, ly - 4 * k, lw + 26 * k, lh + 8 * k, idx))

            gc.SetFont(ui_font(8.5, 500), wx.Colour(C["faint"]))
            cite = "MSpektra"
            lines = self._wrap(gc, cite, W * 0.60 - M, 2)
            ch = gc.GetTextExtent("Xg")[1]
            for j, ln in enumerate(lines):
                gc.DrawText(ln, M + 2 * k, H - M * 0.6 - ch * (len(lines) - j) * 1.1)

        def _tile(self, gc, i, x, y, w, h):
            k = self.k
            title, desc, _id, group = self.tools[i]
            self._rects.append((x, y, w, h, i))
            state = 2 if i == self._down else 1 if i == self._hover else 0
            y = self._glass(gc, x, y, w, h, 22 * k, state, focus=(i == self._focus))
            s = round(52 * k)
            px, py = round(x + 22 * k), round(y + 22 * k)
            self._badge(gc, px, py, s, "brand" if title == "Deconvolute" else "accent")
            m = round(s * 0.18)
            self._glyph(gc, title, px + m, py + m, s - 2 * m,
                        wx.Colour("#1F2A55") if title == "Deconvolute" else wx.Colour(255, 255, 255))
            ty = py + s + 16 * k
            gc.SetFont(ui_font(14, 700), wx.Colour(C["text"]))
            ttw, tth = gc.GetTextExtent(title)
            gc.DrawText(title, px, ty)
            if state == 1:
                ax, ay = x + w - 28 * k, py + s / 2.0
                gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(wx.Colour(C["accent"])).Width(2.0 * k)
                                       .Cap(wx.CAP_ROUND).Join(wx.JOIN_ROUND)))
                gc.StrokeLines([wx.Point2D(ax - 5 * k, ay - 6 * k), wx.Point2D(ax + 1 * k, ay),
                                wx.Point2D(ax - 5 * k, ay + 6 * k)])
            gc.SetFont(ui_font(9.5, 400), wx.Colour(C["muted"]))
            maxw = w - 44 * k
            lines = self._wrap(gc, desc.replace("\n", " "), maxw, 3)
            lh = gc.GetTextExtent("Xg")[1]
            yy = ty + tth + 6 * k
            for ln in lines:  # (line spacing as used for the tile height)
                if gc.GetTextExtent(ln)[0] > maxw:
                    while ln and gc.GetTextExtent(ln + "\u2026")[0] > maxw:
                        ln = ln[:-1]
                    ln = ln.rstrip() + "\u2026"
                gc.DrawText(ln, px, yy)
                yy += lh * 1.18

    return GlassLauncher


def restyle_launcher(frame):
    import wx
    old = [c for c in frame.GetChildren() if isinstance(c, wx.Panel)]
    if not old:
        return
    oldpanel = old[0]
    classic = []
    unidec_id = None
    for c in oldpanel.GetChildren():
        if isinstance(c, wx.Button):
            lines = [ln.strip() for ln in c.GetLabel().splitlines()]
            title = lines[0] if lines else c.GetLabel()
            desc = " ".join(ln for ln in lines[1:] if ln)
            if title == "UniDec":
                unidec_id = c.GetId()
            classic.append((title, desc, c.GetId()))
    if not classic:
        return
    order = {t: i for i, t in enumerate(_LAUNCH_ORDER)}
    classic.sort(key=lambda b: order.get(b[0], 99))
    more = [c for c in classic if c[0] != "UniDec"]

    def open_lcms(e=None):
        import unilcms
        unilcms.open_window()

    def open_hrms(e=None):
        import hrms
        hrms.open_window()

    tools = []
    for title, desc, key in WORKSPACES:
        if key == "deconv":
            if unidec_id is None:
                continue
            tools.append((title, desc, unidec_id, "main"))
            continue
        tid = wx.NewIdRef()
        frame.Bind(wx.EVT_BUTTON, open_lcms if key == "lcms" else open_hrms, id=tid)
        setattr(frame, "_udp_%s_id" % key, tid)
        tools.append((title, desc, int(tid), "main"))
    version = str(getattr(getattr(frame, "eng", None), "version", "") or "")
    canvas = _glass_class()(frame, tools, getattr(frame, "imagepath", None), version, more=more)
    oldpanel.Hide()
    fs = wx.BoxSizer(wx.VERTICAL)
    fs.Add(canvas, 1, wx.EXPAND)
    frame.SetSizer(fs, True)
    frame.SetBackgroundColour(wx.Colour(C["bg"]))
    frame.SetTitle(APP_NAME)
    set_app_icon(frame)
    fs.Fit(frame)
    frame.Layout()
    frame.Centre()
    canvas.SetFocus()


# --------------------------------------------------------------------------
# fast start: the start screen needs only wx; the libraries of the windows
# are loaded in the background while it is shown
# --------------------------------------------------------------------------
UNIDEC_VERSION = "8.2.1"
CLASSIC_TOOLS = [  # (title, description, button number in UniDec's Launcher)
    ("Batch MS analysis", "Batch process and visualize MS spectra", 4),
    ("Charge-detection deconvolution", "Deconvolution of charge detection MS", 9),
    ("CD-MS chromatography", "Deconvolution of CD-MS chromatograms", 11),
    ("IsoDec", "Deconvolution of isotope patterns", 12),
    ("LC-MS deconvolution", "Deconvolution of chromatograms", 8),
    ("Data Collector", "Visualize multiple spectra, extract trends and Kd's", 2),
    ("UltraMeta Data Collector", "Visualize multiple HDF5 data sets, fit trends", 7),
    ("Import Wizard", "Batch convert Waters Raw to Txt", 3),
    ("HDF5 Import Wizard", "Import data into HDF5 for batch analysis", 6),
    ("Batch processing pipeline", "Batch processing workflow", 10),
    ("Analysis API console", "Script the analysis with a console", 5),
]
# Libraries without windows only (safe to import outside the main thread).
# Only what LCMS Analysis and HRMS Analysis import themselves: NumPy, the parts
# of Matplotlib a plot needs (the font list first: on the first start on a
# computer it is built here, which takes a few seconds), and the SciPy
# modules of the LCMS peak integration. The deconvolution methods run in the
# worker process (deconv_worker.py), and the Deconvolute window and the
# classic tools open in their own process, so UniDec's libraries (numba,
# pandas, h5py ... about 10 s on a cold start) are never loaded into the
# process of the start screen and the analysis windows. Importing is work on
# Python's global lock: a second thread does not make it faster, it only
# slows the window being built or drawn on the main thread meanwhile. A
# click on a tile waits for the module being imported at that moment, so
# the SciPy packages scipy.signal pulls in are listed one by one (short
# waits) and the rarely needed ones come last.
PRELOAD_POSTRUN = ("numpy", "matplotlib", "matplotlib.font_manager", "matplotlib.colors", "matplotlib.transforms",
                   "matplotlib.ticker", "matplotlib.text", "matplotlib.patches", "matplotlib.collections",
                   "matplotlib.backend_bases", "matplotlib.figure", "matplotlib.backends.backend_agg",
                   "scipy.special", "scipy.linalg", "scipy.fft", "scipy.sparse", "scipy.ndimage", "scipy.spatial",
                   "scipy.optimize", "scipy.interpolate", "scipy.stats", "scipy.signal")
# the libraries of the Deconvolute window, loaded in the background only
# when it is opened inside this process ("deconvolute_own_process": false)
PRELOAD_DECONV = ("scipy.interpolate", "scipy.stats", "scipy.sparse", "scipy.special", "scipy.fft", "scipy.linalg",
                  "scipy.optimize", "scipy.constants", "llvmlite.binding", "numba", "pandas", "dateutil.parser",
                  "PIL.Image", "PIL.PngImagePlugin", "pyparsing", "fontTools", "mpl_toolkits.mplot3d",
                  "mpl_toolkits.axes_grid1", "matplotlib.widgets", "matplotlib.tri", "matplotlib.style",
                  "matplotlib.offsetbox", "matplotlib.legend", "pubsub.pub", "h5py", "regex", "lxml.etree",
                  "lxml.html", "pymzml", "pyteomics.mzxml", "pyteomics.mass", "cffi", "jinja2", "mpld3", "networkx",
                  "natsort", "wheezy.template", "pyimzml.ImzMLParser")
_PRE = {"pause_until": 0.0, "started": False, "done": False, "hold": 0, "lock": None}


def _pre_lock():
    if _PRE["lock"] is None:
        import threading
        _PRE["lock"] = threading.Lock()
    return _PRE["lock"]


def pause_preload(seconds=5.0):
    """The background loading waits a few seconds (a window was built)."""
    import time
    _PRE["pause_until"] = max(_PRE["pause_until"], time.time() + seconds)


class hold_preload(object):
    """While a window is built on the main thread, the background loading
    stops; entering waits until the module it is loading is finished, so
    the two threads never import at the same time. The hold flag is set
    before waiting for the lock: the loading thread checks it before taking
    the lock for the next module, so entering waits for one module only.
    (Taking the lock first, to set the flag, made a click on a tile wait
    for the whole list: the loading thread releases and takes the lock
    again at once between two modules, and a thread already waiting for a
    lock rarely gets it in between.)
    wait=False: no wait for the module being loaded (a file dialog of
    Windows: it opens at once, the loading stops after that module)."""

    def __init__(self, wait=True):
        self.wait = wait

    def __enter__(self):
        _PRE["hold"] += 1  # only the main thread changes it
        if self.wait:
            lock = _pre_lock()
            lock.acquire()
            lock.release()
        return self

    def __exit__(self, *a):
        _PRE["hold"] = max(0, _PRE["hold"] - 1)
        pause_preload(3)
        return False


def start_preload(delay=0.0, deconv=False):
    """Loads the libraries of LCMS/HRMS Analysis in a background thread
    (once); with deconv=True also those of the Deconvolute window, for a
    start screen that opens it inside its own process."""
    if _PRE["started"]:
        return
    _PRE["started"] = True
    import threading
    import time
    groups = [("postrun", PRELOAD_POSTRUN)]
    if deconv:
        groups.append(("deconvolution", PRELOAD_DECONV))

    def work():
        import importlib
        import sys
        time.sleep(delay)
        lock = _pre_lock()
        for label, names in groups:
            t0 = time.perf_counter()
            for name in names:
                while True:
                    while _PRE["hold"] > 0 or time.time() < _PRE["pause_until"]:
                        time.sleep(0.1)
                    with lock:
                        if _PRE["hold"] > 0:
                            continue
                        if name not in sys.modules:
                            try:
                                importlib.import_module(name)
                            except Exception:
                                pass
                        break
            _log("%s libraries loaded in the background in %.1f s" % (label, time.perf_counter() - t0))
        _PRE["done"] = True

    threading.Thread(target=work, name="preload", daemon=True).start()


def _process_has_window(pid):
    """True when the process shows a visible top level window (Windows)."""
    import ctypes
    from ctypes import wintypes
    user32 = ctypes.windll.user32
    found = []

    @ctypes.WINFUNCTYPE(ctypes.c_bool, wintypes.HWND, wintypes.LPARAM)
    def visit(hwnd, lparam):
        owner = wintypes.DWORD()
        user32.GetWindowThreadProcessId(hwnd, ctypes.byref(owner))
        if owner.value == pid and user32.IsWindowVisible(hwnd):
            found.append(hwnd)
            return False
        return True

    user32.EnumWindows(visit, 0)
    return bool(found)


def start_screen(ensure_unidec=None, root=None, spawn=None):
    """The start screen (LCMS Analysis, HRMS Analysis, Deconvolute) in its own
    wx.App; returns when every window is closed. spawn(args) starts MS
    Analysis again in a new process (launch_unidec.spawn): the Deconvolute
    window and the classic tools then open in their own process and UniDec's
    libraries never load into this one. Without it (or with
    "deconvolute_own_process": false in the settings) they open inside this
    process, as before, with their libraries loaded in the background."""
    import time
    import importlib
    import wx
    t0 = time.perf_counter()
    app = wx.App(False)
    app.SetAppName(APP_NAME)
    frame = wx.Frame(None, title=APP_NAME)
    frame.SetBackgroundColour(wx.Colour(C["bg"]))
    try:
        own_process = spawn is not None and _load().get("deconvolute_own_process", True) is not False
    except Exception:
        own_process = spawn is not None

    def need_unidec():
        if ensure_unidec is not None:
            ensure_unidec()
        else:
            install_unidec()

    # One window at a time: while fn builds a window the main thread is busy
    # and clicks queue up; they would open the same window again (or worse)
    # right after it. So a click is ignored while a window is being built
    # and for a moment afterwards (the queued ones), also after a failure.
    guard = {"building": False, "until": 0.0}
    loading = {}  # title -> process of a window that is opening in its own process

    def busy(fn):
        """fn builds the window; it may return a function to run afterwards,
        outside the busy cursor (UniDec's windows run their own event loop)."""
        def run(e=None):
            if guard["building"] or time.monotonic() < guard["until"]:
                _log("click ignored: a window is being opened")
                return
            guard["building"] = True
            then = None
            try:
                with wx.BusyCursor(), hold_preload():
                    then = fn()
            except Exception:
                import traceback
                tb = traceback.format_exc()
                print(tb)
                wx.MessageBox("The window could not be opened:\n\n" + tb[-1200:], APP_NAME, wx.ICON_ERROR)
            finally:
                guard["building"] = False
                guard["until"] = time.monotonic() + 0.6
            if callable(then):
                try:
                    then()
                except Exception:
                    import traceback
                    tb = traceback.format_exc()
                    print(tb)
                    wx.MessageBox("The window could not be opened:\n\n" + tb[-1200:], APP_NAME, wx.ICON_ERROR)
        return run

    def open_lcms():
        import unilcms
        unilcms.open_window()

    def open_hrms():
        import hrms
        hrms.open_window()

    def launch(args, title, fallback):
        """The window in a new process. The busy cursor stays until that
        process shows a window; if it ends without one, fallback (the
        in-process way) runs instead. Returns False when no process could
        be started, so the caller falls back at once. A second click on the
        same tile while its process is still loading is ignored."""
        if title in loading:
            _log("%s is being opened already" % title)
            return True
        try:
            proc = spawn(args)
        except Exception as e:
            _log("%s not started as its own process: %s" % (title, e))
            return False
        if proc is None:
            return False
        t1 = time.perf_counter()
        _log("%s started as its own process (pid %d)" % (title, proc.pid))
        loading[title] = proc
        wx.BeginBusyCursor()
        try:
            frame.SetTitle("%s: opening %s" % (APP_NAME, title))
        except Exception:
            pass
        state = {"done": False}

        def finish():
            if not state["done"]:
                state["done"] = True
                loading.pop(title, None)
                try:
                    wx.EndBusyCursor()
                    frame.SetTitle(APP_NAME)
                except Exception:
                    pass

        def poll():
            if state["done"]:
                return
            code = proc.poll()
            if code is not None:
                finish()
                if code not in (0, 3):  # 3: its libraries failed to load, reported there already
                    _log("%s process ended with code %s before showing a window; opening it here" % (title, code))
                    busy(fallback)()
                return
            try:
                shown = _process_has_window(proc.pid)
            except Exception:
                shown = True
            if shown:
                _log("%s window shown after %.1f s" % (title, time.perf_counter() - t1))
                finish()
            elif time.perf_counter() - t1 > 120:
                finish()
            else:
                wx.CallLater(250, poll)
        wx.CallLater(250, poll)
        return True

    def deconv_here():
        need_unidec()
        from unidec.GUniDec import UniDecApp
        print("Launching UniDec")
        app2 = UniDecApp()
        return app2.start

    def open_deconv():
        if own_process and launch(["--deconv"], "Deconvolute", deconv_here):
            return None
        return deconv_here()

    def classic(n):
        def here():
            need_unidec()
            L = importlib.import_module("unidec.Launcher")
            return lambda: getattr(L.Lview, "button%d" % n)(None)

        title = next((t for t, d, k in CLASSIC_TOOLS if k == n), "Classic tool %d" % n)

        def run():
            if own_process and launch(["--classic", str(n)], title, here):
                return None
            return here()
        return run

    handlers = {"lcms": open_lcms, "hrms": open_hrms, "deconv": open_deconv}
    tools = []
    for title, desc, key in WORKSPACES:
        tid = wx.NewIdRef()
        frame.Bind(wx.EVT_BUTTON, busy(handlers[key]), id=tid)
        tools.append((title, desc, int(tid), "main"))
    more = []
    for title, desc, n in CLASSIC_TOOLS:
        tid = wx.NewIdRef()
        frame.Bind(wx.EVT_BUTTON, busy(classic(n)), id=tid)
        more.append((title, desc, int(tid)))
    canvas = _glass_class()(frame, tools, None, UNIDEC_VERSION, more=more)
    fs = wx.BoxSizer(wx.VERTICAL)
    fs.Add(canvas, 1, wx.EXPAND)
    frame.SetSizer(fs)
    set_app_icon(frame)
    fs.Fit(frame)
    frame.Centre()
    try:  # before Show: the taskbar button has the icon of this MSpektra.exe from the start
        import unidec_ui_addons
        unidec_ui_addons.install_relaunch(root, frames=[frame])
    except Exception as e:
        _log("taskbar:", e)
    frame.Show()
    canvas.SetFocus()
    _log("start screen shown in %.1f s" % (time.perf_counter() - t0))

    def after_show():
        start_preload(deconv=not own_process)
    wx.CallLater(200, after_show)
    app.MainLoop()
    return 0


def set_app_icon(win):
    import wx
    ico = os.path.join(os.path.dirname(os.path.abspath(__file__)), "msanalysis.ico")
    if os.path.isfile(ico):
        try:
            bundle = wx.IconBundle(ico, wx.BITMAP_TYPE_ICO)
            win.SetIcons(bundle)
        except Exception:
            try:
                win.SetIcon(wx.Icon(ico, wx.BITMAP_TYPE_ICO))
            except Exception:
                pass


def _install_launcher():
    import importlib
    try:
        L = importlib.import_module("unidec.Launcher")
    except Exception as e:
        _log("launcher:", e)
        return
    orig = L.Lview.__init__

    def __init__(self, parent):
        orig(self, parent)
        try:
            restyle_launcher(self)
        except Exception as e:
            _log("launcher restyle failed:", e)

    L.Lview.__init__ = __init__


# --------------------------------------------------------------------------
# main entry points
# --------------------------------------------------------------------------
def apply_main_theme(frame):
    """Theme a UniDec main window just before it is shown."""
    import wx
    if getattr(frame, "_udp_themed", False):
        return
    frame._udp_themed = True
    if _ST["dpi_aware"]:
        scale_fixed_sizes(frame)
        try:
            ctl = getattr(frame, "controls", None)
            s = ctl.GetSizer() if ctl is not None else None
            f = _scale(frame)
            if s is not None and f > 1.01:
                imflag = getattr(frame.config, "imflag", 0) or 0
                s.SetMinSize(wx.Size(int((250 + 10 * imflag) * f), 0))
        except Exception:
            pass
    theme_tree(frame)
    try:
        if frame.GetTitle().strip() == "UniDec":
            frame.SetTitle("Deconvolute")
        set_app_icon(frame)
    except Exception:
        pass
    # UniDec and MetaUniDec bind their buttons on the main window, so the new
    # toolbar can reuse the button IDs; other windows keep their buttons
    ctl_mod = type(getattr(frame, "controls", None)).__module__
    toolbar_ok = ctl_mod in ("unidec.modules.gui_elements.ud_controls",
                             "unidec.metaunidec.gui_elements.ud_cont_meta") and \
        all(getattr(frame.controls, a, None) is not None
            for a in ("openbutton", "procbutton", "udbutton", "ppbutton", "autobutton"))
    try:
        from unidec.modules.mainwindow import Mainwindow
        is_unidec = isinstance(frame, Mainwindow)
    except Exception:
        is_unidec = False
    steps = []
    if toolbar_ok:
        steps += [_build_toolbar, _restyle_controls]
    if is_unidec:
        steps.append(_make_cards)
    for step in steps:
        try:
            step(frame)
        except Exception as e:
            _log(step.__name__, "failed:", e)
    relayout_foldpanels(frame)
    frame.Layout()


def _install_generic_show():
    """Every frame and dialog, when it is shown: inside the work area of its
    display (fit_to_screen), again after it moved to a display of another
    scale; on scaled displays the fixed sizes of UniDec's own windows are
    scaled first (main windows are themed and placed by apply_main_theme and
    _install_launch)."""
    import wx
    native = tuple(c for c in (getattr(wx, n, None) for n in (
        "MessageDialog", "FileDialog", "DirDialog", "ColourDialog", "FontDialog", "PrintDialog",
        "PageSetupDialog", "ProgressDialog", "FindReplaceDialog")) if isinstance(c, type))

    def refit(win):
        try:
            if win and win.IsShown():
                fit_to_screen(win, grow=False)
        except RuntimeError:
            pass

    def fix(win):
        if isinstance(win, native):
            return  # drawn by Windows, placed by it
        first = not getattr(win, "_udp_fitted", False)
        themed = getattr(win, "_udp_themed", False)
        if first and not themed and not getattr(win, "_udp_scaled", False) and _from_unidec(win) and \
                _ST["dpi_aware"] and _ST["system_scale"] > 1.01:
            scale_fixed_sizes(win)
        fit_to_screen(win, grow=first and not themed)
        if first:
            win._udp_fitted = True
            if win.GetSizer() is not None:
                win.Layout()
            try:  # moved to a display of another scale: Windows resizes it, keep it inside
                win.Bind(wx.EVT_DPI_CHANGED, lambda e, w=win: (e.Skip(), wx.CallAfter(refit, w)))
            except Exception:
                pass

    for cls in (wx.Frame, wx.Dialog):
        orig = cls.Show

        def make(orig):
            def Show(self, show=True):
                if show:
                    try:
                        fix(self)
                    except Exception as e:
                        _log("window fit:", e)
                return orig(self, show)
            return Show

        cls.Show = make(orig)
    orig_modal = wx.Dialog.ShowModal

    def ShowModal(self):
        try:
            fix(self)
        except Exception as e:
            _log("window fit:", e)
        return orig_modal(self)

    wx.Dialog.ShowModal = ShowModal


def _install_plot_size():
    """Plot windows start at their 96-dpi pixel size on scaled displays."""
    if not _ST["dpi_aware"] or _ST["system_scale"] <= 1.01:
        return
    from unidec.modules.plotting import PlottingWindow as pw
    base = pw.PlottingWindowBase
    orig = base.__init__
    f = _ST["system_scale"]

    def __init__(self, *args, **kwargs):
        try:
            fs = kwargs.get("figsize") or (6. * 0.9, 5. * 0.9)
            kwargs["figsize"] = (float(fs[0]) * f, float(fs[1]) * f)
        except Exception:
            pass
        return orig(self, *args, **kwargs)

    base.__init__ = __init__


def configure(settings_path=None):
    """The light part (only wx): settings file, bundled fonts, size fixes on
    scaled displays. Enough for the start screen, LCMS and HRMS Analysis."""
    if _ST.get("configured"):
        return
    _ST["configured"] = True
    _ST["settings_path"] = settings_path
    try:
        register_fonts(os.path.join(os.path.dirname(os.path.abspath(__file__)), "fonts"))
    except Exception as e:
        _log("fonts:", e)
    for name, func in (("display size", _patch_display_size), ("other windows", _install_generic_show)):
        try:
            func()
        except Exception as e:
            _log(name, "NOT enabled:", e)


def install_unidec():
    """The UniDec windows (imports UniDec's interface modules, so it is done
    only when the Deconvolute window or a classic tool is opened)."""
    if _ST["installed"]:
        return
    _ST["installed"] = True
    steps = (("section headers", _install_captions), ("plot cards", _install_resize_plots),
             ("plot size", _install_plot_size), ("peak list", _install_peaklist),
             ("window set-up", _install_launch), ("last folder", _install_last_folder))
    for name, func in steps:
        try:
            func()
        except Exception as e:
            _log(name, "NOT enabled:", e)
    _log("modern interface enabled (DPI aware: %s, scale %.2f)" % (_ST["dpi_aware"], _ST["system_scale"]))


def install(settings_path=None):
    configure(settings_path)
    install_unidec()

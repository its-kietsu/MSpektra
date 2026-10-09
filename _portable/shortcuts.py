"""Keyboard shortcuts of the Postrun windows of MS Analysis: one registry.

Every shortcut is an entry (keys, action, where, scope, handler). The registry builds
the accelerator table of each window (Ctrl combinations that never collide with
typing), dispatches the other keys from the window's EVT_CHAR_HOOK (it comes before
the focused control gets the key), adds the keys to the tooltips of the tool bar and
to menu entries, and fills the F1 window (a searchable list of every shortcut).

Scopes
* window: anywhere in the window, also while the cursor is in a field (Ctrl+O, Ctrl+W,
  Ctrl+R, Ctrl+Tab, Ctrl+1/2/3, Ctrl+E, Ctrl+D, Ctrl+I, F1, F9).
* undo: Ctrl+Z, Ctrl+Y, Ctrl+Shift+Z undo or redo a step of the window, except in a
  text field, where they undo the typing in that field (as Windows does).
* view: keys of the view shown (Home, + and -, arrows, Delete, Esc) when the focus is
  on a plot or the view itself, never in a field or a list.
* plot: single letters for the tools, only while a plot (or the tool bar) has the
  focus (click a plot once); typing a letter in a field never switches a tool.
* native and mouse entries are handled by their controls (lists, dialogs, plots) and
  are listed for F1 only.
"""
import sys
import wx

WINDOWS = {"lcms": "LCMS Postrun", "hrms": "HRMS Postrun"}
VIEWS = {"ms": "Mass spectrometry view", "pda": "PDA view", "compare": "Compare view"}

CTRL, SHIFT, ALT = wx.MOD_CONTROL, wx.MOD_SHIFT, wx.MOD_ALT


class Shortcut(object):
    def __init__(self, keys, action, scope, handler=None, combos=(), views=None, windows=("lcms", "hrms"),
                 accel=False, tool=None, where=None, menu=None):
        self.keys, self.action, self.scope, self.handler = keys, action, scope, handler
        self.combos = list(combos)  # (modifiers, key code)
        self.views = views  # None: every view; else a set of "ms", "pda", "compare"
        self.windows = windows
        self.accel = accel  # in the accelerator table of the window (else dispatched from EVT_CHAR_HOOK)
        self.tool = tool  # the tool of the tool bar it chooses (its tooltip shows the key)
        self.where = where  # text for F1 (default: from scope and views)
        self.menu = menu  # menu entries (label start) that show this key

    def where_text(self):
        if self.where:
            return self.where
        if self.views is None:
            v = "every view"
        else:
            v = ", ".join(VIEWS[k] for k in ("ms", "pda", "compare") if k in self.views)
        if self.scope == "plot":
            return v + " (a plot clicked)"
        if self.scope == "view":
            return v + " (not in a field)"
        if self.scope == "undo":
            return v + " (in a field: the typing there)"
        return v


# ------------------------------------------------------------------ context
class Ctx(object):
    """What a key acts on: the window, the file shown, the view shown, the plot with the focus."""

    def __init__(self, frame):
        self.frame = frame
        self.doc = getattr(frame, "active", None)
        self.focus = wx.Window.FindFocus()
        self.kind = focus_kind(self.focus)
        self.tab = view_of(frame)
        self.view = view_kind(self.tab)
        self.card = plot_of(self.focus)
        if self.card is not None and self.tab is not None:
            try:
                if self.card not in self.tab.all_cards() and not _inside(self.card, self.tab):
                    self.card = None
            except Exception:
                pass

    def last_card(self):
        """The plot with the focus, else the plot clicked last in this view."""
        if self.card is not None:
            return self.card
        c = self.frame.__dict__.get("_last_card")
        try:
            if c and self.tab is not None and _inside(c, self.tab) and c.IsShownOnScreen():
                return c
        except Exception:
            pass
        return None


def _inside(win, parent):
    w = win
    while w is not None:
        if w is parent:
            return True
        try:
            w = w.GetParent()
        except Exception:
            return False
    return False


def window_kind(frame):
    n = type(frame).__name__
    return "hrms" if n.startswith("HRMS") else "lcms"


def view_of(frame):
    """The view shown: the Compare view while it is shown, else the page of the file shown."""
    try:
        tab = frame.__dict__.get("_compare_tab")
        if tab is not None and tab.IsShown():
            return tab
        doc = frame.active
        i = doc.book.GetSelection()
        pg = doc.pages[i][1] if 0 <= i < len(doc.pages) else None
        if pg is not None and getattr(pg, "is_compare", False):
            return frame.compare_tab(create=False)
        return pg
    except Exception:
        return None


def view_kind(tab):
    if tab is None:
        return None
    names = [c.__name__ for c in type(tab).__mro__]
    if "CompareTab" in names:
        return "compare"
    if "PDATab" in names:
        return "pda"
    if "MSTab" in names:
        return "ms"
    return None


def plot_of(w):
    """The PlotCard of a focused window (the matplotlib canvas inside it)."""
    import unilcms as U
    for _ in range(4):
        if w is None:
            return None
        if isinstance(w, U.PlotCard):
            return w
        try:
            w = w.GetParent()
        except Exception:
            return None
    return None


def focus_kind(w):
    """text (typing goes there), list (its own keys), control (check box, button), plot, view, none."""
    if w is None:
        return "none"
    try:
        p = w.GetParent()
    except Exception:
        p = None
    textish = (wx.TextCtrl, wx.SpinCtrl, wx.ComboBox, wx.SearchCtrl)
    if hasattr(wx, "SpinCtrlDouble"):
        textish = textish + (wx.SpinCtrlDouble,)
    if isinstance(w, textish) or isinstance(p, textish):
        return "text"
    if isinstance(w, (wx.Choice, wx.ListBox, wx.CheckListBox, wx.ListCtrl, wx.TreeCtrl)):
        return "list"
    names = [c.__name__ for c in type(w).__mro__]
    if "FilesPanel" in names or "PeakTable" in names or "TableCard" in names:
        return "list"
    if "FigureCanvasWxAgg" in names or "FigureCanvasWx" in names:
        return "plot"
    if isinstance(w, (wx.CheckBox, wx.RadioButton, wx.Button, wx.ToggleButton)) or "ColourPickerCtrl" in names \
            or "FlatButton" in names:
        return "control"
    return "view"


# ------------------------------------------------------------------ handlers
def _undo(ctx):
    import undo
    undo.do_undo(ctx.frame)
    return True


def _redo(ctx):
    import undo
    undo.do_redo(ctx.frame)
    return True


def _open(ctx):
    ctx.frame.on_open(None)  # the raw data dialog; Open in the top bar also has the saved analyses


def _close(ctx):
    ctx.frame.close_doc(ctx.frame.active)


def _report(ctx):
    ctx.frame.make_report()


def _next(step):
    def h(ctx):
        ctx.frame.next_doc(step)
    return h


def _show_view(i):
    def h(ctx):
        doc = ctx.frame.active
        if doc is None or not getattr(doc, "pages", None) or i >= len(doc.pages):
            return False
        en = list(getattr(doc, "seg_enabled", None) or [True] * len(doc.pages))
        if i < len(en) and not en[i]:
            ctx.frame.status_msg("%s: this file has no data for it" % doc.pages[i][0])
            return True
        doc.show_page(doc.pages[i][1])
        return True
    return h


def _deconvolute(ctx):
    ctx.frame.deconvolute()


def _integrate(ctx):
    if ctx.view in ("ms", "pda"):
        ctx.tab.tool_action("auto")
        return True
    return False


def _side(ctx):
    ctx.frame.set_side(None)
    return True


def _help(ctx):
    show_window(ctx.frame)
    return True


def _export(ctx):
    """The data of the plot with the focus (or clicked last); else the main plot of the view."""
    tab, card = ctx.tab, ctx.last_card()
    if tab is None:
        return False
    try:
        if ctx.view == "compare":
            for pol, c in getattr(tab, "mzcards", {}).items():
                if card is c:
                    tab.export_spectrum(pol)
                    return True
            tab.on_export_traces()
            return True
        for v in tab.views:
            if card is v.card:
                tab.export_chrom(v)
                return True
        if ctx.view == "ms":
            for col in tab.cols:
                if card is col["spec_card"]:
                    tab.on_export_spectrum(col["e"])
                    return True
            tab.on_export_spectrum(None)
            return True
        if ctx.view == "pda":
            tab.on_export_uv()
            return True
    except Exception as ex:
        ctx.frame.status_msg("Export failed: %s" % ex)
        return True
    return False


def _full(ctx):
    if ctx.tab is None:
        return False
    ctx.tab.tool_action("full")
    if ctx.card is not None and ctx.card not in ctx.tab.all_cards():
        ctx.card.reset_view()  # (a plot of a deconvolution result)
    return True


def _zoom(f):
    def h(ctx):
        tab = ctx.tab
        if tab is None:
            return False
        cards = [ctx.card] if ctx.card is not None else [c for c in tab.all_cards()[:1]]
        for card in cards:
            if card is None or not card.full:
                continue
            x0, x1 = card.ax.get_xlim()
            c, h_ = 0.5 * (x0 + x1), 0.5 * (x1 - x0) * f
            a, b = c - h_, c + h_
            lo, hi = card.full
            a, b = max(a, lo), min(b, hi)
            if b <= a:
                continue
            card.ax.set_xlim(a, b)
            card.ylock = None
            if hasattr(card, "set_cross"):  # the PDA map: its crosshair follows
                card.set_cross(card.t, card.wl, redraw=False)
            else:
                card.autoscale_y()
            card.canvas.draw_idle()
            if card.on_view:
                card.on_view()
        return True
    return h


def _step(n, dj=0):
    def h(ctx):
        tab = ctx.tab
        if ctx.view == "pda":
            if ctx.card is not None and hasattr(ctx.card, "on_step"):
                return True  # the map stepped its crosshair itself (matplotlib has the key first); taken here so
                # that the key does not move the focus to the next control (dialog navigation)
            tab.step_cross(n, dj)
            return True
        if ctx.view != "ms" or dj or getattr(tab, "data", None) is None:
            return False
        import numpy as np
        d = tab.data
        if getattr(tab, "_averaging", False):
            return True
        idx = d.event_scans(tab.cols[0]["e"])
        rt = np.sort(np.asarray(d.rt[idx], float))
        if not len(rt):
            return True
        t = tab.pick_t
        if t is None and tab.avg:
            k = int(np.searchsorted(rt, tab.avg[1] if n > 0 else tab.avg[0])) + (0 if n > 0 else -1)
        elif t is None:
            k = 0 if n > 0 else len(rt) - 1
        else:
            k = int(np.argmin(np.abs(rt - t))) + n
        k = min(max(k, 0), len(rt) - 1)
        tab.chrom_pick(float(rt[k]), None, None)
        return True
    return h


def _delete_peak(ctx):
    tab = ctx.tab
    if ctx.view not in ("ms", "pda") or tab.sel_peak is None:
        return False
    p = tab.peaks[tab.sel_peak] if 0 <= tab.sel_peak < len(tab.peaks) else None
    tab.delete_peak(tab.sel_peak)
    if p is not None:
        tab.status("Peak at %.3f min deleted" % p["rt"])
    return True


def _escape(ctx):
    tab = ctx.tab
    if tab is None:
        return False
    if getattr(tab, "tool", "select") != "select":
        tab.tools.set_mode("select")
        tab.set_tool("select")
        return True
    cols = list(getattr(tab, "cols", None) or []) + list((getattr(tab, "mzcols", None) or {}).values())
    for col in cols:  # a measurement waiting for its second peak
        if col.get("m1") is not None:
            col.update(m1=None, measure=None)
            card = col.get("spec_card")
            if ctx.view == "compare":
                for p in ("+", "-"):
                    tab.plot_mz(p, keep_view=True)
            elif card is not None:
                tab.plot_spec(col, keep_view=True)
            tab.status("Measurement cancelled")
            return True
    return False


def _tool(key):
    def h(ctx):
        tab = ctx.tab
        if tab is None or not hasattr(tab, "tools"):
            return False
        b = tab.tools.buttons.get(key)
        if b is None or not b.IsShown():
            return False
        if key == "mzspec":
            tab.tool_action("mzspec")
            return True
        tab.tools.set_mode(key)
        tab.set_tool(key)
        return True
    return h


def _k(c):
    return ord(c.upper())


# ------------------------------------------------------------------ registry
ALL = {"ms", "pda", "compare"}
MSPDA = {"ms", "pda"}
REGISTRY = [
    Shortcut("Ctrl+Shift+S", "Save the analysis", "window", lambda ctx: __import__('ms_project').save(ctx.frame),
             [(CTRL | SHIFT, _k("s"))], accel=True),
    # window: Ctrl combinations (accelerator table), also while typing in a field
    Shortcut("Ctrl+O", "Open a LabSolutions .lcd file", "window", _open, [(CTRL, _k("o"))], accel=True,
             windows=("lcms",), menu=("Open raw data",)),
    Shortcut("Ctrl+O", "Open a Bruker .d folder", "window", _open, [(CTRL, _k("o"))], accel=True,
             windows=("hrms",), menu=("Open a Bruker .d folder",)),
    Shortcut("Ctrl+W", "Close the file shown", "window", _close, [(CTRL, _k("w"))], accel=True),
    Shortcut("Ctrl+R", "Create a report", "window", _report, [(CTRL, _k("r"))], accel=True),
    Shortcut("Ctrl+Tab, Ctrl+Page Down", "Next open file", "window", _next(1),
             [(CTRL, wx.WXK_TAB), (CTRL, wx.WXK_PAGEDOWN), (CTRL, wx.WXK_NUMPAD_PAGEDOWN)], accel=True),
    Shortcut("Ctrl+Shift+Tab, Ctrl+Page Up", "Previous open file", "window", _next(-1),
             [(CTRL | SHIFT, wx.WXK_TAB), (CTRL, wx.WXK_PAGEUP), (CTRL, wx.WXK_NUMPAD_PAGEUP)], accel=True),
    Shortcut("Ctrl+1", "Mass spectrometry view", "window", _show_view(0), [(CTRL, _k("1"))], accel=True),
    Shortcut("Ctrl+2", "PDA view", "window", _show_view(1), [(CTRL, _k("2"))], accel=True, windows=("lcms",)),
    Shortcut("Ctrl+3", "Compare view", "window", _show_view(2), [(CTRL, _k("3"))], accel=True, windows=("lcms",)),
    Shortcut("Ctrl+E", "Export the data of the plot clicked", "window", _export,
             [(CTRL, _k("e"))], accel=True, menu=("Export chromatogram", "Export spectrum", "Export the traces")),
    Shortcut("Ctrl+D", "Deconvolute the spectrum shown", "window", _deconvolute,
             [(CTRL, _k("d"))], accel=True, menu=("Deconvolute this spectrum",)),
    Shortcut("Ctrl+I", "Integrate the peaks in the visible window", "window", _integrate,
             [(CTRL, _k("i"))], views=MSPDA, accel=True),
    # undo and redo (in a text field: the typing in that field)
    Shortcut("Ctrl+Z", "Undo the last step", "undo", _undo, [(CTRL, _k("z"))]),
    Shortcut("Ctrl+Y, Ctrl+Shift+Z", "Redo the step undone", "undo", _redo,
             [(CTRL, _k("y")), (CTRL | SHIFT, _k("z"))]),
    Shortcut("F1", "This list of keyboard shortcuts", "window", _help, [(0, wx.WXK_F1)]),
    Shortcut("F9", "Show or hide the side panel", "window", _side, [(0, wx.WXK_F9)]),
    # view keys (not in a field or a list)
    Shortcut("Esc", "Select tool; cancels an open measurement", "view", _escape,
             [(0, wx.WXK_ESCAPE)]),
    Shortcut("Home", "Full view of every plot", "view", _full, [(0, wx.WXK_HOME), (0, wx.WXK_NUMPAD_HOME)],
             menu=("Full view",)),
    Shortcut("+", "Zoom in on the time or m/z axis of the plot clicked", "view", _zoom(0.8),
             [(0, ord("+")), (SHIFT, ord("=")), (0, ord("=")), (0, wx.WXK_NUMPAD_ADD), (SHIFT, ord("+"))]),
    Shortcut("-", "Zoom out on the time or m/z axis of the plot clicked", "view", _zoom(1.25),
             [(0, ord("-")), (0, wx.WXK_NUMPAD_SUBTRACT)]),
    Shortcut("Left / Right", "Spectra at the previous or next scan (Shift: 10 scans)", "view", None,
             [], views={"ms"}),
    Shortcut("Left / Right, Up / Down", "Move the time or the wavelength slider (Shift: 10 steps)", "view",
             None, [], views={"pda"}),
    Shortcut("Delete", "Delete the integrated peak chosen", "view", _delete_peak,
             [(0, wx.WXK_DELETE), (0, wx.WXK_NUMPAD_DELETE)], views=MSPDA),
    # (handled by the plot itself: matplotlib gets the key first, see unilcms.PlotCard._hotkey)
    Shortcut("Ctrl+C", "Copy the plot clicked as an image", "plot", None,
             [(CTRL, _k("c"))], where="every plot (a plot clicked)"),
    Shortcut("Ctrl+S", "Save the plot clicked as an image", "plot", None,
             [(CTRL, _k("s"))], where="every plot (a plot clicked)"),
    # tools: single keys while a plot (or the tool bar) has the focus
    Shortcut("S", "Select tool", "plot", _tool("select"), [(0, _k("s"))], tool="select"),
    Shortcut("Z", "Zoom tool", "plot", _tool("zoom"), [(0, _k("z"))], tool="zoom"),
    Shortcut("P", "Pan tool", "plot", _tool("pan"), [(0, _k("p"))], tool="pan"),
    Shortcut("B", "Background tool", "plot", _tool("background"), [(0, _k("b"))],
             views=MSPDA, tool="background"),
    Shortcut("I", "Integration tool Click", "plot", _tool("click_peak"), [(0, _k("i"))],
             views=MSPDA, tool="click_peak"),
    Shortcut("D", "Integration tool Drag", "plot", _tool("drag_peak"), [(0, _k("d"))],
             views=MSPDA, tool="drag_peak"),
    Shortcut("K", "Integration tool Split", "plot", _tool("split"), [(0, _k("k"))],
             views=MSPDA, tool="split"),
    Shortcut("R", "Integration tool Delete", "plot", _tool("delete"), [(0, _k("r"))],
             views=MSPDA, tool="delete"),
    Shortcut("X", "Mass chromatogram tool", "plot", _tool("xic"), [(0, _k("x"))],
             views={"ms", "compare"}, tool="xic"),
    Shortcut("U", "Measure tool", "plot", _tool("measure"), [(0, _k("u"))],
             views={"ms", "compare"}, tool="measure"),
    Shortcut("L", "Label tool", "plot", _tool("label"), [(0, _k("l"))],
             tool="label"),
    Shortcut("F", "Formula tool", "plot", _tool("formula"), [(0, _k("f"))], views={"ms"},
             windows=("hrms",), tool="formula"),
    Shortcut("A", "Align tool", "plot", _tool("align"), [(0, _k("a"))],
             views={"compare"}, tool="align"),
    Shortcut("G", "Guide line tool", "plot", _tool("guide"), [(0, _k("g"))], views={"compare"}, tool="guide"),
    Shortcut("M", "m/z tool on or off", "plot", _tool("mzspec"), [(0, _k("m"))],
             views={"compare"}, tool="mzspec"),
]
# arrows (their combinations are added here: one entry each for F1 above)
_ARROWS = [
    Shortcut("Left", "", "view", _step(-1), [(0, wx.WXK_LEFT), (0, wx.WXK_NUMPAD_LEFT)], views=MSPDA),
    Shortcut("Right", "", "view", _step(1), [(0, wx.WXK_RIGHT), (0, wx.WXK_NUMPAD_RIGHT)], views=MSPDA),
    Shortcut("Shift+Left", "", "view", _step(-10), [(SHIFT, wx.WXK_LEFT), (SHIFT, wx.WXK_NUMPAD_LEFT)], views=MSPDA),
    Shortcut("Shift+Right", "", "view", _step(10), [(SHIFT, wx.WXK_RIGHT), (SHIFT, wx.WXK_NUMPAD_RIGHT)],
             views=MSPDA),
    Shortcut("Up", "", "view", _step(0, 1), [(0, wx.WXK_UP), (0, wx.WXK_NUMPAD_UP)], views={"pda"}),
    Shortcut("Down", "", "view", _step(0, -1), [(0, wx.WXK_DOWN), (0, wx.WXK_NUMPAD_DOWN)], views={"pda"}),
    Shortcut("Shift+Up", "", "view", _step(0, 10), [(SHIFT, wx.WXK_UP), (SHIFT, wx.WXK_NUMPAD_UP)],
             views={"pda"}),
    Shortcut("Shift+Down", "", "view", _step(0, -10), [(SHIFT, wx.WXK_DOWN), (SHIFT, wx.WXK_NUMPAD_DOWN)],
             views={"pda"}),
]
# keys handled by their own controls, and the mouse: listed in F1
NATIVE = [
    ("Delete, Backspace", "Delete the peak chosen in the peak table", "peak table (a row chosen)"),
    ("Space, Delete", "Untick the file chosen (Space ticks and unticks it)", "Compare view, list of files"),
    ("Up / Down", "Previous or next deconvolution result", "list of open files"),
    ("Enter", "Show the result at the cursor in the table", "list of open files"),
    ("Space", "Hide the result or show it again", "list of open files"),
    ("Delete", "Close the result", "list of open files"),
    ("Left / Right", "Fold or unfold the results of a file", "list of open files"),
    ("Ctrl+A", "Choose every row", "tables of the deconvolution"),
    ("Delete", "Leave out the calibrant point chosen", "calibration window"),
    ("Esc", "Close the window (Cancel)", "dialogs (settings, graph properties, calibration)"),
    ("Enter", "OK or Apply", "dialogs"),
    ("Ctrl+Z", "Undo the typing", "a text field"),
]
MOUSE = [
    ("drag", "Average a time range (chromatograms); zoom the m/z axis (spectra)", "Mass spectrometry and PDA views"),
    ("Shift + drag", "Background range", "chromatograms; Compare view with m/z on"),
    ("Ctrl + drag", "Zoom a box (any tool)", "every plot"),
    ("double click", "Full view of the plot", "every plot"),
    ("wheel", "Scroll the column of tiles", "every view"),
    ("Ctrl + wheel", "Zoom the time or m/z axis", "every plot"),
    ("Ctrl + Shift + wheel", "Zoom the wavelength axis", "PDA map"),
    ("Alt + click", "One scan at that time (m/z tool)", "Compare view"),
    ("right click", "Menu of the plot", "every plot"),
    ("drag a gap", "Resize a tile", "every view"),
    ("double click a gap", "Fit every tile in the window", "every view"),
]


def _match(sc, mods, code):
    for m, c in sc.combos:
        if m == mods and c == code:
            return True
    return False


def _applies(sc, ctx, wk):
    if wk not in sc.windows:
        return False
    if sc.views is not None and ctx.view not in sc.views:
        return False
    return True


_VKSCAN = {}


def _layout_char(e, mods):
    """'+' or '-' when this key with these modifiers types it on the keyboard layout in use (Windows): the key
    code of EVT_CHAR_HOOK is the key, not the character (a Swiss or German keyboard types + with Shift+1 or
    its own key). None otherwise."""
    if not sys.platform.startswith("win"):
        return None
    try:
        import ctypes
        user32 = ctypes.windll.user32
        if not _VKSCAN:
            user32.GetKeyboardLayout.restype = ctypes.c_void_p
            user32.GetKeyboardLayout.argtypes = [ctypes.c_uint32]
            user32.VkKeyScanExW.restype = ctypes.c_short
            user32.VkKeyScanExW.argtypes = [ctypes.c_wchar, ctypes.c_void_p]
            _VKSCAN["init"] = True
        hkl = user32.GetKeyboardLayout(0)
        raw = e.GetRawKeyCode()
        for ch in "+-":
            key = (hkl, ch)
            if key not in _VKSCAN:
                _VKSCAN[key] = user32.VkKeyScanExW(ch, hkl) & 0xFFFF  # low byte: the key, high byte: Shift 1,
            v = _VKSCAN[key]                                          # Ctrl 2, Alt 4
            if v == 0xFFFF:
                continue
            st = (v >> 8) & 0xFF
            need = (SHIFT if st & 1 else 0) | (CTRL if st & 2 else 0) | (ALT if st & 4 else 0)
            if raw == (v & 0xFF) and mods == need:
                return ch
    except Exception:
        return None
    return None


def on_char_hook(frame, e):
    """EVT_CHAR_HOOK of a Postrun window: the key goes to its shortcut (if any) before the focused control."""
    try:
        # a key of another top level window of it (the Polymer and Kinetics windows are wx.Frame children, whose
        # char hook comes up here too) is theirs: its grids, fields and plots get it, never the Postrun views
        w = e.GetEventObject() if isinstance(e.GetEventObject(), wx.Window) else wx.Window.FindFocus()
        if w is not None and wx.GetTopLevelParent(w) is not frame:
            e.Skip()
            return
    except Exception:
        pass
    try:
        code, mods = e.GetKeyCode(), e.GetModifiers()
        if code in (wx.WXK_CONTROL, wx.WXK_SHIFT, wx.WXK_ALT, wx.WXK_WINDOWS_LEFT, wx.WXK_WINDOWS_RIGHT):
            e.Skip()
            return
        mods &= (CTRL | SHIFT | ALT)
        ch = _layout_char(e, mods)
        if ch is not None:  # + and - as the keyboard types them (the numeric keypad has its own codes)
            code, mods = ord(ch), 0
        ctx = Ctx(frame)
        wk = window_kind(frame)
        for sc in REGISTRY + _ARROWS:
            if sc.accel or sc.handler is None or not _match(sc, mods, code) or not _applies(sc, ctx, wk):
                continue
            if sc.scope == "undo" and ctx.kind == "text":
                continue  # the field undoes its own typing
            if sc.scope == "view" and ctx.kind not in ("plot", "view", "none"):
                continue
            if sc.scope == "plot" and (ctx.kind not in ("plot", "view", "none") or
                                       (sc.tool is None and ctx.card is None)):
                continue
            if sc.handler(ctx) is False:
                continue
            if sc.scope != "undo":
                import undo
                undo.schedule(frame)
            return
        if ctx.kind == "plot" and code in _NAV_KEYS and not mods & (CTRL | ALT):
            return  # (an arrow key on a plot: never the dialog navigation to another control)
    except Exception as ex:
        print("[shortcuts] %r" % ex)
    e.Skip()


_NAV_KEYS = (wx.WXK_LEFT, wx.WXK_RIGHT, wx.WXK_UP, wx.WXK_DOWN, wx.WXK_NUMPAD_LEFT, wx.WXK_NUMPAD_RIGHT,
             wx.WXK_NUMPAD_UP, wx.WXK_NUMPAD_DOWN)


def install(frame):
    """The accelerator table (Ctrl combinations) and the key dispatch of a Postrun window."""
    wk = window_kind(frame)
    acc = []
    for sc in REGISTRY:
        if not sc.accel or wk not in sc.windows:
            continue
        for mods, code in sc.combos:
            i = wx.NewIdRef()

            def run(e, sc=sc):
                ctx = Ctx(frame)
                if sc.views is not None and ctx.view not in sc.views:
                    return
                try:
                    sc.handler(ctx)
                finally:
                    import undo
                    undo.schedule(frame)
            frame.Bind(wx.EVT_MENU, run, id=i)
            flags = (wx.ACCEL_CTRL if mods & CTRL else 0) | (wx.ACCEL_SHIFT if mods & SHIFT else 0) | \
                (wx.ACCEL_ALT if mods & ALT else 0)
            acc.append((flags or wx.ACCEL_NORMAL, code, i))
    frame.SetAcceleratorTable(wx.AcceleratorTable(acc))
    frame.Bind(wx.EVT_CHAR_HOOK, lambda e: on_char_hook(frame, e))


def plot_clicked(card):
    """A click on a plot: it takes the keyboard focus (tool letters, Ctrl+C) and is the plot of Ctrl+E."""
    try:
        if wx.Window.FindFocus() is not card.canvas:
            card.canvas.SetFocus()
        top = wx.GetTopLevelParent(card)
        top.__dict__["_last_card"] = card
    except Exception:
        pass


# ------------------------------------------------------------------ tooltips and menus
def keys_for_tool(tab, key):
    wk = window_kind(wx.GetTopLevelParent(tab))
    vk = view_kind(tab)
    for sc in REGISTRY:
        if sc.tool == key and wk in sc.windows and (sc.views is None or vk in sc.views):
            return sc.keys
    return {"full": "Home", "auto": "Ctrl+I", "deconvolute": "Ctrl+D"}.get(key)


def decorate_tools(tab):
    """The keys of the tools in the tooltips of the tool bar (e.g. "Zoom: drag a box (Z)")."""
    tools = getattr(tab, "tools", None)
    if tools is None:
        return
    for key, b in tools.buttons.items():
        k = keys_for_tool(tab, key)
        if not k:
            continue
        try:
            tip = b.GetToolTipText() or ""
        except Exception:
            continue
        if "(%s" % k in tip or (" %s)" % k) in tip:
            continue
        if tip.endswith("(Esc)"):
            tip = tip[:-5] + "(%s or Esc)" % k
        else:
            tip = "%s (%s)" % (tip, k) if tip else k
        b.SetToolTip(tip)


def menu_label(label):
    """A menu entry with the key of its shortcut after a tab (shown at the right of the menu)."""
    if not label or "\t" in label:
        return label
    for sc in REGISTRY:
        for start in sc.menu or ():
            if label.startswith(start):
                return "%s\t%s" % (label, sc.keys.split(",")[0])
    return label


# ------------------------------------------------------------------ F1 window
def rows(wk=None):
    """(window, keys, action, where) of every shortcut (wk: one window)."""
    out = []
    for sc in REGISTRY:
        if not sc.action:
            continue
        for w in sc.windows:
            if wk is None or w == wk:
                out.append((WINDOWS[w], sc.keys, sc.action, sc.where_text()))
    for w in (WINDOWS if wk is None else [wk]):
        for keys, act, where in NATIVE:
            if "Compare" in where and w == "hrms":
                continue
            if "calibration" in where and w == "lcms":
                continue
            out.append((WINDOWS[w], keys, act, where))
        for keys, act, where in MOUSE:
            if ("PDA" in where or "Compare" in where) and w == "hrms":
                continue
            out.append((WINDOWS[w], keys, act, where + " (mouse)"))
    return out


class ShortcutsWindow(wx.Frame):
    """F1: every shortcut, by window, with a search field."""

    def __init__(self, parent):
        wx.Frame.__init__(self, parent, title="Keyboard shortcuts",
                          style=wx.DEFAULT_FRAME_STYLE | wx.FRAME_FLOAT_ON_PARENT)
        self.wk = window_kind(parent)
        p = wx.Panel(self)
        self.search = wx.SearchCtrl(p, style=wx.TE_PROCESS_ENTER)
        self.search.ShowCancelButton(True)
        self.search.SetDescriptiveText("Search")
        self.which = wx.Choice(p, choices=["This window (%s)" % WINDOWS[self.wk], "Every window"])
        self.which.SetSelection(0)
        self.list = wx.ListCtrl(p, style=wx.LC_REPORT | wx.LC_SINGLE_SEL | wx.BORDER_THEME)
        for i, (name, w) in enumerate((("Keys", 170), ("Action", 420), ("Where", 300), ("Window", 110))):
            self.list.InsertColumn(i, name, width=self.FromDIP(w))
        top = wx.BoxSizer(wx.HORIZONTAL)
        top.Add(self.search, 1, wx.RIGHT, self.FromDIP(8))
        top.Add(self.which, 0)
        s = wx.BoxSizer(wx.VERTICAL)
        s.Add(top, 0, wx.EXPAND | wx.ALL, self.FromDIP(10))
        s.Add(self.list, 1, wx.EXPAND | wx.LEFT | wx.RIGHT, self.FromDIP(10))
        s.AddSpacer(self.FromDIP(10))
        p.SetSizer(s)
        self.search.Bind(wx.EVT_TEXT, lambda e: self.fill())
        self.search.Bind(wx.EVT_SEARCHCTRL_CANCEL_BTN, lambda e: (self.search.SetValue(""), self.fill()))
        self.which.Bind(wx.EVT_CHOICE, lambda e: self.fill())
        self.Bind(wx.EVT_CHAR_HOOK, self._key)
        self.SetSize(self.FromDIP(wx.Size(1040, 640)))
        self.CentreOnParent()
        self.fill()
        self.search.SetFocus()

    def _key(self, e):
        if e.GetKeyCode() == wx.WXK_ESCAPE:
            self.Close()
        else:
            e.Skip()

    def fill(self):
        q = self.search.GetValue().strip().lower()
        rs = rows(self.wk if self.which.GetSelection() == 0 else None)
        if q:
            words = q.split()
            rs = [r for r in rs if all(w in " ".join(r).lower() for w in words)]
        self.list.Freeze()
        try:
            self.list.DeleteAllItems()
            for w, keys, act, where in rs:
                i = self.list.InsertItem(self.list.GetItemCount(), keys)
                self.list.SetItem(i, 1, act)
                self.list.SetItem(i, 2, where)
                self.list.SetItem(i, 3, w)
        finally:
            self.list.Thaw()
        self.shown = rs


def show_window(frame):
    w = frame.__dict__.get("_shortcuts_window")
    try:
        if w:
            w.Raise()
            w.search.SetFocus()
            return w
    except RuntimeError:
        pass
    w = ShortcutsWindow(frame)
    frame.__dict__["_shortcuts_window"] = w
    w.Show()
    return w

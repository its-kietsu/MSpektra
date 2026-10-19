"""The menu bar and the toolbar of LCMS and HRMS Analysis (ribbon.py): which commands they show and call.

The top bar (Open, Save, Undo, Redo, the views, the chips, Deconvolute, Report, the side panel, help) and the
tool strip of every view are still made, hidden; the ribbon calls them and shows their state, so that the
shortcuts, menus and undo keep working as before. Settings: "ribbon": false gives the old top bar back."""
import wx

import unidec_theme as T
import ribbon as R

SHORT = {"Mass spectrometry": "MS"}  # labels of the view buttons
VIEW_ICONS = {"Mass spectrometry": "ms", "MS": "ms", "PDA (UV/Vis)": "pda", "Compare": "compare"}
MS, PDA, CMP, HR = ("ms", "hrms"), ("pda",), ("compare",), ("hrms",)
CHROM_VIEWS = ("ms", "hrms", "pda")


def enabled(frame):
    return T._load().get("ribbon", True) is not False


# ------------------------------------------------------------------ the view shown and its tools
def current_view(fr):
    import shortcuts
    return shortcuts.view_of(fr)


def kind(fr):
    """"ms", "hrms", "pda" or "compare": the view shown."""
    import shortcuts
    v = current_view(fr)
    k = shortcuts.view_kind(v)
    if k == "ms" and hasattr(v, "cal"):
        return "hrms"
    return k


def strip(fr):
    v = current_view(fr)
    return getattr(v, "tools", None) if v is not None else None


def _btn(fr, key):
    s = strip(fr)
    try:
        return s.buttons.get(key) if s is not None else None
    except RuntimeError:
        return None


def _safe(f, *a):
    try:
        return f(*a)
    except Exception:
        return False


def tool_item(fr, key, label, icon=None, size="large"):
    """A tool of the view shown (its hidden tool strip); shown only when that view has it."""
    def run(b):
        s = strip(fr)
        if s is not None and key in s.buttons:
            s._clicked(key)

    def tip():
        b = _btn(fr, key)
        return b.GetToolTipText() if b is not None else label

    def shown():
        b = _btn(fr, key)
        return b is not None and b.IsShown()

    def active():
        b = _btn(fr, key)
        return bool(b is not None and b.active)

    return R.Item(key, label, icon or key, run, size=size, tip=tip, shown=shown, active=active)


def menu_at(button, entries):
    """A menu of (label, function) under a ribbon button."""
    m = wx.Menu()
    acts = {}
    for label, fn in entries:
        it = m.Append(wx.ID_ANY, label)
        acts[it.GetId()] = fn
    try:
        sel = button.GetPopupMenuSelectionFromUser(m, wx.Point(0, button.GetSize().height))
    finally:
        m.Destroy()
    if sel in acts:
        wx.CallAfter(acts[sel])


def _method_item(fr):
    return R.Item("method", "Method", "method", lambda b: __import__("method_presets").show_menu(fr, b), menu=True,
                  tip="Method presets", shown=lambda: "method" in getattr(fr, "chips", {}))


def _deconv_item(fr):
    return R.Item("deconvolute", "Deconvolute", "deconvolve", lambda b: fr.deconvolute(),
                  tip="Deconvolute the spectrum shown (Ctrl+D)", shown=lambda: kind(fr) in MS)


def _report_item(fr):
    return R.Item("report", "Report", "report", lambda b: fr.make_report(), tip="PDF or Word report (Ctrl+R)")


# ------------------------------------------------------------------ tabs
def home_tab(fr):
    import undo
    ub = fr.__dict__.get("_undo_buttons") or (None, None)

    def undo_item(i, key, label):
        b = ub[i]
        return R.Item(key, label, key, lambda _b: (undo.do_undo if i == 0 else undo.do_redo)(fr), size="small",
                      tip=lambda: b.GetToolTipText() if b else label, enabled=lambda: bool(b and b.IsEnabled()))

    return R.Tab("Home", [  # (Open, Save and Save as: in the File menu)
        R.Group("Edit", [undo_item(0, "undo", "Undo"), undo_item(1, "redo", "Redo")]),
        R.Group("Navigate", [tool_item(fr, "select", "Select"), tool_item(fr, "zoom", "Zoom"),
                             tool_item(fr, "pan", "Pan"), tool_item(fr, "full", "Fit")]),
        R.Group("Peaks", [tool_item(fr, "auto", "Auto"), tool_item(fr, "drag_peak", "Drag"),
                          tool_item(fr, "click_peak", "Click"), tool_item(fr, "split", "Split"),
                          tool_item(fr, "delete", "Delete"), tool_item(fr, "clear", "Clear")]),
        R.Group("Chromatogram", [tool_item(fr, "background", "Background"), tool_item(fr, "xic", "Mass chrom.")]),
        R.Group("Compare", [tool_item(fr, "align", "Align"), tool_item(fr, "guide", "Guide line"),
                            tool_item(fr, "area", "Area"), tool_item(fr, "mzspec", "m/z")]),
        R.Group("Annotate", [tool_item(fr, "measure", "Measure"), tool_item(fr, "label", "Label"),
                             tool_item(fr, "formula", "Formula")]),
        R.Group("Process", [
            _method_item(fr), tool_item(fr, "calibrate", "Calibrate"), _deconv_item(fr),
            R.Item("polymer", "Polymer", "polymer",
                   lambda b: menu_at(b, __import__("ms_polymer_ui").menu_items(current_view(fr))), menu=True,
                   tip="Mass distribution, repeat units, chain lengths", shown=lambda: _btn(fr, "polymer") is not None),
            _report_item(fr)])])


def view_group(fr):
    """MS, PDA and Compare: large buttons pinned at the start of the ribbon (the view switch of the top bar)."""
    seg = fr.seg

    def item(i):
        lab = seg.labels[i] if i < len(seg.labels) else ""

        def run(b):
            if i < len(seg.labels) and seg.enabled[i] and i != seg.active:
                seg.active = i
                seg.Refresh()
                seg.on_change(i)

        return R.Item("view%d" % i, SHORT.get(lab, lab), VIEW_ICONS.get(lab, "ms"), run,
                      tip=lab, shown=lambda: i < len(seg.labels) and len(seg.labels) > 1,
                      enabled=lambda: i < len(seg.enabled) and seg.enabled[i], active=lambda: seg.active == i)

    return R.Group("Data", [item(i) for i in range(len(seg.labels))])


def tabs(fr):
    return [home_tab(fr)]


# ------------------------------------------------------------------ the menu bar (as in Bruker DataAnalysis)
def _cmd(label, fn, ok=True, keys=None):
    """A menu entry: greyed when not ok (the view shown has no such command)."""
    return (label + ("\t" + keys if keys else ""), fn if ok else None)


def _view_cmd(fr, label, method, args=(), keys=None, kinds=None):
    """A command of the view shown; left out of the menu in other kinds of views."""
    if kinds is not None and kind(fr) not in kinds:
        return None
    v = current_view(fr)
    f = getattr(v, method, None) if v is not None else None
    ok = f is not None and (kinds is None or kind(fr) in kinds)
    return _cmd(label, (lambda: f(*args)) if ok else None, ok, keys)


def _tool(fr, key, label, keys=None):
    """A tool of the view shown, ticked while it is chosen."""
    b = _btn(fr, key)
    if b is None or not b.IsShown():
        return None  # (not a tool of this view)
    lab = label + ("\t" + keys if keys else "")
    if b.mode:
        return (lab, lambda: strip(fr)._clicked(key), bool(b.active))
    return (lab, lambda: strip(fr)._clicked(key))


def menu_file(fr):
    import ms_project
    items = [(lab, (lambda f=f: wx.CallAfter(f))) for lab, f in fr.raw_open_items()]
    if items:
        items[0] = (items[0][0] + "\tCtrl+O", items[0][1])
    has = bool(getattr(fr.active, "path", None))
    return items + [
        _cmd("Open a saved analysis…", lambda: ms_project.open_file(fr)),
        (None, None),
        _cmd("Save", lambda: ms_project.save(fr), has, "Ctrl+Shift+S"),
        _cmd("Save as…", lambda: ms_project.save_as(fr), has),
        (None, None),
        _cmd("Report…", lambda: fr.make_report(), has, "Ctrl+R"),
        (None, None),
        _cmd("Close the file", lambda: fr.close_doc(fr.active), has, "Ctrl+W"),
        _cmd("Close every file", lambda: fr.close_all(), has),
        (None, None),
        _cmd("Exit", lambda: fr.Close())]


def menu_edit(fr):
    import undo
    m = undo.manager(fr)
    u = m.undo_label() if m is not None else None
    r = m.redo_label() if m is not None else None
    import unilcms
    return [
        _cmd("Undo" + (" " + u if u else ""), lambda: undo.do_undo(fr), u is not None, "Ctrl+Z"),
        _cmd("Redo" + (" " + r if r else ""), lambda: undo.do_redo(fr), r is not None, "Ctrl+Y"),
        (None, None),
        _cmd("Image size for copy and save…", lambda: unilcms.ask_export_settings(fr))]


def menu_view(fr):
    seg = fr.seg
    out = []
    if len(seg.labels) > 1:
        for i, lab in enumerate(seg.labels):
            def go(i=i):
                if seg.enabled[i] and i != seg.active:
                    seg.active = i
                    seg.Refresh()
                    seg.on_change(i)
            out.append((lab + "\tCtrl+%d" % (i + 1), go if seg.enabled[i] else None, seg.active == i))
        out.append((None, None))
    side = fr.__dict__.get("_side_shown")
    if kind(fr) == "compare":
        side = bool(_safe(lambda: current_view(fr).side.IsShown()))
    fp = fr.files

    def files_toggle():
        fp.collapsed = not fp.collapsed
        T._save({"files_collapsed": fp.collapsed})
        fp._apply_width()
        fp.GetParent().Layout()
        fp.Refresh()

    return out + [
        ("Side panel\tF9", lambda: fr.set_side(None), bool(side)),
        ("File list", files_toggle, not fp.collapsed),
        (None, None),
        _tool(fr, "full", "Full view", "Home"),
        _view_cmd(fr, "Graph properties…", "on_graph_props", kinds=CMP),
        (None, None),
        _cmd("Next file", lambda: fr.next_doc(1), len(fr.docs) > 1, "Ctrl+Tab"),
        _cmd("Previous file", lambda: fr.next_doc(-1), len(fr.docs) > 1, "Ctrl+Shift+Tab")]


def menu_chromatogram(fr):
    k = kind(fr)
    if k == "compare":
        return [
            _tool(fr, "align", "Align", "A"), _tool(fr, "guide", "Guide line", "G"),
            (None, None),
            _view_cmd(fr, "Export the traces (CSV)…", "on_export_traces"),
            _view_cmd(fr, "Export the table (CSV)…", "on_export_table")]
    return [
        _view_cmd(fr, "Average a time range…", "ask_average"),
        (None, None),
        _tool(fr, "background", "Background range", "B"),
        _view_cmd(fr, "Clear the background range", "clear_bg", kinds=MS),
        (None, None),
        _tool(fr, "xic", "Mass chromatogram tool", "X"),
        _view_cmd(fr, "Add a mass chromatogram…", "ask_xic", (None,), kinds=MS),
        _view_cmd(fr, "Remove every mass chromatogram", "clear_xic", (None,), kinds=MS),
        _view_cmd(fr, "Find mass…", "ask_find_mass", kinds=("ms",)),
        (None, None),
        _view_cmd(fr, "Go to a time and wavelength…", "ask_goto", kinds=PDA),
        _view_cmd(fr, "Map and chromatogram settings…", "ask_settings", kinds=PDA)]


def menu_peaks(fr):
    if kind(fr) == "compare":  # the peak area of a region in every run
        return [
            _tool(fr, "area", "Area of a region"),
            _view_cmd(fr, "Use the region", "on_area_fields"),
            _view_cmd(fr, "Clear the region", "clear_area"),
            (None, None),
            _view_cmd(fr, "Area vs X…", "plot_area_graph"),
            _view_cmd(fr, "Kinetic fitting…", "open_kinetics"),
            _view_cmd(fr, "Peak areas to Excel…", "on_export_areas_xlsx")]
    return [
        _tool(fr, "auto", "Integrate", "Ctrl+I"),
        _view_cmd(fr, "Integrate with options…", "ask_integrate", kinds=CHROM_VIEWS),
        _view_cmd(fr, "Add the selected range as a peak", "on_add_peak", kinds=CHROM_VIEWS),
        (None, None),
        _tool(fr, "drag_peak", "Drag a peak", "D"), _tool(fr, "click_peak", "Click a peak", "I"),
        _tool(fr, "split", "Split a peak", "K"), _tool(fr, "delete", "Delete a peak", "R"),
        _tool(fr, "clear", "Remove every peak"),
        (None, None),
        _view_cmd(fr, "Export the peak table…", "on_export_table", kinds=CHROM_VIEWS)]


def menu_spectrum(fr):
    v = current_view(fr)
    hr = kind(fr) == "hrms"
    cal = getattr(v, "cal", None)
    out = []
    if kind(fr) == "compare":  # (the spectra of the Compare view: with m/z on)
        out = [_tool(fr, "mzspec", "m/z spectra", "M"), _tool(fr, "xic", "Mass chromatogram tool", "X")]
    out += [_tool(fr, "measure", "Measure", "U"), _tool(fr, "label", "Label", "L")]
    if hr:
        out += [_tool(fr, "formula", "Formula", "F"), (None, None),
                _tool(fr, "calibrate", "Calibrate…"),
                _cmd("Calibration details…", lambda: cal.open_details(), cal is not None),
                _cmd("Load a calibration…", lambda: cal.load_cal(), cal is not None),
                _cmd("Remove the calibration", lambda: cal.remove(), cal is not None)]
    return out + [
        (None, None),
        _view_cmd(fr, "Export the spectrum…", "on_export_spectrum", kinds=MS) if kind(fr) != "pda" else
        _view_cmd(fr, "Export the UV spectrum…", "on_export_uv", kinds=PDA),
        _view_cmd(fr, "Copy the spectrum data", "copy_spec", (None,), kinds=MS)]


def menu_deconvolution(fr):
    v = current_view(fr)
    dec = getattr(v, "dec", None)
    ms = kind(fr) in MS
    out = [_cmd("Deconvolute…", lambda: fr.deconvolute(), ms, "Ctrl+D"),
           _cmd("Mass shifts…", lambda: dec.shifts.open_dialog(), ms and dec is not None),
           _cmd("Close the results", lambda: dec.close(), ms and dec is not None and bool(dec.results))]
    if getattr(v, "FULL_WINDOW", False):
        out.append(_view_cmd(fr, "Open in the Deconvolute window", "open_full", kinds=("ms",)))
    if kind(fr) == "hrms":
        import ms_polymer_ui
        out += [(None, None), ("Polymer", [(lab, fn) for lab, fn in ms_polymer_ui.menu_items(v)])]
    return out


def menu_method(fr):
    import method_presets
    return method_presets.menu_items(fr)


def menu_help(fr):
    import shortcuts
    out = [_cmd("Keyboard shortcuts", lambda: shortcuts.show_window(fr), keys="F1"),
           _cmd("About " + fr.TITLE, lambda: fr.on_help())]
    try:
        import updater
        prev = updater.previous()
        if prev:
            out += [(None, None), _cmd("Back to MSpektra %s" % prev, updater.back_to_previous)]
    except Exception:
        pass
    return out


def menus(fr):
    m = [("File", menu_file), ("Edit", menu_edit), ("View", menu_view), ("Chromatogram", menu_chromatogram),
         ("Spectrum", menu_spectrum), ("Peaks", menu_peaks), ("Deconvolution", menu_deconvolution)]
    if "method" in getattr(fr, "chips", {}):
        m.append(("Method", menu_method))
    m.append(("Help", menu_help))
    return [(t, (lambda f=f: [e for e in f(fr) if e is not None])) for t, f in m]


# ------------------------------------------------------------------ in the window
def install(fr, root, tb, vs):
    """The ribbon in place of the top bar tb (which stays, hidden) in the sizer vs of the window."""
    rb = R.Ribbon(root, tabs(fr), view_group(fr), menus=menus(fr))
    fr.__dict__["ribbon"] = rb
    row = rb.tabrow
    for key, chip in list(getattr(fr, "chips", {}).items()):  # the chips: at the right of the tab row
        try:
            tb.sizer.Detach(chip)
        except Exception:
            pass
        chip.Reparent(row)
        row.add(chip, 3)
    row.sizer.AddSpacer(fr.FromDIP(8))
    fr._top_labels = {}  # (no top bar to fit)
    tb.Hide()
    vs.Insert(0, rb, 0, wx.EXPAND)
    return rb

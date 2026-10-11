"""The ribbon of LCMS and HRMS Analysis (ribbon.py): which commands it shows and what they call.

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


def view_item(fr, key, label, icon, method, kinds, args=(), tip=None, size="large", cond=None):
    """A command of the view shown: view.method(*args) (or method(view) for a function), shown for the kinds of
    views given."""
    def run(b):
        v = current_view(fr)
        if isinstance(method, str):
            f = getattr(v, method, None)
            if f is not None:
                f(*args)
        else:
            method(v)

    def shown():
        if kind(fr) not in kinds:
            return False
        v = current_view(fr)
        if isinstance(method, str) and not hasattr(v, method):
            return False
        return cond is None or bool(_safe(cond, v))

    return R.Item(key, label, icon, run, size=size, tip=tip or label, shown=shown)


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


class _Ev:
    def __init__(self, obj):
        self.obj = obj

    def GetEventObject(self):
        return self.obj


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

    return R.Tab("Home", [
        R.Group("File", [
            R.Item("open", "Open", "folder", lambda b: fr.on_open_menu(_Ev(b)), menu=True,
                   tip="Open raw data or a saved analysis (Ctrl+O)"),
            R.Item("save", "Save", "save", lambda b: __import__("ms_project").save(fr, None),
                   tip="Save the analysis (Ctrl+Shift+S)"),
            R.Item("save_as", "Save as", "save_as", lambda b: __import__("ms_project").save_as(fr, None),
                   tip="Save the analysis in another file")]),
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


def peaks_tab(fr):
    return R.Tab("Peaks", [
        R.Group("Integrate", [
            tool_item(fr, "auto", "Auto"),
            view_item(fr, "ask_integrate", "Integrate…", "auto", "ask_integrate", CHROM_VIEWS,
                      tip="Integrate the peaks in the visible window, with options"),
            view_item(fr, "add_peak", "Add range", "drag_peak", "on_add_peak", CHROM_VIEWS,
                      tip="Add the selected range as a peak")]),
        R.Group("Edit peaks", [tool_item(fr, "drag_peak", "Drag"), tool_item(fr, "click_peak", "Click"),
                               tool_item(fr, "split", "Split"), tool_item(fr, "delete", "Delete"),
                               tool_item(fr, "clear", "Clear")]),
        R.Group("Background", [tool_item(fr, "background", "Background"),
                               view_item(fr, "clear_bg", "Clear", "clear", "clear_bg", MS,
                                         tip="Clear the background range")]),
        R.Group("Peak table", [view_item(fr, "export_peaks", "Export…", "export", "on_export_table", CHROM_VIEWS,
                                         tip="Export the peak table")]),
        R.Group("Peak area", [
            tool_item(fr, "area", "Area"),
            view_item(fr, "use_region", "Use region", "table", "on_area_fields", CMP),
            view_item(fr, "clear_area", "Clear region", "clear", "clear_area", CMP),
            view_item(fr, "area_graph", "Area vs X…", "xic", "plot_area_graph", CMP),
            view_item(fr, "kinetics", "Kinetics…", "pda", "open_kinetics", CMP, tip="Kinetic fitting"),
            view_item(fr, "areas_xlsx", "Excel…", "export", "on_export_areas_xlsx", CMP,
                      tip="Peak areas to Excel")])])


def spectrum_tab(fr):
    return R.Tab("Spectrum", [
        R.Group("Average", [view_item(fr, "ask_average", "Average…", "ms", "ask_average", CHROM_VIEWS,
                                      tip="Average a time range")]),
        R.Group("Mass chromatograms", [
            tool_item(fr, "xic", "Mass chrom."),
            view_item(fr, "ask_xic", "Add…", "xic", "ask_xic", MS, args=(None,), tip="Add a mass chromatogram"),
            view_item(fr, "clear_xic", "Remove all", "clear", "clear_xic", MS, args=(None,),
                      tip="Remove every mass chromatogram")]),
        R.Group("m/z", [tool_item(fr, "mzspec", "m/z")]),
        R.Group("Annotate", [tool_item(fr, "measure", "Measure"), tool_item(fr, "label", "Label"),
                             tool_item(fr, "formula", "Formula")]),
        R.Group("Calibration", [
            tool_item(fr, "calibrate", "Calibrate"),
            view_item(fr, "cal_details", "Details…", "table", lambda v: v.cal.open_details(), HR,
                      tip="Calibration details"),
            view_item(fr, "cal_load", "Load…", "folder", lambda v: v.cal.load_cal(), HR, tip="Load a calibration"),
            view_item(fr, "cal_remove", "Remove", "clear", lambda v: v.cal.remove(), HR,
                      tip="Remove the calibration")]),
        R.Group("PDA", [view_item(fr, "goto", "Go to…", "select", "ask_goto", PDA,
                                  tip="Go to a time and wavelength"),
                        view_item(fr, "map_settings", "Settings…", "settings", "ask_settings", PDA,
                                  tip="Map and chromatogram settings")]),
        R.Group("Export", [
            view_item(fr, "export_spec", "Export…", "export", "on_export_spectrum", MS, tip="Export the spectrum"),
            view_item(fr, "export_uv", "Export…", "export", "on_export_uv", PDA, tip="Export the UV spectrum"),
            view_item(fr, "copy_spec", "Copy data", "copy", "copy_spec", MS, args=(None,),
                      tip="Copy the spectrum data")])])


def deconv_tab(fr):
    def has_dec(v):
        return getattr(v, "dec", None) is not None

    return R.Tab("Deconvolution", [
        R.Group("Deconvolute", [_deconv_item(fr), _method_item(fr)]),
        R.Group("Results", [
            view_item(fr, "shifts", "Mass shifts…", "measure", lambda v: v.dec.shifts.open_dialog(), MS,
                      cond=has_dec),
            view_item(fr, "full_window", "Full window", "full", "open_full", ("ms",),
                      cond=lambda v: getattr(v, "FULL_WINDOW", False), tip="Open this spectrum in the Deconvolute window"),
            view_item(fr, "close_dec", "Close results", "clear", lambda v: v.dec.close(), MS, cond=has_dec)]),
        R.Group("Polymer", [
            view_item(fr, "poly%d" % i, lab, ic, (lambda v, i=i: __import__("ms_polymer_ui").open_tools(v, i)), HR)
            for i, (lab, ic) in enumerate([("Distribution…", "ms"), ("Repeat unit…", "measure"),
                                           ("Chain lengths…", "label")])])])


def report_tab(fr):
    return R.Tab("Report", [
        R.Group("Report", [_report_item(fr)]),
        R.Group("Export", [
            view_item(fr, "exp_peaks", "Peak table…", "table", "on_export_table", CHROM_VIEWS,
                      tip="Export the peak table"),
            view_item(fr, "exp_spec", "Spectrum…", "export", "on_export_spectrum", MS, tip="Export the spectrum"),
            view_item(fr, "exp_uv", "UV spectrum…", "export", "on_export_uv", PDA, tip="Export the UV spectrum"),
            view_item(fr, "exp_traces", "Traces (CSV)…", "compare", "on_export_traces", CMP),
            view_item(fr, "exp_table", "Table (CSV)…", "table", "on_export_table", CMP),
            view_item(fr, "exp_xlsx", "Areas (Excel)…", "export", "on_export_areas_xlsx", CMP)]),
        R.Group("Images", [R.Item("img_size", "Image size…", "image",
                                  lambda b: __import__("unilcms").ask_export_settings(fr),
                                  tip="Size of copied and saved images")])])


def view_tab(fr):
    def side_on():
        if kind(fr) == "compare":
            return bool(_safe(lambda: current_view(fr).side.IsShown()))
        return bool(fr.__dict__.get("_side_shown"))

    def files_toggle(b):
        fp = fr.files
        fp.collapsed = not fp.collapsed
        T._save({"files_collapsed": fp.collapsed})
        fp._apply_width()
        fp.GetParent().Layout()
        fp.Refresh()

    return R.Tab("View", [
        R.Group("Panels", [
            R.Item("side", "Side panel", "panel", lambda b: fr.set_side(None), tip="Show or hide the side panel (F9)",
                   active=side_on),
            R.Item("files", "File list", "list", files_toggle, tip="Narrow or widen the list of files",
                   active=lambda: not fr.files.collapsed)]),
        R.Group("Plots", [tool_item(fr, "full", "Fit"),
                          view_item(fr, "graph_props", "Graph…", "settings", "on_graph_props", CMP,
                                    tip="Graph properties")]),
        R.Group("Files", [
            R.Item("prev_doc", "Previous", "prev", lambda b: fr.next_doc(-1), size="small",
                   tip="Previous file (Ctrl+PgUp)", enabled=lambda: len(fr.docs) > 1),
            R.Item("next_doc", "Next", "next", lambda b: fr.next_doc(1), size="small", tip="Next file (Ctrl+PgDn)",
                   enabled=lambda: len(fr.docs) > 1),
            R.Item("close_doc", "Close", "clear", lambda b: fr.close_doc(fr.active), size="small",
                   tip="Close the file shown (Ctrl+W)", enabled=lambda: bool(getattr(fr.active, "path", None)))]),
        R.Group("Help", [
            R.Item("keys", "Shortcuts", "list", lambda b: __import__("shortcuts").show_window(fr),
                   tip="Every keyboard shortcut (F1)"),
            R.Item("about", "About", "help", lambda b: fr.on_help(), tip="About " + fr.TITLE)])])


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
    return [home_tab(fr), peaks_tab(fr), spectrum_tab(fr), deconv_tab(fr), report_tab(fr), view_tab(fr)]


# ------------------------------------------------------------------ in the window
def install(fr, root, tb, vs):
    """The ribbon in place of the top bar tb (which stays, hidden) in the sizer vs of the window."""
    rb = R.Ribbon(root, tabs(fr), view_group(fr))
    fr.__dict__["ribbon"] = rb
    row = rb.tabrow
    for key, chip in list(getattr(fr, "chips", {}).items()):  # the chips: at the right of the tab row
        try:
            tb.sizer.Detach(chip)
        except Exception:
            pass
        chip.Reparent(row)
        row.add(chip, 3)
    FB = T._cls("FlatButton")
    hb = FB(row, wx.ID_ANY, "", icon="help", kind="ghost", tooltip="About " + fr.TITLE, height=30, padx=8)
    hb.Bind(wx.EVT_BUTTON, fr.on_help)
    row.add(hb, 6)
    row.sizer.AddSpacer(fr.FromDIP(8))
    fr._top_labels = {}  # (no top bar to fit)
    tb.Hide()
    vs.Insert(0, rb, 0, wx.EXPAND)
    return rb

# -*- coding: utf-8 -*-
"""Graph properties of a plot (right click > Graph properties in the Compare view of LCMS Analysis):
fonts and text sizes, axis titles and ranges, tick steps and direction, frame, grid, the line of
each trace, legend, a title and a panel letter. The properties belong to the plot they were set on
and last for the session (they are not saved).

props(): the defaults, which give the plot exactly as it is drawn without them. apply_axes() puts
the text, tick, frame and grid properties on a matplotlib axes after the owner has drawn it;
the owner uses the others (ranges, trace lines, label sizes, legend) where it draws those.
GraphPropsDialog edits a copy and calls on_apply(props) for Apply and OK (Cancel puts the
properties of the moment the dialog opened back)."""
import copy

import wx

import unidec_theme as T
import unilcms as U
from unidec_theme import ui_font

FONTS = [("", "As the program"), ("Arial", "Arial"), ("Helvetica", "Helvetica"), ("Calibri", "Calibri"),
         ("Segoe UI", "Segoe UI"), ("Verdana", "Verdana"), ("Times New Roman", "Times New Roman"),
         ("Cambria", "Cambria"), ("Georgia", "Georgia"), ("Courier New", "Courier New")]
# a font missing on the computer: the next of its list (matplotlib takes the first it has)
FALLBACK = {"Arial": ["Liberation Sans", "DejaVu Sans"], "Helvetica": ["Arial", "Liberation Sans", "DejaVu Sans"],
            "Calibri": ["Carlito", "Arial", "DejaVu Sans"], "Segoe UI": ["Arial", "DejaVu Sans"],
            "Verdana": ["DejaVu Sans"], "Times New Roman": ["Liberation Serif", "DejaVu Serif"],
            "Cambria": ["Caladea", "Times New Roman", "DejaVu Serif"], "Georgia": ["Times New Roman", "DejaVu Serif"],
            "Courier New": ["Liberation Mono", "DejaVu Sans Mono"]}
FRAMES = [("box", "Box (four sides)"), ("open", "Left and bottom"), ("bottom", "Bottom only"), ("none", "None")]
TICKDIR = [("out", "Outside"), ("in", "Inside"), ("inout", "Across the frame")]
GRIDS = [("none", "None"), ("x", "Vertical lines"), ("y", "Horizontal lines"), ("both", "Both")]
STYLES = [("-", "Solid"), ("--", "Dashed"), (":", "Dotted"), ("-.", "Dash dot")]
LEGEND_LOC = [("best", "Best place"), ("upper right", "Top right"), ("upper left", "Top left"),
              ("lower right", "Bottom right"), ("lower left", "Bottom left"), ("center right", "Middle right")]
BLACK = "#000000"
MAX_TICKS = 50  # a tick step giving more ticks than this in the view is replaced by the automatic ticks


def props():
    """The defaults: the plot as drawn without graph properties (sizes in points)."""
    return {"font": "Arial",
            "size_ticks": 9.0, "size_titles": 10.0, "size_names": 8.0, "size_rt": 7.0, "size_guides": 7.5,
            "size_bar": 8.5, "bold_titles": False, "bold_names": False, "bold_rt": False,
            "xtitle": "", "ytitle": "", "x0": None, "x1": None, "y0": None, "y1": None, "xstep": None, "ystep": None,
            "minor": True, "tickdir": "out", "frame": "open", "frame_lw": 0.8, "grid": "none",
            "legend_loc": "best", "legend_frame": False, "fill_alpha": 16,
            "title": "", "size_title": 11.0, "bold_title": True, "letter": "", "size_letter": 12.0,
            "traces": {}}  # per trace (key of the owner): {"lw": pt or None (the common width), "ls": "-"}


def is_default(p):
    d = props()
    return all(p.get(k) == v for k, v in d.items())


_HAVE = {}


def _have(name):
    """Is the font installed? (matplotlib logged a warning for every missing font of a list, on every draw)"""
    if name not in _HAVE:
        try:
            from matplotlib import font_manager as fm
            fm.findfont(fm.FontProperties(family=name), fallback_to_default=False)
            _HAVE[name] = True
        except Exception:
            _HAVE[name] = False
    return _HAVE[name]


def family(p):
    """matplotlib font family list of the properties (the installed ones of the font and its stand-ins),
    None for the program's font."""
    f = p.get("font") or ""
    if not f:
        return None
    fam = [f] + FALLBACK.get(f, [])
    return [x for x in fam if _have(x)] or fam[-1:]


_FOUND = {}


def first_family(p):
    """The first font of family(p) that this computer has (tick numbers take one family only), None for the
    program's font."""
    fam = family(p)
    if not fam:
        return None
    key = tuple(fam)
    if key not in _FOUND:
        _FOUND[key] = fam[-1]
        try:
            from matplotlib import font_manager as fm
            for f in fam:
                try:
                    fm.findfont(fm.FontProperties(family=f), fallback_to_default=False)
                    _FOUND[key] = f
                    break
                except Exception:
                    continue
        except Exception:
            pass
    return _FOUND[key]


def ink(p, default):
    return default or BLACK


def font_kw(p, size_key, bold_key=None):
    """fontsize, fontweight and fontfamily keywords of a text of the plot."""
    kw = {"fontsize": float(p.get(size_key, 8.0))}
    if bold_key and p.get(bold_key):
        kw["fontweight"] = "bold"
    fam = family(p)
    if fam:
        kw["fontfamily"] = fam
    return kw


def _step_locator(step, auto):
    """Major ticks every step, unless that gives more than MAX_TICKS in the view (a tiny step drew
    thousands of ticks, which took minutes or ran out of memory) or none at all (a step larger than
    the view): then the automatic ticks of auto."""
    from matplotlib.ticker import MultipleLocator

    class StepLocator(MultipleLocator):
        def __call__(self):
            vmin, vmax = sorted(self.axis.get_view_interval())
            n = (vmax - vmin) / step
            if not n <= MAX_TICKS:
                return auto()
            locs = self.tick_values(vmin, vmax)
            if not any(vmin <= v <= vmax for v in locs):
                return auto()
            return locs

    return StepLocator(step)


def apply_axes(ax, p, default_ink, y_values=True):
    """Text, ticks, frame, grid, title and panel letter of the properties on the axes (after the owner
    drew it and set its axis titles). y_values: the y axis shows numbers (not a scale bar or nothing)."""
    from matplotlib.ticker import NullLocator, AutoMinorLocator
    col = ink(p, default_ink)
    fam = family(p)
    lw = float(p.get("frame_lw", 0.8))
    tick_kw = dict(which="major", labelsize=float(p["size_ticks"]), labelcolor=col)
    ax.tick_params(**tick_kw)
    ax.tick_params(which="both", direction=p.get("tickdir", "out"), color=col, width=lw)
    # the tick numbers keep a family set once (also through ax.cla()): set every time, the program's own
    # (rcParams) without a font chosen
    import matplotlib
    rc = matplotlib.rcParams["font.family"]
    tick_fam = first_family(p) if fam else (rc[0] if isinstance(rc, (list, tuple)) and rc else str(rc))
    try:
        ax.tick_params(which="both", labelfontfamily=tick_fam)
    except (TypeError, ValueError, AttributeError):  # an older matplotlib: the labels one by one
        for t in ax.get_xticklabels() + ax.get_yticklabels():
            t.set_fontfamily(fam or tick_fam)
    for t in (ax.xaxis.get_offset_text(), ax.yaxis.get_offset_text()):
        t.set_color(col)
        if fam:
            t.set_fontfamily(fam)
    for axis, key in ((ax.xaxis, "xtitle"), (ax.yaxis, "ytitle")):
        lab = axis.label
        if p.get(key):
            lab.set_text(p[key])
        lab.set_fontsize(float(p["size_titles"]))
        lab.set_fontweight("bold" if p.get("bold_titles") else "normal")
        lab.set_color(col)
        if fam:
            lab.set_fontfamily(fam)
    if not p.get("minor", True):
        ax.xaxis.set_minor_locator(NullLocator())
        ax.yaxis.set_minor_locator(NullLocator())
    if p.get("xstep"):
        ax.xaxis.set_major_locator(_step_locator(float(p["xstep"]), ax.xaxis.get_major_locator()))
        if p.get("minor", True):
            ax.xaxis.set_minor_locator(AutoMinorLocator())
    if p.get("ystep") and y_values:
        ax.yaxis.set_major_locator(_step_locator(float(p["ystep"]), ax.yaxis.get_major_locator()))
        if p.get("minor", True):
            ax.yaxis.set_minor_locator(AutoMinorLocator())
    fr = p.get("frame", "box")
    show = {"box": ("left", "right", "top", "bottom"), "open": ("left", "bottom"), "bottom": ("bottom",),
            "none": ()}.get(fr, ("left", "right", "top", "bottom"))
    for k, sp in ax.spines.items():
        sp.set_visible(k in show)
        sp.set_linewidth(lw)
        sp.set_color(col)
    if fr in ("bottom", "none"):
        ax.tick_params(axis="y", which="both", left=False)
    if fr == "none":
        ax.tick_params(axis="x", which="both", bottom=False)
    g = p.get("grid", "none")
    if g != "none":
        ax.grid(True, axis={"x": "x", "y": "y", "both": "both"}[g], which="major", color="#C9D0DA", lw=0.5,
                ls=(0, (2, 2)), zorder=0)
        ax.set_axisbelow(True)
    else:
        ax.grid(False)
    apply_titles(ax, p, default_ink)


def apply_titles(ax, p, default_ink, pad=6.0):
    """The title (centred) and the panel letter (bold, left) above the frame, pad points above it (more
    when other labels sit there, e.g. the times of guide lines)."""
    col = ink(p, default_ink)
    title_kw = font_kw(p, "size_title", "bold_title")
    ax.set_title(p.get("title") or "", loc="center", color=col, pad=pad, **title_kw)
    letter_kw = font_kw(p, "size_letter")
    letter_kw["fontweight"] = "bold"
    ax.set_title(p.get("letter") or "", loc="left", color=col, pad=pad, **letter_kw)


def fit_title(ax, r, fitted=False, min_size=6.0, max_lines=3):
    """A title wider than the frame (a long title, a narrow copied image) goes on more lines, broken
    at spaces, and gets smaller type where it would take more than max_lines (down to min_size points);
    centred, it also keeps clear of the panel letter at the left. It ran out of the figure and over the
    letter, and the copied image lost its layout. fitted: the frame is fitted to the texts around it
    afterwards (copied and saved images, the report: fit_axes). Call it after the texts left and right
    of the frame are placed."""
    t = ax.title
    txt = t.get_text()
    if not txt:
        return
    fig = ax.figure
    pt = fig.dpi / 72.0
    ab = ax.get_window_extent(r)
    gap = 2.0 * (ax._left_title.get_window_extent(r).width + 6.0 * pt) if ax._left_title.get_text() else 0.0
    if t.get_window_extent(r).width <= ab.width - gap:
        return
    room = ab.width
    if fitted:
        # the width the frame gets: the figure less the texts left and right of it, not the frame of the
        # moment, which the long title had narrowed already (before it was fitted)
        vis = [(a, a.get_visible()) for a in (ax.title, ax._left_title)]
        try:
            for a, _ in vis:
                a.set_visible(False)
            tb = ax.get_tightbbox(r)
            room = max(ab.width, fig.bbox.width - max(0.0, ab.x0 - tb.x0) - max(0.0, tb.x1 - ab.x1) - 6.0 * pt)
        except Exception:
            pass
        finally:
            for a, v in vis:
                a.set_visible(v)
    room = max(room - gap, 40.0 * pt)
    import textwrap
    while True:
        t.set_text(txt)
        lines = 1
        w = t.get_window_extent(r).width
        if w > room and " " in txt.strip():
            n = max(8, int(len(txt) * room / w))
            while True:
                parts = textwrap.wrap(txt, n, break_long_words=False)
                t.set_text("\n".join(parts))
                lines = len(parts)
                if t.get_window_extent(r).width <= room or n <= 8:
                    break
                n = max(8, int(n * 0.9))
        if (t.get_window_extent(r).width <= room and lines <= max_lines) or t.get_fontsize() <= min_size:
            break
        t.set_fontsize(max(min_size, t.get_fontsize() - 0.5))


def title_room(p, pad=6.0, ax=None, r=None):
    """Height (inches) above the frame taken by the title and the panel letter, 0 without them (measured
    when the axes and a renderer are given: a title on more lines)."""
    hs = [float(p["size_title"]) for _ in [0] if p.get("title")] + [float(p["size_letter"]) for _ in [0]
                                                                     if p.get("letter")]
    if hs and ax is not None and r is not None:
        try:
            pt = ax.figure.dpi / 72.0
            hs += [t.get_window_extent(r).height / pt / 1.3 for t in (ax.title, ax._left_title) if t.get_text()]
        except Exception:
            pass
    return (pad + 1.3 * max(hs)) / 72.0 + 0.04 if hs else 0.0


def _num(text, lo=None, hi=None):
    """A number typed (comma as decimal point too), None when empty or not a number."""
    t = (text or "").strip().replace(",", ".")
    if not t:
        return None
    try:
        v = float(t)
    except ValueError:
        return None
    if v != v or v in (float("inf"), float("-inf")):  # "nan" or "inf" typed: not a number to use
        return None
    if lo is not None:
        v = max(lo, v)
    if hi is not None:
        v = min(hi, v)
    return v


def _fmt(v):
    return "" if v is None else "%g" % v


def _sel(choices, key):
    for i, (k, _) in enumerate(choices):
        if k == key:
            return i
    return 0


class GraphPropsDialog(wx.Dialog):
    """Graph properties of one plot. traces: [(key, name, colour, width or None, style)] of the
    traces drawn; common_lw: the line width of every trace (side panel). on_apply(props, colours,
    common_lw) draws the plot with them (colours: {key: colour} changed in the dialog; None: the colours
    of the moment the dialog opened, for Cancel)."""

    PAGES = ("Text", "Axes", "Traces", "Title")

    def __init__(self, parent, p, traces, common_lw, on_apply, auto_xtitle="", auto_ytitle=""):
        wx.Dialog.__init__(self, wx.GetTopLevelParent(parent), title="Graph properties",
                           style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        self.SetBackgroundColour(wx.Colour(U.C["panel"]))
        self.orig = copy.deepcopy(p)
        self.orig_colours = {k: c for k, _, c, _, _ in traces}
        self.cur_colours = dict(self.orig_colours)
        self.orig_lw = common_lw
        self.p = copy.deepcopy(p)
        self.traces = traces
        self.on_apply = on_apply
        self.applied = False
        nb = self.nb = wx.Notebook(self)
        nb.SetFont(ui_font(9, 400))
        self._text_page(nb)
        self._axes_page(nb, auto_xtitle, auto_ytitle)
        self._traces_page(nb, common_lw)
        self._title_page(nb)
        vs = wx.BoxSizer(wx.VERTICAL)
        vs.Add(nb, 1, wx.EXPAND | wx.ALL, self.FromDIP(8))
        reset = U.flat(self, "Reset all", handler=self.on_reset)
        vs.Add(U.dialog_buttons(self, U.flat(self, "Cancel", handler=self.on_cancel),
                                U.flat(self, "Apply", handler=lambda e: self.apply()),
                                U.flat(self, "OK", "primary", handler=self.on_ok), left=(reset,)), 0, wx.EXPAND)
        self.SetSizerAndFit(vs)
        try:  # never taller than the screen (a small laptop screen at 125 %): OK and Cancel stay in reach
            area = T.display_area(self)
            if self.GetSize().height > area.height:
                self.SetMinSize(wx.Size(self.GetSize().width, area.height))  # (the fit set it as the minimum)
                self.SetSize(wx.Size(self.GetSize().width, area.height))
        except Exception:
            pass
        self.SetMinSize(self.GetSize())
        self.Bind(wx.EVT_CLOSE, self.on_cancel)
        self.Bind(wx.EVT_CHAR_HOOK, self._on_key)
        # Enter in a field: OK, as in the other dialogs of the program
        for c in [spec[1] for spec in self._ctl.values()] + [w for _, _, w, _ in self._trace_ctl]:
            if isinstance(c, wx.TextCtrl):
                c.Bind(wx.EVT_TEXT_ENTER, self.on_ok)
        self.CentreOnParent()

    def _on_key(self, e):
        if e.GetKeyCode() == wx.WXK_ESCAPE:  # Esc: Cancel (the dialog has no standard Cancel button)
            self.on_cancel()
        else:
            e.Skip()

    # ------------------------------------------------------------------ pages
    def _page(self, nb, title):
        fm = U.Form(nb, label_w=150)
        nb.AddPage(fm, title)
        return fm

    def _size(self, fm, label, key, bold_key=None, tip=None):
        c = fm.text(_fmt(self.p[key]), 46, tooltip=tip)
        items = [c, "pt"]
        b = None
        if bold_key:
            b = wx.CheckBox(fm, label="bold")
            b.SetFont(ui_font(9, 400))
            b.SetValue(bool(self.p.get(bold_key)))
            items.append(b)
        fm.row(label, *items)
        self._ctl[key] = ("num", c, 4.0, 40.0)
        if bold_key:
            self._ctl[bold_key] = ("check", b)

    def _text_page(self, nb):
        self._ctl = {}
        fm = self._page(nb, "Text")
        fm.section("Font", T.GROUP["blue"])
        c = fm.choice([v for _, v in FONTS], 190)
        c.SetSelection(_sel(FONTS, self.p["font"]))
        fm.row("Font", c)
        self._ctl["font"] = ("choice", c, FONTS)
        fm.section("Sizes", T.GROUP["green"])
        self._size(fm, "Axis numbers", "size_ticks")
        self._size(fm, "Axis titles", "size_titles", "bold_titles")
        self._size(fm, "Trace names", "size_names", "bold_names", tip="Also the legend")
        self._size(fm, "Retention times", "size_rt", "bold_rt")
        self._size(fm, "Guide line times", "size_guides")
        self._size(fm, "Scale bar", "size_bar")

    def _axes_page(self, nb, auto_xtitle, auto_ytitle):
        fm = self._page(nb, "Axes")
        fm.section("Titles", T.GROUP["blue"])
        c = fm.text(self.p["xtitle"], 190, tooltip="Empty: %s" % (auto_xtitle or "automatic"))
        fm.row("x axis title", c)
        self._ctl["xtitle"] = ("text", c)
        c = fm.text(self.p["ytitle"], 190, tooltip="Empty: %s" % (auto_ytitle or "automatic"))
        fm.row("y axis title", c)
        self._ctl["ytitle"] = ("text", c)
        fm.section("Ranges and ticks", T.GROUP["green"])
        a, b = fm.text(_fmt(self.p["x0"]), 58, "Empty: automatic"), fm.text(_fmt(self.p["x1"]), 58, "Empty: automatic")
        fm.row("x range", a, "to", b)
        self._ctl["x0"], self._ctl["x1"] = ("num", a, None, None), ("num", b, None, None)
        a, b = fm.text(_fmt(self.p["y0"]), 58, "Empty: automatic"), fm.text(_fmt(self.p["y1"]), 58, "Empty: automatic")
        fm.row("y range", a, "to", b)
        self._ctl["y0"], self._ctl["y1"] = ("num", a, None, None), ("num", b, None, None)
        a, b = fm.text(_fmt(self.p["xstep"]), 58, "Empty: automatic"), fm.text(_fmt(self.p["ystep"]), 58,
                                                                                "Empty: automatic")
        fm.row("Tick step x, y", a, ",", b)
        self._ctl["xstep"], self._ctl["ystep"] = ("num", a, 1e-9, None), ("num", b, 1e-9, None)
        c = fm.choice([v for _, v in TICKDIR], 150)
        c.SetSelection(_sel(TICKDIR, self.p["tickdir"]))
        fm.row("Ticks", c)
        self._ctl["tickdir"] = ("choice", c, TICKDIR)
        self._ctl["minor"] = ("check", fm.check("Minor ticks", self.p["minor"]))
        fm.section("Frame", T.GROUP["red"])
        c = fm.choice([v for _, v in FRAMES], 150)
        c.SetSelection(_sel(FRAMES, self.p["frame"]))
        fm.row("Frame", c)
        self._ctl["frame"] = ("choice", c, FRAMES)
        c = fm.text(_fmt(self.p["frame_lw"]), 46)
        fm.row("Frame line", c, "pt")
        self._ctl["frame_lw"] = ("num", c, 0.1, 4.0)
        c = fm.choice([v for _, v in GRIDS], 150)
        c.SetSelection(_sel(GRIDS, self.p["grid"]))
        fm.row("Grid", c)
        self._ctl["grid"] = ("choice", c, GRIDS)

    def _traces_page(self, nb, common_lw):
        sw = wx.ScrolledWindow(nb, style=wx.VSCROLL)
        sw.SetBackgroundColour(wx.Colour(U.C["panel"]))
        sw.SetScrollRate(0, self.FromDIP(10))
        fm = U.Form(sw, label_w=150)
        s = wx.BoxSizer(wx.VERTICAL)
        s.Add(fm, 1, wx.EXPAND)
        sw.SetSizer(s)
        nb.AddPage(sw, "Traces")
        fm.section("Every trace", T.GROUP["blue"])
        c = fm.text(_fmt(common_lw), 46, tooltip="Traces without their own width")
        fm.row("Line width", c, "pt")
        self._ctl["_common_lw"] = ("num", c, 0.2, 4.0)
        c = fm.text(_fmt(self.p["fill_alpha"]), 46)
        fm.row("Fill opacity", c, "%")
        self._ctl["fill_alpha"] = ("num", c, 0.0, 100.0)
        c = fm.choice([v for _, v in LEGEND_LOC], 150)
        c.SetSelection(_sel(LEGEND_LOC, self.p["legend_loc"]))
        fm.row("Legend", c)
        self._ctl["legend_loc"] = ("choice", c, LEGEND_LOC)
        self._ctl["legend_frame"] = ("check", fm.check("Frame around the legend", self.p["legend_frame"]))
        fm.section("Each trace", T.GROUP["green"])
        self._trace_ctl = []
        if not self.traces:
            fm.note("No traces")
        for key, name, colour, lw, ls in self.traces:
            pk = wx.ColourPickerCtrl(fm, colour=wx.Colour(colour))
            pk.SetToolTip("Colour of %s" % name)
            w = fm.text(_fmt(lw), 40, tooltip="Empty: the common width")
            st = fm.choice([v for _, v in STYLES], 92)
            st.SetSelection(_sel(STYLES, ls or "-"))
            short = name if len(name) <= 24 else name[:22] + "…"
            fm.row(short, pk, w, "pt", st)
            self._trace_ctl.append((key, pk, w, st))
        fm.Fit()
        sw.FitInside()
        rows = min(len(self.traces), 7)
        sw.SetMinSize(wx.Size(-1, self.FromDIP(250 + 34 * rows)))

    def _title_page(self, nb):
        fm = self._page(nb, "Title")
        fm.section("Title", T.GROUP["blue"])
        c = fm.text(self.p["title"], 220)
        fm.row("Title", c)
        self._ctl["title"] = ("text", c)
        self._size(fm, "Title size", "size_title", "bold_title")
        fm.section("Panel letter", T.GROUP["green"])
        c = fm.text(self.p["letter"], 46, tooltip="E.g. A or (b)")
        fm.row("Letter", c)
        self._ctl["letter"] = ("text", c)
        self._size(fm, "Letter size", "size_letter")

    # ------------------------------------------------------------------ values
    def read(self):
        """(props, colours, common width) of the controls."""
        p = copy.deepcopy(self.p)
        common = self.orig_lw
        for k, spec in self._ctl.items():
            kind, c = spec[0], spec[1]
            if kind == "num":
                v = _num(c.GetValue(), spec[2], spec[3])
                if k == "_common_lw":
                    common = v if v is not None else self.orig_lw
                elif k.startswith("size_") or k in ("frame_lw", "fill_alpha"):
                    p[k] = v if v is not None else props()[k]
                else:
                    p[k] = v
            elif kind == "text":
                p[k] = c.GetValue().strip()
            elif kind == "check":
                p[k] = bool(c.GetValue())
            elif kind == "choice":
                p[k] = spec[2][max(0, c.GetSelection())][0]
        for a, b in (("x0", "x1"), ("y0", "y1")):
            if p[a] is not None and p[b] is not None and p[b] < p[a]:
                p[a], p[b] = p[b], p[a]
            if p[a] is not None and p[b] is not None and p[a] == p[b]:
                p[a] = p[b] = None  # an empty range: automatic
        colours = {}
        tr = {}
        for key, pk, w, st in self._trace_ctl:
            colours[key] = pk.GetColour().GetAsString(wx.C2S_HTML_SYNTAX)
            lw = _num(w.GetValue(), 0.2, 6.0)
            ls = STYLES[max(0, st.GetSelection())][0]
            if lw is not None or ls != "-":
                tr[key] = {"lw": lw, "ls": ls}
        # traces not drawn now (left out or another signal) keep their lines
        keep = {k: v for k, v in self.p.get("traces", {}).items() if k not in [t[0] for t in self.traces]}
        keep.update(tr)
        p["traces"] = keep
        return p, colours, common

    def apply(self):
        p, colours, common = self.read()
        changed = {k: c for k, c in colours.items() if c.lower() != (self.cur_colours.get(k) or "").lower()}
        self.on_apply(p, changed, common)
        self.cur_colours.update(changed)
        self.applied = True
        return p

    def on_ok(self, e=None):
        self.apply()
        self.EndModal(wx.ID_OK)

    def on_cancel(self, e=None):
        if self.applied:  # the plot as it was when the dialog opened
            self.on_apply(copy.deepcopy(self.orig), None, self.orig_lw)
        self.EndModal(wx.ID_CANCEL)

    def on_reset(self, e=None):
        d = props()
        d["traces"] = {}
        for k, spec in self._ctl.items():
            kind, c = spec[0], spec[1]
            if k == "_common_lw":
                continue
            if kind in ("num", "text"):
                c.SetValue(_fmt(d[k]) if kind == "num" else d[k])
            elif kind == "check":
                c.SetValue(bool(d[k]))
            elif kind == "choice":
                c.SetSelection(_sel(spec[2], d[k]))
        for key, pk, w, st in self._trace_ctl:
            w.SetValue("")
            st.SetSelection(0)
            if self.orig_colours.get(key):  # colours changed in this dialog: back to those it opened with
                pk.SetColour(wx.Colour(self.orig_colours[key]))
        self.p["traces"] = {}

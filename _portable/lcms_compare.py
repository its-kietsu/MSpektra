# -*- coding: utf-8 -*-
"""Compare view of LCMS Postrun: the chromatograms of several open files stacked on top of each
other (or overlaid), each processed the same way: PDA at a wavelength, PDA max plot, TIC, base
peak or mass chromatogram; smoothing, baseline, time window, alignment on a peak, scale. The view
belongs to the window (not to one file): every open file is listed and can be ticked, ordered,
renamed and coloured. The numbers come from lcms_compare_core (no wx)."""
import os
import datetime
import time

import numpy as np
import wx

import unidec_theme as T
from unidec_theme import C, ui_font
import unilcms as U
import lcms_compare_core as K
import lcms_integrate as LI
import plot_props as PP

SETTINGS_KEY = "compare"
U.TOOL_ICONS.setdefault("align", '<path d="M12 3v18M5 8l3 4-3 4M19 8l-3 4 3 4"/>')
U.TOOL_ICONS.setdefault("guide", '<path d="M12 3v3M12 9v3M12 15v3M12 21v0"/><path d="M4 20c3 0 4-9 6-9s3 9 5 9 3-4 5-4"/>')
U.TOOL_ICONS.setdefault("mzspec", '<path d="M3 20h18M6 20v-6M10 20V5M14 20v-9M18 20v-4"/>')
U.TOOL_ICONS.setdefault("area", U.TOOL_ICONS.get("integrate", U.TOOL_ICONS["drag_peak"]))
SPEC_AVG = [("peak", "Over the peak"), ("scan", "One scan at the peak top")]
GRADIENT = ("#0B2A6B", "#8DB8F2")  # dark to light: a series (time points, fractions)


class ComparePage(wx.Panel):
    """Place of the Compare view in the views of one file: the view itself is shared by the window
    and shown instead of the file's views while this page is chosen."""
    is_compare = True

    def __init__(self, parent, frame):
        wx.Panel.__init__(self, parent, style=wx.BORDER_NONE)
        self.SetBackgroundColour(C["bg"])
        self.frame = frame


def _mix(a, b, f):
    a, b = a.lstrip("#"), b.lstrip("#")
    ca = [int(a[i:i + 2], 16) for i in (0, 2, 4)]
    cb = [int(b[i:i + 2], 16) for i in (0, 2, 4)]
    return "#%02X%02X%02X" % tuple(int(round(x + (y - x) * f)) for x, y in zip(ca, cb))


def _tex(text):
    """Text for matplotlib: a $ in a sample name would start a formula."""
    return str(text).replace("$", r"\$")


def _ink(col):
    """Colour of the text that belongs to a trace drawn in col: light colours (the light end of the
    gradient) are darkened, so that a label stays readable on white."""
    try:
        c = col.lstrip("#")
        r, g, b = [int(c[i:i + 2], 16) / 255.0 for i in (0, 2, 4)]
    except (AttributeError, ValueError, IndexError):
        return col
    return _mix(col, "#1B2330", 0.45) if 0.2126 * r + 0.7152 * g + 0.0722 * b > 0.55 else col


def _spread(cx, w, gap, lo, hi):
    """Centres for labels of widths w wanted at the centres cx (sorted, display units) so that they
    do not overlap: touching labels form a group centred on their wanted places, kept within lo to
    hi."""
    groups = [[i] for i in range(len(cx))]

    def place(g):
        tot = sum(w[i] for i in g) + gap * (len(g) - 1)
        start = sum(cx[i] for i in g) / len(g) - tot / 2.0
        start = max(lo, min(start, hi - tot)) if hi - lo > tot else lo
        out, x = [], start
        for i in g:
            out.append(x + w[i] / 2.0)
            x += w[i] + gap
        return start, start + tot, out

    merged = True
    while merged and len(groups) > 1:
        merged = False
        for k in range(len(groups) - 1):
            if place(groups[k])[1] + gap > place(groups[k + 1])[0]:
                groups[k:k + 2] = [groups[k] + groups[k + 1]]
                merged = True
                break
    pos = [0.0] * len(cx)
    for g in groups:
        for i, p in zip(g, place(g)[2]):
            pos[i] = p
    return pos


class _MzRow(U.SplitBox):
    """Row of the two spectrum tiles of the m/z tool, below a line with the suggested neutral mass (bold)
    and a note (muted), each shortened to the width (the tooltip of the line has all)."""

    def __init__(self, parent):
        U.SplitBox.__init__(self, parent, wx.HORIZONTAL)
        self.head_h = self.FromDIP(26)
        self.head = ("", "")
        self.title = "Mass spectra (m/z)"

    def set_head(self, bold, rest):
        if (bold, rest) != self.head:
            self.head = (bold, rest)
            try:
                self.SetToolTip((bold + ("   " + rest if rest else "")).strip())
            except RuntimeError:
                pass
            self.RefreshRect(wx.Rect(0, 0, self.GetClientSize()[0], self.head_h))

    def relayout(self):  # the tiles below the line
        w, h = self.GetClientSize()
        hh = getattr(self, "head_h", 0)
        pos = 0
        for win, size in zip(self.panes, self._sizes()):
            win.SetSize(pos, hh, max(size, 1), max(h - hh, 1))
            pos += size + self.gap
        self.Refresh()

    def _paint(self, e):
        dc = wx.PaintDC(self)
        if not hasattr(self, "head_h"):
            return
        gc = T.crisp(dc)
        w, h = self.GetClientSize()
        x = self.FromDIP(14)
        room = w - x - self.FromDIP(10)
        for text, font, col in ((self.head[0], ui_font(9.5, 700), C["accent_text"]),
                                (self.head[1], ui_font(9, 400), C["muted"])):
            if not text or room < self.FromDIP(30):
                continue
            gc.SetFont(font, wx.Colour(col))
            t = text
            if gc.GetTextExtent(t)[0] > room:
                while t and gc.GetTextExtent(t + "…")[0] > room:
                    t = t[:-1]
                t = (t.rstrip() + "…") if t else ""
            tw, th = gc.GetTextExtent(t)
            gc.DrawText(t, x, (self.head_h - th) / 2.0 - self.FromDIP(1))
            x += tw + self.FromDIP(12)
            room -= tw + self.FromDIP(12)
        if self._hover is None or self._drag is not None:  # the grip of the gap under the mouse (SplitBox)
            return
        sizes = self._sizes()
        pos = sum(sizes[:self._hover + 1]) + self.gap * self._hover
        L, t_ = self.FromDIP(40), max(2, self.FromDIP(4))
        dc.SetPen(wx.TRANSPARENT_PEN)
        dc.SetBrush(wx.Brush(wx.Colour(C["line2"])))
        dc.DrawRoundedRectangle(pos + (self.gap - t_) // 2, self.head_h + (h - self.head_h - L) // 2, t_, L, t_ // 2)


def _key(choices, label):
    return [k for k, v in choices][label] if 0 <= label < len(choices) else choices[0][0]


def _index(choices, key):
    for i, (k, v) in enumerate(choices):
        if k == key:
            return i
    return 0


def _fnum(text, default=None):
    try:
        t = str(text).strip().replace(",", ".")
        return float(t) if t else default
    except ValueError:
        return default


class CompareTab(U.TabBase):
    ylabel = "Absorbance (mAU)"
    TOOL_HINTS = dict(U.TabBase.TOOL_HINTS)
    TOOL_HINTS.update({
        "align": "Align: click a peak to line up every trace there",
        "guide": "Guide line: click to add or remove a dashed line",
        "area": "Area: drag over a peak to integrate every trace",
        "mzspec": "m/z: click a peak of a trace for its mass spectra",
        # the spectrum tools of the Mass spectrometry view, on the spectrum tiles of the m/z tool (shown while it is on)
        "xic": "Mass chromatogram: click an ion in a spectrum below",
        "measure": "Measure: click two ions in a spectrum below",
        "label": "Label: click an ion to pin or unpin its label"})
    SPEC_TOOLS = ("xic", "measure", "label")
    PIN_FMT, PICK_MIN, MZ_FMT = U.MSTab.PIN_FMT, U.MSTab.PICK_MIN, U.MSTab.MZ_FMT  # as the Mass spectrometry view

    def __init__(self, parent, frame):
        U.TabBase.__init__(self, parent, frame)
        self.s = K.default_settings()
        # the spectra of the m/z tool: averaged over the peak or one scan, the spectrum before the peak
        # subtracted, the number of ions labelled
        # the background before the peak subtracted by default, as in the MS view (solvent and TFA ions,
        # present all through the run, otherwise make most of the spectrum and of the suggested mass)
        self.s.update({"spec_avg": "peak", "spec_bg": True, "spec_labels": 8})
        saved = T._load().get(SETTINGS_KEY)
        if isinstance(saved, dict):
            saved = K.migrate(saved)
            for k in ("signal", "wl", "bw", "own_wl", "polarity", "mz", "mz_win", "mz_ppm", "smooth", "baseline",
                      "rolling_min", "align_win", "layout", "scale", "spacing", "skew", "reverse", "colours", "lw",
                      "fill", "labels", "label_text", "rt_labels", "rt_min", "yaxis", "off_spacing", "off_skew",
                      "integ_thr", "spec_avg", "spec_bg", "spec_labels"):
                if k in saved:
                    self.s[k] = saved[k]
            if self.s["layout"] not in [k for k, _ in K.LAYOUTS]:
                self.s["layout"] = "stacked"
            if self.s["spec_avg"] not in [k for k, _ in SPEC_AVG]:
                self.s["spec_avg"] = "peak"
            try:
                self.s["spec_labels"] = min(max(int(self.s["spec_labels"]), 0), 30)
            except (TypeError, ValueError):
                self.s["spec_labels"] = 8
        self.props = PP.props()  # graph properties of this plot (right click > Graph properties; not saved)
        self.entries = []  # {"doc", "include", "label" (typed or None), "colour" (or None), "shift" (min)}
        self.blank_doc = None  # a blank run subtracted from the others
        self.sel = None  # index of the entry chosen in the list
        self.drawn = []  # the traces drawn: entry, proc, offset, colour, label, units
        self._label_arts = []  # names, retention times and guide times placed by _place_labels
        self._table = None
        self._peak_table = None
        self._area_table = None
        self._cache = {}
        self._busy = False
        self.area_range = None  # displayed time range integrated in every trace
        self.area_results = []
        self._area_frames = []
        self._area_ref_key = None  # file chosen as the 100 % reference (kept while its trace is left out)
        self._area_x_error = None  # X values that could not be used (the areas stay, X empty)
        self._area_x_used = ""  # X values and X axis label of the last areas (the undo snapshot)
        self._area_xlabel_used = "Reaction time"
        self._wl_from_view = False  # the first wavelength comes from the PDA view (files_changed)
        self.scroller = U.RowScroller(self)
        self.rows = U.StackBox(self.scroller)
        self.rows.shape = None  # the plot fills the window, the table below it
        self.card = U.PlotCard(self.rows, "Comparison", "", mode="range")
        self.card.on_range = self.on_compare_range
        self.card.on_menu = self.plot_menu
        self.card.on_pick = self.on_pick
        self.card.on_view = lambda: self.replot(keep_view=True, quick=True)
        self.card.readout_fmt = self.readout
        self.card.tool_getter = lambda: self.tool
        self.card.tools_supported = {"align", "guide", "area"}
        self.card.on_tool = self.on_tool
        self.card.image_name = self.image_name
        self.card.margins = [0.84, 0.52, 0.2, 0.14]
        self.card.relabel_fn = self._place_labels  # copied and saved images: labels placed for their size
        self.card.canvas.mpl_connect("resize_event", lambda e: wx.CallAfter(self._relabel_resized))
        import deconv_tab as D
        self.table = D.TableCard(self.rows, "Peaks of each trace", [("Trace", 160)], on_select=self.on_row,
                                 empty="No peaks")
        self.table.extra_menu = lambda: [("Export the table (CSV)…", self.on_export_table)]
        # m/z tool (a toggle in the tool bar): the ESI+ and ESI- spectra of a peak clicked on a trace, in a row of
        # two tiles between the plot and the table, shown only while the tool is on
        self.mz_on = False
        self.mzpick = None  # the peak picked (pick_peak): run, times (MS), spectra, suggested masses
        self._mz_weight = 3.6
        self.mzrow = _MzRow(self.rows)
        self.mzcards = {}
        # the spectrum of each tile as the Mass spectrometry view keeps its spectra (unilcms helpers): pinned labels,
        # measurement
        self.mzcols = {}
        self._spec_style = None  # Display (right click a tile): 0 profile, 1 sticks; None: as the MS view of the file
        for pol in ("+", "-"):
            c = U.PlotCard(self.mzrow, self._pol_name(pol), "", mode="zoom")
            c.on_menu = lambda x, y, p=pol: self.spec_menu(p, x, y)
            c.on_view = lambda p=pol: self.plot_mz(p, keep_view=True)
            c.readout_fmt = lambda x, y: U.MSTab.READ_FMT % (x, U._g(y))
            c.image_name = lambda p=pol: self._spec_name(p)
            c.tool_getter = lambda: self.tool
            c.tools_supported = set(self.SPEC_TOOLS)
            c.on_tool = lambda tool, kind, x, y, p=pol: self.spec_tool(p, tool, x, y)
            c.relabel_fn = lambda p=pol: self._unclash_labels(p)  # copied and saved images: labels for their size
            c.canvas.mpl_connect("resize_event", lambda e, p=pol: wx.CallAfter(self._relabel_tile, p))
            c.Bind(wx.EVT_SIZE, lambda e, p=pol: (e.Skip(), self._fit_mz_title(p)))
            self.mzcards[pol] = c
            self.mzcols[pol] = {"e": None, "spec": None, "spec_card": c, "pins": [], "measure": None, "m1": None,
                                "desc": "", "range": None}
        self._mz_press = None  # (x, y, key) of the last left button press on the plot: the trace a drag starts on
        self.card.canvas.mpl_connect("button_press_event", self._on_press)
        self.mzrow.set_panes([self.mzcards["+"], self.mzcards["-"]], [1, 1])
        self.mzrow.Hide()
        self._build_side()
        self.rows.set_panes([self.card, self.table], [5, 1.5])
        self.scroller.set_box(self.rows)
        self.build_tools([])
        self._show_spec_tools()
        hs = wx.BoxSizer(wx.HORIZONTAL)
        left = wx.BoxSizer(wx.VERTICAL)
        left.Add(self.tools, 0, wx.EXPAND | wx.BOTTOM, self.FromDIP(10))
        left.Add(self.scroller, 1, wx.EXPAND)
        hs.Add(left, 1, wx.EXPAND | wx.ALL, self.FromDIP(12))
        hs.Add(self.side, 0, wx.EXPAND)
        self.SetSizer(hs)
        self.side.Show(bool(getattr(frame, "_side_shown", True)))
        self._tips()
        self._fill_controls()
        self.replot(keep_view=False)

    # ------------------------------------------------------------------ side panel
    def tool_items(self, spectrum_tools):
        return [("mode", "select", "", "select", "Select (Esc)"),
                ("mode", "zoom", "", "zoom", "Zoom: drag a box"),
                ("mode", "pan", "", "pan", "Pan: drag to move the view"),
                ("action", "full", "", "full", "Full view"),
                ("sep", ""),
                ("mode", "align", "Align", "align", self.TOOL_HINTS["align"]),
                ("mode", "guide", "Guide line", "guide", self.TOOL_HINTS["guide"]),
                ("mode", "area", "Area", "area", self.TOOL_HINTS["area"]),
                ("action", "mzspec", "m/z", "mzspec", self.TOOL_HINTS["mzspec"]),
                # the spectrum tools of the Mass spectrometry view, for the spectrum tiles (shown while m/z is on)
                ("sep", ""),
                ("mode", "xic", "Mass chrom.", "xic", self.TOOL_HINTS["xic"]),
                ("mode", "measure", "Measure", "measure", self.TOOL_HINTS["measure"]),
                ("mode", "label", "Label", "label", self.TOOL_HINTS["label"])]

    def tool_action(self, key):
        if key == "mzspec":  # a toggle, not a mouse mode: the tool chosen stays
            self.set_mz(not self.mz_on)
            return
        U.TabBase.tool_action(self, key)

    def set_tool(self, key):
        U.TabBase.set_tool(self, key)
        if key == "select":  # the hint of this view (the one of the file views speaks of chromatograms)
            self.status("Select: click a peak for its spectra, Ctrl + drag to zoom" if self.mz_on else
                        "Select: drag to integrate every trace, Ctrl + drag to zoom")

    def _show_spec_tools(self):
        """The spectrum tools (and the line before them) in the tool bar while m/z is on; a spectrum tool chosen
        when m/z goes off gives way to Select."""
        tools = getattr(self, "tools", None)
        if tools is None:
            return
        on = bool(self.mz_on)
        try:
            kids = list(tools.GetChildren())
            first = tools.buttons.get(self.SPEC_TOOLS[0])
            if first in kids and kids.index(first) > 0:
                kids[kids.index(first) - 1].Show(on)  # the separator
            for k in self.SPEC_TOOLS:
                b = tools.buttons.get(k)
                if b is not None:
                    b.Show(on)
            sel = tools.buttons.get("select")
            if sel is not None:
                sel.SetToolTip("Select: click a peak for its spectra (Esc)" if on else
                               "Select: drag to integrate every trace (Esc)")
            if not on and self.tool in self.SPEC_TOOLS:
                tools.set_mode("select")
                self.set_tool("select")
            import shortcuts
            shortcuts.decorate_tools(self)  # (the tooltip of Select was set again above)
            tools.Layout()
            tools.Refresh()
            self.Layout()
        except RuntimeError:
            pass  # closed meanwhile

    def all_cards(self):
        return [self.card] + ([self.mzcards["+"], self.mzcards["-"]] if self.mz_on else [])

    def _build_side(self):
        sp = self.side
        sp.section("Files", T.GROUP["blue"])
        self.flist = wx.CheckListBox(sp, size=wx.Size(-1, self.FromDIP(150)))
        # File names must not enlarge the virtual width of the whole sidebar.
        self.flist.SetMinSize(wx.Size(10, self.FromDIP(150)))
        self.flist.SetFont(ui_font(9, 400))
        self.flist.SetToolTip("Tick the files to compare")
        sp.full(self.flist)
        sp.buttons(U.flat(sp, "Up", handler=lambda e: self.move(-1)), U.flat(sp, "Down", handler=lambda e: self.move(1)),
                   U.flat(sp, "Reverse", handler=lambda e: self.reverse_order()))
        self.sort = sp.choice(["Sort…", "by name", "by acquisition time", "by sample name", "as opened"])
        self.sort.SetMinSize(wx.Size(self.FromDIP(150), -1))
        sp.row("Order", self.sort)
        self.f_label = sp.text("", 150, tooltip="Name in the plot (empty: automatic)")
        sp.row("Label", self.f_label)
        self.f_colour = wx.ColourPickerCtrl(sp, colour=wx.Colour(U.PAL[0]))
        self.f_colour_reset = U.flat(sp, "Auto", handler=lambda e: self.set_entry_colour(None))
        sp.row("Colour", self.f_colour, self.f_colour_reset)
        self.f_shift = sp.text("0", 64, tooltip="Added to the alignment")
        sp.row("Time shift", self.f_shift, "min")
        self.blank = sp.choice(["none"])
        self.blank.SetMinSize(wx.Size(self.FromDIP(150), -1))
        sp.row("Subtract blank", self.blank)

        sp.section("Signal", T.GROUP["yellow"])
        self.signal = sp.choice([v for k, v in K.SIGNALS])
        self.signal.SetMinSize(wx.Size(self.FromDIP(180), -1))
        sp.row("Trace", self.signal)
        self.wl = sp.text("", 64)
        self.wl_pick = sp.choice(["common…"] + ["%d nm" % w for w in K.COMMON_WL])
        self.wl_pick.SetMinSize(wx.Size(self.FromDIP(96), -1))
        sp.row("Wavelength", self.wl, "nm", self.wl_pick)
        self.bw = sp.text("4", 64, tooltip="Absorbance averaged over this range")
        sp.row("Bandwidth", self.bw, "nm")
        self.own_wl = sp.check("Each file at the wavelength of its PDA view", False)
        sp.buttons(U.flat(sp, "λmax of the main peak", handler=self.on_lambda_max,
                          tooltip="Of the tallest peak of the first file"))
        self.pol = sp.choice(["Positive ions (+)", "Negative ions (−)"])
        self.pol.SetMinSize(wx.Size(self.FromDIP(150), -1))
        sp.row("MS scans", self.pol)
        self.mz = sp.text("", 72)
        self.mz_win = sp.text("0.5", 52)
        self.mz_ppm = wx.CheckBox(sp, label="ppm")
        self.mz_ppm.SetFont(ui_font(9, 400))
        sp.row("m/z", self.mz, "±", self.mz_win, self.mz_ppm)

        sp.section("Processing", T.GROUP["green"])
        self.smooth = sp.spin(0, 0, 51)
        sp.row("Smoothing", self.smooth, "points")
        self.base = sp.choice([v for k, v in K.BASELINES])
        self.base.SetMinSize(wx.Size(self.FromDIP(180), -1))
        sp.row("Baseline", self.base)
        self.roll = sp.text("1.0", 52, tooltip="Wider than the peaks")
        sp.row("Rolling window", self.roll, "min")
        self.t0, self.t1 = sp.text("", 58), sp.text("", 58)
        sp.row("Time window", self.t0, "to", self.t1, "min")
        self.align = sp.text("", 58, tooltip="Every trace lines up its peak here")
        self.align_win = sp.text("0.3", 46)
        sp.row("Align on peak at", self.align, "±", self.align_win, "min")

        sp.section("Display", T.GROUP["red"])
        self.layout_c = sp.choice([v for k, v in K.LAYOUTS])
        sp.row("Layout", self.layout_c)
        self.scale = sp.choice([v for k, v in K.SCALES])
        self.scale.SetMinSize(wx.Size(self.FromDIP(180), -1))
        sp.row("Scale", self.scale)
        self.ref_t = sp.text("", 58, tooltip="Each trace at 100 % at its peak here")
        sp.row("Reference time", self.ref_t, "min")
        self.spacing = sp.spin(100, 0, 400)
        sp.row("Spacing", self.spacing, "% of the tallest")
        self.skew = sp.spin(0, 0, 30)
        sp.row("Skew", self.skew, "% of the time span")
        self.reverse = sp.check("First file at the bottom", False)
        self.colours = sp.choice([v for k, v in K.COLOURS])
        sp.row("Colours", self.colours)
        self.lw = sp.text("0.8", 46)
        sp.row("Line width", self.lw, "pt")
        self.fill = sp.check("Fill under the traces", False)
        self.labels = sp.choice([v for k, v in K.LABELS])
        sp.row("Labels", self.labels)
        self.label_text = sp.choice([v for k, v in K.LABEL_TEXT])
        self.label_text.SetMinSize(wx.Size(self.FromDIP(180), -1))
        sp.row("Label text", self.label_text)
        self.rt_labels = sp.choice([v for k, v in K.RT_LABELS])
        self.rt_labels.SetMinSize(wx.Size(self.FromDIP(180), -1))
        sp.row("Retention times", self.rt_labels)
        self.rt_min = sp.text("10", 46, tooltip="Smaller peaks get no label")
        sp.row("Threshold", self.rt_min, "%")
        self.guides = sp.text("", 150, tooltip="Times separated by commas")
        sp.row("Guide lines", self.guides, "min")
        self.yaxis = sp.choice([v for k, v in K.YAXES])
        self.yaxis.SetMinSize(wx.Size(self.FromDIP(150), -1))
        sp.row("Y axis", self.yaxis)

        sp.section("Peak area / kinetics", T.GROUP["purple"] if "purple" in T.GROUP else T.GROUP["blue"])
        self.area_lo, self.area_hi = sp.text("", 58), sp.text("", 58)
        sp.row("Selected region", self.area_lo, "to", self.area_hi, "min")
        self.area_ref = sp.choice(["largest area (automatic)"])
        self.area_ref.SetMinSize(wx.Size(self.FromDIP(180), -1))
        sp.row("100 % reference", self.area_ref)
        self.area_x = sp.text("", 150, tooltip="One per trace in plotted order, e.g. 0, 5, 10")
        sp.row("Times / values", self.area_x)
        self.area_xlabel = sp.text("Reaction time", 150)
        sp.row("X axis", self.area_xlabel)
        self.area_y = sp.choice(["Actual area", "Relative area (%)"])
        sp.row("Plot", self.area_y)
        sp.buttons(U.flat(sp, "Use region", "primary", icon="integrate", handler=self.on_area_fields),
                   U.flat(sp, "Clear", handler=self.clear_area))
        sp.buttons(U.flat(sp, "Area vs X…", icon="integrate", handler=self.plot_area_graph),
                   U.flat(sp, "Excel…", icon="export", handler=self.on_export_areas_xlsx))
        sp.buttons(U.flat(sp, "Kinetic fitting…", icon="integrate", handler=self.open_kinetics))

        sp.section("Spectra (m/z tool)", T.GROUP["blue"])
        self.spec_avg = sp.choice([v for k, v in SPEC_AVG])
        self.spec_avg.SetMinSize(wx.Size(self.FromDIP(180), -1))
        sp.row("Averaging", self.spec_avg)
        self.spec_bg = sp.check("Subtract the spectrum just before the peak", True)
        self.spec_nlab = sp.spin(8, 0, 30)
        sp.row("Ions labelled", self.spec_nlab, "strongest")

        sp.section("Export", T.GROUP["yellow"])
        sp.buttons(U.flat(sp, "Traces (CSV)…", icon="export", handler=self.on_export_traces),
                   U.flat(sp, "Table (CSV)…", icon="export", handler=self.on_export_table))

        self.flist.Bind(wx.EVT_CHECKLISTBOX, self.on_tick)
        self.flist.Bind(wx.EVT_KEY_DOWN, self._list_key)
        self.flist.Bind(wx.EVT_LISTBOX, lambda e: self.choose(self.flist.GetSelection()))
        self.sort.Bind(wx.EVT_CHOICE, self.on_sort)
        self.blank.Bind(wx.EVT_CHOICE, self.on_blank)
        for c in (self.f_label, self.f_shift):
            c.Bind(wx.EVT_TEXT_ENTER, self.on_entry_fields)
            c.Bind(wx.EVT_KILL_FOCUS, lambda e: (self.on_entry_fields(), e.Skip()))
        self.f_colour.Bind(wx.EVT_COLOURPICKER_CHANGED, lambda e: self.set_entry_colour(
            self.f_colour.GetColour().GetAsString(wx.C2S_HTML_SYNTAX)))
        self.wl_pick.Bind(wx.EVT_CHOICE, self.on_wl_pick)
        for c in (self.signal, self.pol, self.base, self.layout_c, self.scale, self.colours, self.labels,
                  self.label_text, self.rt_labels, self.yaxis):
            c.Bind(wx.EVT_CHOICE, self.on_change)
        for c in (self.own_wl, self.mz_ppm, self.reverse, self.fill):
            c.Bind(wx.EVT_CHECKBOX, self.on_change)
        for c in (self.smooth, self.spacing, self.skew):
            c.Bind(wx.EVT_SPINCTRL, self.on_change)
        self.spec_avg.Bind(wx.EVT_CHOICE, self.on_spec_opts)
        self.spec_bg.Bind(wx.EVT_CHECKBOX, self.on_spec_opts)
        self.spec_nlab.Bind(wx.EVT_SPINCTRL, self.on_spec_opts)
        for c in (self.wl, self.bw, self.mz, self.mz_win, self.roll, self.t0, self.t1, self.align, self.align_win,
                   self.ref_t, self.lw, self.rt_min, self.guides):
            c.Bind(wx.EVT_TEXT_ENTER, self.on_change)
            c.Bind(wx.EVT_KILL_FOCUS, lambda e: (self.on_change(), e.Skip()))
        for c in (self.area_lo, self.area_hi):
            c.Bind(wx.EVT_TEXT_ENTER, self.on_area_fields)
        self.area_ref.Bind(wx.EVT_CHOICE, self.on_area_ref)
        self.area_x.Bind(wx.EVT_TEXT_ENTER, lambda e: self.refresh_area())
        self.area_x.Bind(wx.EVT_TEXT, lambda e: (self._bind_area_x(self.drawn), e.Skip()))  # typed: in the order shown
        self.area_xlabel.Bind(wx.EVT_TEXT_ENTER, lambda e: self.refresh_area())
        self.area_y.Bind(wx.EVT_CHOICE, lambda e: self.refresh_area())

    def _tips(self):
        """Tooltips of the controls made without one."""
        for c, tip in (
                (self.wl_pick, "Common wavelengths"),
                (self.mz_ppm, "Window in ppm"),
                (self.smooth, "Savitzky Golay points (0: none)"),
                (self.t0, "Empty: start of the run"),
                (self.t1, "Empty: end of the run"),
                (self.align_win, "Search window (min)"),
                (self.rt_labels, "Times of each run itself, before any shift")):
            try:
                c.SetToolTip(tip)
            except Exception:
                pass

    def _list_key(self, e):
        """Delete in the list of files unticks the chosen file (Space ticks and unticks it)."""
        i = self.flist.GetSelection()
        if e.GetKeyCode() in (wx.WXK_DELETE, wx.WXK_NUMPAD_DELETE) and 0 <= i < len(self.entries):
            if self.entries[i]["include"]:
                self.entries[i]["include"] = False
                self.flist.Check(i, False)
                self.replot(keep_view=True)
        else:
            e.Skip()

    def _fill_controls(self):
        s = self.s
        self._busy = True
        try:
            self.signal.SetSelection(_index(K.SIGNALS, s["signal"]))
            self.wl.SetValue("" if s["wl"] is None else "%g" % s["wl"])
            self.bw.SetValue("%g" % s["bw"])
            self.own_wl.SetValue(bool(s["own_wl"]))
            self.pol.SetSelection(0 if s["polarity"] == "+" else 1)
            self.mz.SetValue("" if s["mz"] is None else "%g" % s["mz"])
            self.mz_win.SetValue("%g" % s["mz_win"])
            self.mz_ppm.SetValue(bool(s["mz_ppm"]))
            self.smooth.SetValue(int(s["smooth"]))
            self.base.SetSelection(_index(K.BASELINES, s["baseline"]))
            self.roll.SetValue("%g" % s["rolling_min"])
            self.align_win.SetValue("%g" % s["align_win"])
            self.layout_c.SetSelection(_index(K.LAYOUTS, s["layout"]))
            self.scale.SetSelection(_index(K.SCALES, s["scale"]))
            self.spacing.SetValue(int(round(K.spacing_of(s))))
            self.skew.SetValue(int(round(K.skew_of(s))))
            self.reverse.SetValue(bool(s["reverse"]))
            self.colours.SetSelection(_index(K.COLOURS, s["colours"]))
            self.lw.SetValue("%g" % s["lw"])
            self.fill.SetValue(bool(s["fill"]))
            self.labels.SetSelection(_index(K.LABELS, s["labels"]))
            self.label_text.SetSelection(_index(K.LABEL_TEXT, s["label_text"]))
            self.rt_labels.SetSelection(_index(K.RT_LABELS, s["rt_labels"]))
            self.rt_min.SetValue("%g" % s["rt_min"])
            self.yaxis.SetSelection(_index(K.YAXES, s["yaxis"]))
            self.spec_avg.SetSelection(_index(SPEC_AVG, s["spec_avg"]))
            self.spec_bg.SetValue(bool(s["spec_bg"]))
            self.spec_nlab.SetValue(int(s["spec_labels"]))
        finally:
            self._busy = False
        self._enable()

    def _enable(self):
        k = self.s["signal"]
        pda = k in ("pda", "max")
        for c in (self.wl, self.wl_pick, self.bw, self.own_wl):
            c.Enable(k == "pda")
        for c in (self.pol,):
            c.Enable(not pda)
        for c in (self.mz, self.mz_win, self.mz_ppm):
            c.Enable(k == "xic")
        self.roll.Enable(self.s["baseline"] == "rolling")
        self.ref_t.Enable(self.s["scale"] == "ref")
        st = K.offset_layout(self.s)
        for c in (self.spacing, self.skew, self.reverse):
            c.Enable(st)
        self.rt_min.Enable(self.s["rt_labels"] == "peaks")

    def _read_controls(self):
        s = self.s
        s["signal"] = _key(K.SIGNALS, self.signal.GetSelection())
        s["wl"] = _fnum(self.wl.GetValue(), s["wl"] if self.wl.GetValue().strip() else None)
        s["bw"] = max(_fnum(self.bw.GetValue(), 4.0), 0.0)
        s["own_wl"] = self.own_wl.GetValue()
        s["polarity"] = "+" if self.pol.GetSelection() == 0 else "-"
        s["mz"] = _fnum(self.mz.GetValue(), None)
        s["mz_win"] = abs(_fnum(self.mz_win.GetValue(), 0.5) or 0.5)
        s["mz_ppm"] = self.mz_ppm.GetValue()
        s["smooth"] = int(self.smooth.GetValue())
        s["baseline"] = _key(K.BASELINES, self.base.GetSelection())
        s["rolling_min"] = max(_fnum(self.roll.GetValue(), 1.0), 0.01)
        a, b = _fnum(self.t0.GetValue(), None), _fnum(self.t1.GetValue(), None)
        if a is not None and b is not None and b < a:
            a, b = b, a
        s["t0"], s["t1"] = a, b
        s["align"] = _fnum(self.align.GetValue(), None)
        s["align_win"] = max(_fnum(self.align_win.GetValue(), 0.3), 0.01)
        lay_old = s["layout"]
        s["layout"] = _key(K.LAYOUTS, self.layout_c.GetSelection())
        s["scale"] = _key(K.SCALES, self.scale.GetSelection())
        s["ref_t"] = _fnum(self.ref_t.GetValue(), None)
        # spacing and skew of the layout shown in the controls (Offset keeps its own pair)
        sk, kk = ("off_spacing", "off_skew") if lay_old == "offset" else ("spacing", "skew")
        s[sk] = int(self.spacing.GetValue())
        s[kk] = int(self.skew.GetValue())
        s["reverse"] = self.reverse.GetValue()
        s["colours"] = _key(K.COLOURS, self.colours.GetSelection())
        s["lw"] = min(max(_fnum(self.lw.GetValue(), 0.8), 0.2), 4.0)
        s["fill"] = self.fill.GetValue()
        s["labels"] = _key(K.LABELS, self.labels.GetSelection())
        s["label_text"] = _key(K.LABEL_TEXT, self.label_text.GetSelection())
        s["rt_labels"] = _key(K.RT_LABELS, self.rt_labels.GetSelection())
        s["rt_min"] = min(max(_fnum(self.rt_min.GetValue(), 10.0), 0.0), 100.0)
        s["guides"] = sorted(set(round(v, 4) for v in U._numlist(self.guides.GetValue())))
        s["yaxis"] = _key(K.YAXES, self.yaxis.GetSelection())

    def _save(self):
        keep = {k: v for k, v in self.s.items() if k not in ("t0", "t1", "align", "ref_t", "guides")}
        T._save({SETTINGS_KEY: keep})

    def method_get(self):
        """Settings of the comparison in a method preset (method_presets.py; not the files, their names,
        colours and shifts, nor the blank)."""
        return {k: list(x) if isinstance(x, list) else x for k, x in self.s.items()}

    def method_set(self, v):
        """The values given (the others stay) as if chosen in the panel; drawn again."""
        old = dict(self.s)
        self.s.update({k: list(x) if isinstance(x, list) else x for k, x in v.items() if k in self.s})
        if v.get("wl") is not None:
            self._wl_from_view = True  # (not replaced by the wavelength of the PDA view)
        self._fill_controls()
        for c, k in ((self.t0, "t0"), (self.t1, "t1"), (self.align, "align"), (self.ref_t, "ref_t")):
            c.ChangeValue("" if self.s[k] is None else "%.10g" % self.s[k])
        self.guides.ChangeValue(", ".join("%g" % g for g in self.s["guides"]))
        self._save()
        self.card.ylock = None
        if self.s["label_text"] != old["label_text"]:
            self._fill_list()
        self.replot(keep_view=(self.s["t0"], self.s["t1"]) == (old["t0"], old["t1"]))
        if self.mzpick is not None:
            if (self.s["spec_avg"], self.s["spec_bg"]) != (old["spec_avg"], old["spec_bg"]):
                self._mz_compute()
            for p in ("+", "-"):
                self.plot_mz(p, keep_view=True)

    def on_change(self, e=None):
        if self._busy:
            return
        old = dict(self.s)
        self._read_controls()
        if self.s == old:
            return
        if self.s["layout"] != old["layout"]:
            if self.s["layout"] == "offset" and self.s["labels"] == "right":
                self.s["labels"] = "outside"  # the names right of the frame, as in a 2D waterfall figure
            elif old["layout"] == "offset" and self.s["labels"] == "outside":
                self.s["labels"] = "right"
            self._fill_controls()  # the spacing and skew of the new layout
        self._enable()
        self._save()
        self.card.ylock = None  # a y range fixed with the zoom box would hide traces moved by the change
        if self.s["label_text"] != old["label_text"]:
            self._fill_list()  # the names in the list and in the blank choice
        self.replot(keep_view=(self.s["t0"], self.s["t1"]) == (old["t0"], old["t1"]))

    def on_wl_pick(self, e=None):
        i = self.wl_pick.GetSelection()
        if i > 0:
            self.wl.SetValue("%d" % K.COMMON_WL[i - 1])
            self.wl_pick.SetSelection(0)
            if self.s["signal"] != "pda":
                self.signal.SetSelection(_index(K.SIGNALS, "pda"))
            self.on_change()

    # ------------------------------------------------------------------ files
    def _docs(self):
        fr = self.frame
        return [d for d in getattr(fr, "docs", []) if d.attrs.get("path") and not d.attrs.get("loading") and
                (d.attrs.get("ms_data") is not None or d.attrs.get("pda_data") is not None)]

    def files_changed(self):
        """The open files changed (opened, closed): the list follows, keeping the choices made."""
        if not self:  # the window was closed meanwhile (called after a file was read)
            return
        docs = self._docs()
        old = {id(en["doc"]): en for en in self.entries}
        ents = [en for en in self.entries if en["doc"] in docs]
        for d in docs:
            if id(d) not in old:
                ents.append({"doc": d, "include": True, "label": None, "colour": None, "shift": 0.0})
        for k in list(self._cache):
            if k[0] not in {id(d) for d in docs}:
                del self._cache[k]
        for k in list(self.props["traces"]):  # the line of a closed file (a new file could get its id)
            if k not in {id(d) for d in docs}:
                del self.props["traces"][k]
        self.entries = ents
        if self.sel is not None and self.sel >= len(self.entries):
            self.sel = None
        self._pick_alive()
        if docs and (self.s["wl"] is None or not self._wl_from_view):
            # the first wavelength: that of the PDA view of the file shown (else of the first file with
            # one), as the user saw it before switching to Compare; the saved one when no file has one
            act = getattr(self.frame, "active", None)
            for d in sorted(docs, key=lambda d: d is not act):
                p = getattr(d.attrs.get("pda"), "wl_sel", None)
                if p:
                    self.s["wl"] = float(round(p))
                    self.wl.SetValue("%g" % self.s["wl"])
                    self._wl_from_view = True
                    break
        self._fill_list()
        self.replot(keep_view=False)

    def request_files_changed(self):
        """One Compare redraw after a group of dropped files finishes loading."""
        if getattr(getattr(self.frame, "_project", None), "restoring", False):
            return  # project recovery installs the saved Compare state itself
        if getattr(self, "_files_update_pending", False):
            return
        self._files_update_pending = True
        wx.CallLater(60, self._flush_files_changed)

    def _flush_files_changed(self):
        if not self:
            return
        fr = self.frame
        if getattr(getattr(fr, "_project", None), "restoring", False):
            self._files_update_pending = False
            return
        if getattr(fr, "_load_queue", None) or any(d.loading for d in fr.docs):
            wx.CallLater(60, self._flush_files_changed)
            return
        self._files_update_pending = False
        self.files_changed()

    def _prune(self):
        """Files closed meanwhile are dropped from the list (a file closed while the Compare view is
        shown: the window chooses another file before the list follows). True when one was."""
        docs = getattr(self.frame, "docs", None)
        if docs is None:
            return False
        keep = [en for en in self.entries if en["doc"] in docs and en["doc"].attrs.get("path")]
        if len(keep) == len(self.entries):
            return False
        cur = self.entries[self.sel] if self.sel is not None and self.sel < len(self.entries) else None
        self.entries = keep
        self.sel = keep.index(cur) if cur in keep else None
        return True

    def info(self, d):
        """Sample information of a file (cached with the file)."""
        si = d.attrs.get("sample_info")
        if si is None:
            si = {}
            try:
                import lcms_data
                si = lcms_data.read_sample_info(d.attrs.get("path") or "") or {}
            except Exception:
                si = {}
            d.attrs["sample_info"] = si
        return si

    def _base_label(self, en, k=None):
        d = en["doc"]
        fn = os.path.splitext(os.path.basename((d.attrs.get("path") or "").rstrip("\\/")))[0]
        sn = (self.info(d).get("sample_name") or "").strip()
        k = k or self.s["label_text"]
        if k == "file" or not sn:
            return fn
        if k == "both":
            return "%s (%s)" % (sn, fn)
        return sn

    def label_of(self, en):
        if en.get("label"):
            return en["label"]
        lab = self._base_label(en)
        if self.s["label_text"] == "sample" and sum(1 for e in self.entries if not e.get("label") and
                                                    self._base_label(e) == lab) > 1:
            return self._base_label(en, "both")  # runs of the same sample: told apart by the file name
        return lab

    def _fill_list(self):
        if self.blank_doc is not None and self.blank_doc not in [en["doc"] for en in self.entries]:
            self.blank_doc = None
        self.blank.Set(["none"] + [self.label_of(en) for en in self.entries])
        self.blank.SetSelection(next((i + 1 for i, en in enumerate(self.entries) if en["doc"] is self.blank_doc), 0))
        self.flist.Freeze()
        try:
            self.flist.Set([self.label_of(en) + ("  (blank)" if en["doc"] is self.blank_doc else "")
                            for en in self.entries])
            for i, en in enumerate(self.entries):
                self.flist.Check(i, bool(en["include"]))
            if self.sel is not None:
                self.flist.SetSelection(self.sel)
        finally:
            self.flist.Thaw()
        self._fill_entry_fields()

    def _fill_entry_fields(self):
        en = self.entries[self.sel] if self.sel is not None and self.sel < len(self.entries) else None
        for c in (self.f_label, self.f_colour, self.f_shift, self.f_colour_reset):
            c.Enable(en is not None)
        self._busy = True
        try:
            self.f_label.SetValue(en.get("label") or "" if en else "")
            self.f_shift.SetValue("%g" % en.get("shift", 0.0) if en else "0")
            col = (en.get("colour") or self._colour_of(self.sel)) if en else U.PAL[0]
            self.f_colour.SetColour(wx.Colour(col))
        finally:
            self._busy = False

    def choose(self, i):
        self.on_entry_fields()  # a label or shift typed for the file chosen so far is kept (focus order)
        self.sel = i if 0 <= i < len(self.entries) else None
        if self.sel is not None and self.flist.GetSelection() != self.sel:
            self.flist.SetSelection(self.sel)
        self._fill_entry_fields()
        self.replot(keep_view=True, quick=True)

    def select_doc(self, doc):
        """The file chosen in the list of open files is chosen here too."""
        for i, en in enumerate(self.entries):
            if en["doc"] is doc:
                if self.sel != i:
                    self.choose(i)
                return

    def on_blank(self, e=None):
        i = self.blank.GetSelection() - 1
        self.blank_doc = self.entries[i]["doc"] if 0 <= i < len(self.entries) else None
        self._fill_list()
        self.replot(keep_view=True)

    def on_tick(self, e):
        i = e.GetInt()
        if 0 <= i < len(self.entries):
            self.entries[i]["include"] = self.flist.IsChecked(i)
            self.replot(keep_view=True)

    def move(self, step):
        i = self.sel
        if i is None or not (0 <= i + step < len(self.entries)):
            return
        self.entries[i], self.entries[i + step] = self.entries[i + step], self.entries[i]
        self.sel = i + step
        self._fill_list()
        self.replot(keep_view=True)

    def reverse_order(self):
        self.entries.reverse()
        if self.sel is not None:
            self.sel = len(self.entries) - 1 - self.sel
        self._fill_list()
        self.replot(keep_view=True)

    def on_sort(self, e=None):
        k = self.sort.GetSelection()
        self.sort.SetSelection(0)
        if k <= 0:
            return
        cur = self.entries[self.sel] if self.sel is not None else None
        if k == 1:
            self.entries.sort(key=lambda en: os.path.basename(en["doc"].attrs.get("path") or "").lower())
        elif k == 2:
            def when(en):
                a = self.info(en["doc"]).get("acquired")
                return a if isinstance(a, datetime.datetime) else datetime.datetime.max
            self.entries.sort(key=when)
        elif k == 3:
            self.entries.sort(key=lambda en: (self.info(en["doc"]).get("sample_name") or "").lower())
        else:
            order = {id(d): n for n, d in enumerate(getattr(self.frame, "docs", []))}
            self.entries.sort(key=lambda en: order.get(id(en["doc"]), 0))
        self.sel = self.entries.index(cur) if cur is not None else None
        self._fill_list()
        self.replot(keep_view=True)

    def on_entry_fields(self, e=None):
        if self._busy or self.sel is None or self.sel >= len(self.entries):
            return
        en = self.entries[self.sel]
        lab = self.f_label.GetValue().strip()
        sh = _fnum(self.f_shift.GetValue(), 0.0) or 0.0
        if (lab or None) == en.get("label") and abs(sh - en.get("shift", 0.0)) < 1e-12:
            return
        en["label"], en["shift"] = lab or None, sh
        self._fill_list()
        self.replot(keep_view=True)

    def set_entry_colour(self, col):
        if self._busy or self.sel is None or self.sel >= len(self.entries):
            return
        self.entries[self.sel]["colour"] = col
        self._fill_entry_fields()
        self.replot(keep_view=True, quick=True)

    def on_row(self, i):
        """A row of the table: its file is chosen."""
        if 0 <= i < len(self.drawn):
            self.choose(self.entries.index(self.drawn[i]["entry"]))

    # ------------------------------------------------------------------ traces
    def _sig_key(self):
        s = self.s
        return (s["signal"], s["wl"], s["bw"], s["own_wl"], s["polarity"], s["mz"], s["mz_win"], s["mz_ppm"])

    def raw(self, en):
        """(t, y, units, description) of a file, or a TraceError (kept per file and signal)."""
        d = en["doc"]
        own = None
        try:
            own = d.attrs.get("pda").wl_sel if d.attrs.get("pda") is not None else None
        except Exception:
            own = None
        key = (id(d),) + self._sig_key() + ((own,) if self.s["own_wl"] else ())
        if key not in self._cache:
            try:
                self._cache[key] = K.raw_trace(d.attrs.get("pda_data"), d.attrs.get("ms_data"), self.s, own)
            except K.TraceError as ex:
                self._cache[key] = ex
            except Exception as ex:  # a reader error: reported, the other files are drawn
                U._log("compare: %s: %s" % (d.attrs.get("file_name"), ex))
                self._cache[key] = K.TraceError("could not be read (%s)" % ex)
            if len(self._cache) > 200:
                self._cache.pop(next(iter(self._cache)))
        return self._cache[key]

    def _colour_of(self, i):
        en = self.entries[i] if i is not None and i < len(self.entries) else None
        if en is not None and en.get("colour"):
            return en["colour"]
        k = self.s["colours"]
        inc = [j for j, e in enumerate(self.entries) if e["include"] and e["doc"] is not self.blank_doc]
        n = max(len(inc), 1)
        pos = inc.index(i) if i in inc else 0
        if k == "gradient":
            return _mix(GRADIENT[0], GRADIENT[1], pos / max(n - 1, 1))
        if k == "black":
            return "#1B2330"
        if k == "dark":
            return U.TRACE
        return U.PAL[pos % len(U.PAL)]

    def compute(self):
        """The traces to draw: processed, stacked; and the notes (files left out, differences)."""
        s = self.s
        out, notes, raws = [], [], []
        blank = None
        if self.blank_doc is not None:
            ben = next((en for en in self.entries if en["doc"] is self.blank_doc), None)
            r = self.raw(ben) if ben is not None else None
            if isinstance(r, K.TraceError):
                notes.append("The blank %s has %s: nothing subtracted." % (self.label_of(ben), r))
            elif r is not None:
                sh = float(ben.get("shift", 0.0) or 0.0)
                blank = (r[0] + sh, r[1])
        for i, en in enumerate(self.entries):
            if not en["include"] or en["doc"] is self.blank_doc:
                continue
            r = self.raw(en)
            if isinstance(r, K.TraceError):
                notes.append("\u2192 %s: %s" % (self.label_of(en), r))  # (marked: a file left out)
                continue
            if blank is not None:
                sh = float(en.get("shift", 0.0) or 0.0)
                r = (r[0], K.subtract_blank(r[0] + sh, r[1], blank[0], blank[1]), r[2], r[3])
            raws.append((i, en, r))
        shifts = [0.0] * len(raws)
        if s["align"] is not None and raws:
            apex = K.alignment_apexes([(r[0], r[1]) for _, _, r in raws], s["align"], s["align_win"], s["smooth"])
            shifts = [0.0 if a is None else s["align"] - a for a in apex]
            for (i, en, r), a in zip(raws, apex):
                if a is None:  # drawn as it is: say so (it looked aligned)
                    notes.append("%s: no peak within %.2f \u00b1 %g min, not aligned" % (self.label_of(en), s["align"],
                                                                                    s["align_win"]))
        for (i, en, r), sh in zip(raws, shifts):
            t, y, units, desc = r
            pr = K.process(t, y, s, sh + float(en.get("shift", 0.0) or 0.0))
            if not pr["ok"]:
                notes.append("\u2192 %s: no data in the time window" % self.label_of(en))
                continue
            if pr.get("ref_missing"):  # it cannot be put in % of a peak it does not have
                notes.append("\u2192 %s: no signal at the reference time %.2f min" % (self.label_of(en), s["ref_t"]))
                continue
            out.append({"i": i, "entry": en, "proc": pr, "units": units, "desc": desc, "align": sh,
                        "label": self.label_of(en), "colour": self._colour_of(i),
                        "scaled": bool(pr.get("scaled"))})  # in % only when a factor was applied
        offs = K.stack([o["proc"] for o in out], K.spacing_of(s), s["reverse"], K.skew_of(s), s["layout"])
        for o, (dx, dy) in zip(out, offs):
            o["dx"], o["dy"] = dx, dy
        infos = []
        for o in out:
            d = o["entry"]["doc"]
            pda = d.attrs.get("pda_data")
            ms = d.attrs.get("ms_data")
            t_end = None
            try:
                t_end = float(pda.times[-1]) if pda is not None else float(np.max(ms.rt))
            except Exception:
                t_end = None
            infos.append({"label": o["label"], "t_end": t_end, "method": self.info(d).get("method_file")})
        notes += K.consistency_notes(infos)
        if s["scale"] == "ref" and s["ref_t"] is None:
            notes.append("Give the reference time for the scale.")
        return out, notes

    # ------------------------------------------------------------------ plot
    def _use_legend(self, n):
        """Names in a legend: chosen, or when names at the ends of the traces would sit on top of
        each other (overlay, or stacked without spacing)."""
        s = self.s
        if s["labels"] == "legend":
            return True
        return s["labels"] in ("right", "left", "outside") and n > 1 and (not K.offset_layout(s) or
                                                                          K.spacing_of(s) == 0)

    def replot(self, keep_view=True, quick=False):
        card, ax = self.card, self.card.ax
        if self._prune():
            self._fill_list()
        self._pick_alive()
        view = ax.get_xlim() if (keep_view and card.full) else None
        card.reset()
        self._label_arts = []
        s = self.s
        P = self.props
        units_all = None
        out, notes = self.compute() if self.entries else ([], [])
        left = [n.lstrip("\u2192 ") for n in notes if n.startswith("\u2192")]  # files left out
        notes = [n.lstrip("\u2192 ") for n in notes]
        self.drawn = out
        self.area_results = self._compute_area_results(out)
        if out:
            # % only when every trace was put in % (no reference time, or a trace with no top above 0: not)
            units = set("%" if o["scaled"] else o["units"] for o in out)
            units_all = units.pop() if len(units) == 1 else ""
        ylab = ("Relative intensity (%)" if units_all == "%" else "Absorbance (mAU)" if units_all == "mAU"
                else "Intensity (counts)" if units_all == "counts" else "Intensity")
        self._auto_titles = ("Time (min)", ylab)
        ymode = K.yaxis_mode(s, len(out)) if out else "axis"
        U.style_axes(ax, "Time (min)", ylab if ymode == "axis" else "")
        card.margins = [0.84, 0.52, 0.2, 0.14]
        PP.apply_axes(ax, P, U.INK, y_values=(ymode == "axis"))
        if ymode != "axis" and P.get("ytitle"):
            ax.yaxis.label.set_text("")  # no axis: no title on it
        if not out:
            msg = ("Open two or more files" if not self.entries else
                   "Tick the files to compare" if not any(e["include"] for e in self.entries)
                   else "No trace to draw: " + "; ".join(notes[:3]))
            ax.text(0.5, 0.5, msg, transform=ax.transAxes, ha="center", va="center", fontsize=9, color=U.MUTED,
                    wrap=True)
            card.full = None
            card.set_title("Comparison", "")
            card.set_legend([])
            card.set_badge(None)
            self.table.set_rows([], "; ".join(notes))
            self._table = None
            self._peak_table = None
            self._area_table = None
            card.draw()
            return
        lw = float(s["lw"])
        lo, hi = None, None
        from matplotlib.transforms import blended_transform_factory
        tr_lab = blended_transform_factory(ax.transAxes, ax.transData)
        offset = s["layout"] == "offset" and len(out) > 1
        top_dy = max([o["dy"] for o in out] + [0.0])
        fill_a = min(max(float(P.get("fill_alpha", 16)), 0.0), 100.0) / 100.0
        for n, o in enumerate(out):
            pr, col = o["proc"], o["colour"]
            x, y = pr["t"] + o["dx"], pr["y"] + o["dy"]
            fill = (col, fill_a) if s["fill"] else None
            tp = P["traces"].get(self._trace_key(o), {})
            if offset:  # the front trace (lowest) over the ones behind it, as in a waterfall
                z = 3 + 0.4 * (1.0 - (o["dy"] / top_dy if top_dy > 0 else 0.0)) + (len(out) - n) * 1e-4
            else:
                z = 3 + (len(out) - n) * 0.01
            ln = U.plot_trace(ax, x, y, fill=fill, fill_base=o["dy"], color=col, lw=float(tp.get("lw") or lw),
                              ls=tp.get("ls") or "-", zorder=z)
            ln._cmp = o
            if n < len(self.area_results):
                self._draw_area(ax, o, self.area_results[n], z)
            if o["i"] == self.sel and len(out) > 1:  # the file chosen in the list: a halo, on screen only
                hl = U.plot_trace(ax, x, y, color=col, lw=lw + 3.5, alpha=0.2, zorder=2.5, solid_capstyle="round")
                hl._ui_only = True
                hl._no_scale = True
            lo = float(x[0]) if lo is None else min(lo, float(x[0]))
            # offset: the time axis is that of the front trace (the traces behind it, moved right, are cut
            # at its end as in a waterfall figure); otherwise every trace is shown to its end
            end = float(pr["t"][-1]) if offset else float(x[-1])
            hi = end if hi is None else max(hi, end)
        if self.area_range is not None:
            a, b = self.area_range
            band = ax.axvspan(a, b, color=C["accent"], alpha=0.07, zorder=1)
            band._no_scale = band._area_region = True
            for t in (a, b):
                edge = ax.axvline(t, color=C["accent"], lw=1.0, ls="--", zorder=5)
                edge._no_scale = edge._area_region = True
        self._draw_mz_mark(ax, out)
        if P.get("x0") is not None or P.get("x1") is not None:  # graph properties: the range set there
            w = hi - lo if hi > lo else 1.0
            lo = P["x0"] if P.get("x0") is not None else lo
            hi = P["x1"] if P.get("x1") is not None else hi
            if hi <= lo:  # one end set beyond the other end of the runs: as wide as the runs from it
                if P.get("x0") is not None:
                    hi = lo + w
                else:
                    lo = hi - w
        card.full = (lo, hi) if (lo is not None and hi is not None and hi > lo) else None
        if view or card.full:
            ax.set_xlim(*(view if view else card.full))
        x0, x1 = ax.get_xlim()
        for g in s["guides"]:  # guide lines (their times above the frame: _place_labels)
            if x0 <= g <= x1:
                gl = ax.axvline(g, color="#7B8594", lw=0.8, ls=(0, (4, 3)), zorder=2)
                gl._no_scale = True
        if s["align"] is not None and x0 <= s["align"] <= x1:
            al = ax.axvline(s["align"], color=C["accent"], lw=0.7, ls=(0, (1, 2)), zorder=2, alpha=0.8)
            al._no_scale = True
        card.label_room_pt = (P["size_rt"] + 15.0) if s["rt_labels"] != "none" else 12.0
        card.y_fixed = (P.get("y0"), P.get("y1"))  # graph properties: the ends of the y range set there
        card.autoscale_y(pad=0.06)
        card.set_legend([])
        if self._use_legend(len(out)):  # in the plot itself, so that copied and saved images have it
            from matplotlib.lines import Line2D
            hs = []
            for o in out:
                tp = P["traces"].get(self._trace_key(o), {})
                hs.append(Line2D([], [], color=o["colour"], lw=1.6, ls=tp.get("ls") or "-"))
            lkw = PP.font_kw(P, "size_names", "bold_names")
            prop = {"size": lkw["fontsize"], "weight": lkw.get("fontweight", "normal")}
            if lkw.get("fontfamily"):
                prop["family"] = lkw["fontfamily"]
            lg = ax.legend(hs, [_tex(o["label"]) for o in out], loc=P.get("legend_loc", "best"), prop=prop,
                           frameon=bool(P.get("legend_frame")), handlelength=1.5, labelspacing=0.25,
                           borderaxespad=0.4)
            if P.get("legend_frame"):
                lg.get_frame().set_linewidth(0.6)
                lg.get_frame().set_edgecolor("#7B8594")
            lg.set_zorder(8)
        if ymode != "axis":  # no y axis: no numbers on it
            ax.set_yticks([])
            ax.yaxis.set_minor_locator(__import__("matplotlib").ticker.NullLocator())
        # stacked: the y axis means nothing for the upper traces; a scale bar says how tall a signal is
        if ymode == "bar":
            unit = units_all or ""
            bar = K.nice(0.5 * max([o["proc"]["top"] for o in out] + [0.0]))
            if bar > 0:  # left of the plot, where the y axis would be (it covered the start of the traces)
                xb = -0.012
                b, = ax.plot([xb, xb], [0.0, bar], transform=tr_lab, color="#1B2330", lw=1.6, solid_capstyle="butt",
                             zorder=7, clip_on=False)
                b._no_scale = True
                # from the foot of the bar upwards: centred on a short bar it ran into the time axis
                bkw = PP.font_kw(P, "size_bar")
                ax.text(xb - 0.007, 0.0, ("%s %s" % (K.bar_text(bar), unit)).strip(), transform=tr_lab, ha="right",
                        va="bottom", rotation=90, color="#1B2330", zorder=7, clip_on=False, **bkw)
        if s["signal"] == "pda" and s["own_wl"]:
            sub = ["PDA at the wavelength of each file (bandwidth %g nm)" % s["bw"]]
        else:
            sub = [K.signal_text(s, s["wl"])]
        if s["smooth"] >= 3:
            sub.append("smoothed (%d points)" % s["smooth"])
        if s["baseline"] != "none":
            sub.append("baseline: %s" % dict(K.BASELINES)[s["baseline"]].split(" (")[0].lower())
        if self.blank_doc is not None:
            sub.append("blank subtracted")
        if s["align"] is not None:
            sub.append("aligned at %.2f min" % s["align"])
        if s["scale"] != "abs":
            sub.append(dict(K.SCALES)[s["scale"]].lower().replace(" (100 %)", ""))
        sub.append("%d of %d files" % (len(out), len(self.entries)))
        card.set_title("Comparison", ", ".join(sub))
        card.set_badge(("%d left out" % len(left)) if left else None, tip="\n".join(left))
        if self.area_range is not None or not quick or not self.table.rows:
            self._fill_table(out, notes)
        else:
            self.table.info = "; ".join(notes) if notes else self._table_note()
            self.table.Refresh()
        self._labels_for_frame()
        card.draw()
        if self.mz_on:
            self._mz_retitle()

    def _place_labels(self):
        """Names of the traces, retention times and the times of the guide lines, placed for the
        present size of the plot so that they do not run into each other. Also called for copied
        and saved images (PlotCard.relabel_fn), whose size and text differ from the tile."""
        ax = self.card.ax
        for a in getattr(self, "_label_arts", []):
            try:
                a.remove()
            except Exception:
                pass
        self._label_arts = []
        out = self.drawn
        if not out:
            return
        try:
            fig = ax.figure
            r = fig.canvas.get_renderer()
            bb = ax.get_window_extent(r)
        except Exception as ex:
            U._log("compare labels: %s" % ex)
            return
        from matplotlib import patheffects
        from matplotlib.transforms import blended_transform_factory
        s = self.s
        P = self.props
        nkw = PP.font_kw(P, "size_names", "bold_names")
        rkw = PP.font_kw(P, "size_rt", "bold_rt")
        gkw = PP.font_kw(P, "size_guides")
        pt = fig.dpi / 72.0
        x0, x1 = ax.get_xlim()
        halo = [patheffects.withStroke(linewidth=2.2, foreground="white")]
        taken = []  # display boxes (x0, y0, x1, y1) of the labels placed

        def keep(a):
            a._dec_label = True  # the text size of copied images applies (PlotCard._render)
            self._label_arts.append(a)
            return a

        def free(bx):
            return not any(t[0] < bx[2] and bx[0] < t[2] and t[1] < bx[3] and bx[1] < t[3] for t in taken)

        lg = ax.get_legend()
        if lg is not None:
            try:
                e = lg.get_window_extent(r)
                taken.append((e.x0, e.y0, e.x1, e.y1))
            except Exception:
                pass
        # times of the guide lines above the frame: each centred on its line or nudged a little (at
        # most a third of its width) beside a close one; lines closer still get a row above (four
        # rows at most: the time of a line that finds no room is left out, the line stays)
        gl = [g for g in sorted(s["guides"]) if x0 <= g <= x1]
        if gl:
            tr_top = blended_transform_factory(ax.transData, ax.transAxes)
            rows = [[] for _ in range(4)]  # per row: [centre, width, wanted centre, artist]
            gap = 3 * pt
            h = 0.0
            for g in gl:
                a = keep(ax.annotate("%.2f" % g, (g, 1.0), xycoords=tr_top, xytext=(0, 2), textcoords="offset points",
                                     ha="center", va="bottom", color="#5B6475", annotation_clip=False,
                                     zorder=6, **gkw))
                e = a.get_window_extent(r)
                w, h = e.width, max(h, e.height)
                cx = float(ax.transData.transform((g, 0.0))[0])
                for k, row in enumerate(rows):
                    if not row:
                        row.append([cx, w, cx, a])
                        break
                    last = row[-1]
                    d = last[0] + last[1] / 2 + gap - (cx - w / 2)
                    if d <= 0:
                        row.append([cx, w, cx, a])
                        break
                    prev = row[-2] if len(row) > 1 else None
                    room = (last[0] - last[1] / 2) - (prev[0] + prev[1] / 2 + gap) if prev else 1e9
                    if d / 2 <= w / 3 and abs(last[0] - d / 2 - last[2]) <= last[1] / 3 and room >= d / 2:
                        last[0] -= d / 2
                        row.append([cx + d / 2, w, cx, a])
                        break
                    if k == len(rows) - 1:  # no room in any row
                        a.remove()
                        self._label_arts.remove(a)
            for k, row in enumerate(rows):
                for c, w, cx, a in row:
                    dy = 2.0 + k * (h / pt + 1.0)
                    a.xyann = ((c - cx) / pt, dy)
                    y0 = bb.y1 + dy * pt
                    taken.append((c - w / 2, y0, c + w / 2, y0 + h))
            n_rows = sum(1 for row in rows if row)
            top = 0.04 + n_rows * (h / fig.dpi + 0.02)  # room above the frame for the rows (inches)
        else:
            top = 0.0
        titled = P.get("title") or P.get("letter")
        if titled:  # the title or panel letter (graph properties) above those times
            pad = 6.0 + (max(0.0, top - 0.04) * 72.0 if gl else 0.0)
            PP.apply_titles(ax, P, U.INK, pad=pad)
            for t in (ax.title, ax._left_title):
                t._dec_label = True  # the text size of copied images applies (PlotCard._render)
        if s["labels"] == "outside" and not self._use_legend(len(out)):
            self._labels_outside(ax, r, bb, x0, x1, pt, keep, taken, nkw)
        if titled:  # a long title on more lines, clear of the letter (its room: without the names beside)
            cw, ch = self.card.canvas.GetClientSize()
            image = abs(fig.bbox.width - cw) > 1.5 or abs(fig.bbox.height - ch) > 1.5  # (or the report)
            PP.fit_title(ax, r, fitted=image)
            top = max(top, PP.title_room(P, pad, ax, r))
        self.card.margins[3] = max(0.14, top)
        # names at the ends of the traces: in the order of the stack, each within the band of its
        # trace (above the trace under the name, below the next baseline), apart from each other
        if s["labels"] in ("right", "left") and not self._use_legend(len(out)):
            right = s["labels"] == "right"
            tr_lab = blended_transform_factory(ax.transAxes, ax.transData)
            levels = sorted(set(round(o["dy"], 9) for o in out))
            items = []
            for n, o in enumerate(out):
                pr = o["proc"]
                xs = pr["t"] + o["dx"]
                m = (xs >= x0) & (xs <= x1)
                if not np.any(m):
                    continue
                a = keep(ax.text(0.995 if right else 0.005, 0.0, _tex(o["label"]), transform=tr_lab,
                                 ha="right" if right else "left", va="bottom", color=_ink(o["colour"]),
                                 zorder=6, clip_on=False, path_effects=halo, **nkw))
                e = a.get_window_extent(r)
                # the time range under the name (its end at 0.5 % of the width from the frame, 3 pt more)
                wd = (e.width + 3 * pt) / max(bb.width, 1.0) * (x1 - x0) + 0.005 * (x1 - x0)
                mm = m & ((xs >= x1 - wd) if right else (xs <= x0 + wd))
                if not np.any(mm):
                    mm = m
                ytop = o["dy"] + float(np.nanmax(pr["y"][mm]))  # the trace under the name
                yd = float(ax.transData.transform((x0, ytop))[1]) + 2 * pt
                base = float(ax.transData.transform((x0, o["dy"]))[1])
                k = levels.index(round(o["dy"], 9))
                if k + 1 < len(levels):  # below the next baseline (it went into the band of the next trace)
                    top = float(ax.transData.transform((x0, levels[k + 1]))[1])
                    if yd + e.height > top - 0.5 * pt:
                        # a peak under the name rises into the next band (a long name, a narrow plot): over
                        # it would not help, so the name sits at the level of the end of the trace
                        me = m & ((xs >= x1 - 0.03 * (x1 - x0)) if right else (xs <= x0 + 0.03 * (x1 - x0)))
                        if np.any(me):
                            yd = float(ax.transData.transform((x0, o["dy"] + float(np.nanmax(pr["y"][me]))))[1]) + 2 * pt
                    yd = min(yd, top - e.height - 0.5 * pt)
                yd = max(yd, base + 0.5 * pt)
                items.append([yd, e.height, a, e, o["dy"], n])
            items.sort(key=lambda it: (it[4], -it[5]))  # bottom of the stack first
            gap = 0.5 * pt
            for k in range(1, len(items)):
                items[k][0] = max(items[k][0], items[k - 1][0] + items[k - 1][1] + gap)
            if items and items[-1][0] + items[-1][1] > bb.y1 - pt:
                items[-1][0] = bb.y1 - pt - items[-1][1]
                for k in range(len(items) - 2, -1, -1):
                    items[k][0] = min(items[k][0], items[k + 1][0] - items[k][1] - gap)
            inv = ax.transData.inverted()
            for yd, hh, a, e, _dy, _n in items:
                a.set_y(float(inv.transform((bb.x0, yd))[1]))
                xa = bb.x1 - 0.005 * bb.width if right else bb.x0 + 0.005 * bb.width
                taken.append((xa - e.width, yd, xa, yd + hh) if right else (xa, yd, xa + e.width, yd + hh))
        # retention times: the time of the run itself (before alignment and shift), above the peak;
        # when another label is there (overlay, close spacing): beside the peak top or higher up;
        # left out when there is no room in the frame
        if s["rt_labels"] != "none":
            import ms_deconv as D
            cand = []
            for o in out:
                pr = o["proc"]
                m = (pr["t"] + o["dx"] >= x0) & (pr["t"] + o["dx"] <= x1)
                if np.sum(m) < 3:
                    continue
                tt, yy = pr["t"][m], pr["y"][m]
                if s["rt_labels"] == "main":
                    j = int(np.argmax(yy))
                    # a peak top, not the edge of the view or a flat trace (a blank minus itself)
                    idx = [j] if 0 < j < len(yy) - 1 and yy[j] > float(np.min(yy)) else []
                else:
                    idx = D.label_maxima_index(tt, yy, 8, (x1 - x0) / 30.0, s["rt_min"] / 100.0)
                shift = o["align"] + float(o["entry"].get("shift", 0.0) or 0.0)
                for j in idx:
                    xy = (float(tt[j] + o["dx"]), float(yy[j] + o["dy"]))
                    cand.append((float(ax.transData.transform(xy)[1]), xy, "%.2f" % (tt[j] - shift), o))
            cand.sort(key=lambda c: c[0])
            for yd, xy, txt, o in cand:
                a = keep(ax.annotate(txt, xy, xytext=(0, 2), textcoords="offset points", ha="center", va="bottom",
                                     color=_ink(o["colour"]), zorder=6, path_effects=halo, **rkw))
                e = a.get_window_extent(r)
                w, h = e.width / pt, e.height / pt
                # above the peak; left or right of its top (under the place above it); higher; beside
                for dx, dy in ((0, 0), (-(w / 2 + 2), -(h + 0.5)), (w / 2 + 2, -(h + 0.5)), (0, h + 0.5),
                               (-(w + 1), 0), (w + 1, 0), (0, 2 * (h + 0.5))):
                    bx = (e.x0 + dx * pt, e.y0 + dy * pt, e.x1 + dx * pt, e.y1 + dy * pt)
                    if bx[0] >= bb.x0 and bx[2] <= bb.x1 and bx[3] <= bb.y1 + 0.5 * pt and free(bx):
                        if dx or dy:
                            a.xyann = (dx, 2 + dy)
                        taken.append(bx)
                        break
                else:  # no room: left out rather than on top of another label or outside the frame
                    a.remove()
                    self._label_arts.remove(a)

    def _labels_outside(self, ax, r, bb, x0, x1, pt, keep, taken, nkw):
        """Names right of the frame, each at the level of the end of its trace (median of the last 3 % of
        the view, not a spike of noise), in its colour, moved apart where they would touch; the right
        margin of the plot grows to hold the longest (as in a 2D waterfall figure)."""
        from matplotlib.transforms import blended_transform_factory
        tr_lab = blended_transform_factory(ax.transAxes, ax.transData)
        items = []
        for o in self.drawn:
            pr = o["proc"]
            xs = pr["t"] + o["dx"]
            m = (xs >= x0) & (xs <= x1)
            if not np.any(m):
                continue
            xe = float(xs[np.flatnonzero(m)[-1]])
            me = m & (xs >= xe - 0.03 * (x1 - x0))
            with np.errstate(all="ignore"):
                v = float(np.nanmedian(pr["y"][me])) if np.any(me) else float("nan")
            if not np.isfinite(v):
                continue
            yl = o["dy"] + v
            a = keep(ax.annotate(_tex(o["label"]), (1.0, yl), xycoords=tr_lab, xytext=(5, 0), textcoords="offset points",
                                 ha="left", va="center", color=_ink(o["colour"]), zorder=6, annotation_clip=False,
                                 **nkw))
            e = a.get_window_extent(r)
            # a long name takes at most 30 % of the width of the figure (a narrow image): on more lines,
            # then smaller (down to 5.5 pt), so that the plot keeps its room
            maxw = 0.3 * ax.figure.bbox.width
            if e.width > maxw and " " in o["label"].strip():
                import textwrap
                n = max(8, int(len(o["label"]) * maxw / e.width) + 2)
                a.set_text(_tex("\n".join(textwrap.wrap(o["label"], n))))
                e = a.get_window_extent(r)
            while e.width > maxw and a.get_fontsize() > 5.5:
                a.set_fontsize(max(5.5, a.get_fontsize() - 0.5))
                e = a.get_window_extent(r)
            yc = float(ax.transData.transform((x0, yl))[1])
            items.append([yc - e.height / 2.0, e.height, a, e, yc])
        if not items:
            return
        items.sort(key=lambda it: it[0])
        gap = 0.8 * pt
        for k in range(1, len(items)):
            items[k][0] = max(items[k][0], items[k - 1][0] + items[k - 1][1] + gap)
        lim_top, lim_bot = bb.y1 + 0.5 * items[-1][1], bb.y0 - 0.5 * items[0][1]
        if items[-1][0] + items[-1][1] > lim_top:  # within the height of the frame (half a line over at most)
            items[-1][0] = lim_top - items[-1][1]
            for k in range(len(items) - 2, -1, -1):
                items[k][0] = min(items[k][0], items[k + 1][0] - items[k][1] - gap)
        if items[0][0] < lim_bot:
            items[0][0] = lim_bot
            for k in range(1, len(items)):
                items[k][0] = max(items[k][0], items[k - 1][0] + items[k - 1][1] + gap)
        wmax = 0.0
        for y0, h, a, e, yc in items:
            a.xyann = (5, (y0 + h / 2.0 - yc) / pt)
            xa = bb.x1 + 5 * pt
            taken.append((xa, y0, xa + e.width, y0 + h))
            wmax = max(wmax, e.width)
        self.card.margins[2] = max(0.2, (wmax + 9 * pt) / ax.figure.dpi + 0.06)

    def _trace_key(self, o):
        """Key of a trace in the graph properties: its file."""
        return id(o["entry"]["doc"])

    def _labels_for_frame(self):
        """Labels placed for the tile, and once more when the room they take (the rows above the frame, a
        title on more lines, the names right of it) moved the frame: names placed for the frame before
        ran into each other."""
        ax = self.card.ax
        self._place_labels()
        pos = ax.get_position().bounds
        self.card.place()
        if max(abs(a - b) for a, b in zip(pos, ax.get_position().bounds)) > 1e-3:
            self._place_labels()

    def _relabel_resized(self):
        try:
            if self.IsShown() and self.drawn:
                self._labels_for_frame()
                self.card.draw()
        except RuntimeError:
            pass  # window closed

    def _table_note(self):
        return ""

    @staticmethod
    def _area_x_words(text):
        return K.split_x(text)

    def _bind_area_x(self, out):
        """The X values as typed (in the plotted order of out) go with their files ("xval" of each entry),
        so that they stay with their runs when the order changes. Values that do not match the traces
        one to one are left as typed (the table says what is missing)."""
        text = self.area_x.GetValue()
        words = self._area_x_words(text)
        self._area_x_shown = text
        if not words:
            for en in self.entries:
                en.pop("xval", None)
        self._area_x_bound = bool(words) and len(words) == len(out)
        if self._area_x_bound:
            for o, w in zip(out, words):
                o["entry"]["xval"] = w

    def _sync_area_x(self, out):
        """X values typed, restored or set elsewhere are bound to the files plotted now; once bound, the
        field follows the plotted order (sorted, moved, reversed, a trace left out or ticked again)."""
        text = self.area_x.GetValue()
        ents = [o["entry"] for o in out]
        if (text != getattr(self, "_area_x_shown", None) or not getattr(self, "_area_x_bound", False)
                or not any("xval" in en for en in ents)):
            self._bind_area_x(out)
            return
        words = [en.get("xval") for en in ents]
        if any(w is None for w in words):  # a run without a value (a file opened since): as typed
            return
        if words != self._area_x_words(text):  # (the text as typed while the order is the same)
            new = ", ".join(words)
            self.area_x.ChangeValue(new)
            self._area_x_shown = new

    def _area_x_values(self):
        """The X values typed (K.parse_x_values: commas, semicolons or new lines between them)."""
        return K.parse_x_values(self.area_x.GetValue())

    def _update_area_reference(self, out):
        """The 100 % reference: automatic, or the run chosen, kept by its file (it stays chosen while its
        trace is left out and is chosen again when it comes back)."""
        keys = [self._trace_key(o) for o in out]
        labels = ["largest area (automatic)"] + ["%d. %s" % (i + 1, o["label"]) for i, o in enumerate(out)]
        if list(self.area_ref.GetItems()) != labels:
            self.area_ref.SetItems(labels)
        self.area_ref.SetSelection(keys.index(self._area_ref_key) + 1 if self._area_ref_key in keys else 0)
        self._area_ref_keys = keys

    def on_area_ref(self, e=None):
        j = self.area_ref.GetSelection() - 1
        keys = getattr(self, "_area_ref_keys", [])
        self._area_ref_key = keys[j] if 0 <= j < len(keys) else None
        self.refresh_area()

    def on_compare_range(self, a, b, background=False):
        """The default drag measures every run; MS averaging is an explicit toggle."""
        self.set_area_range(a, b)

    def set_area_reference(self, o):
        self._area_ref_key = self._trace_key(o)
        self._update_area_reference(self.drawn)
        self.refresh_area()

    def _compute_area_results(self, out):
        self._update_area_reference(out)
        self._sync_area_x(out)
        self._area_x_error = None
        self._area_x_used = self.area_x.GetValue()
        self._area_xlabel_used = self.area_xlabel.GetValue()
        if self.area_range is None:
            return []
        try:
            xs = self._area_x_values()
            if xs and len(xs) != len(out):
                raise ValueError("Give one X value per trace (%d values for %d traces)" % (len(xs), len(out)))
        except ValueError as ex:  # the areas stay; X is left empty
            self._area_x_error = str(ex)
            xs = [None] * len(out)
        entries = [{"label": o["label"], "proc": o["proc"], "units": o["units"], "dx": o["dx"]}
                   for o in out]
        keys = self._area_ref_keys
        ref = keys.index(self._area_ref_key) if self._area_ref_key in keys else None
        return K.region_areas(entries, self.area_range[0], self.area_range[1], ref, xs)

    def _area_status(self):
        if self.area_range is None:
            return
        msg = "Integrated %.4f to %.4f min in %d traces" % (self.area_range[0], self.area_range[1],
                                                          len(self.area_results))
        self.status(msg + ("; X not used: %s" % self._area_x_error if self._area_x_error else ""))

    def set_area_range(self, a, b):
        try:
            a, b = sorted((float(a), float(b)))
        except (TypeError, ValueError):
            self.status("Give the start and end of the region in minutes")
            return
        if not np.isfinite([a, b]).all() or b <= a:
            self.status("The area region must have a start before its end")
            return
        self.area_range = (a, b)
        self.area_lo.SetValue("%.4f" % a)
        self.area_hi.SetValue("%.4f" % b)
        self.replot(keep_view=True, quick=False)
        if self.area_results:
            self._area_status()

    def on_area_fields(self, e=None):
        a, b = _fnum(self.area_lo.GetValue(), None), _fnum(self.area_hi.GetValue(), None)
        if a is None or b is None:
            self.status("Give both ends of the region")
            return
        self.set_area_range(a, b)

    def clear_area(self, e=None):
        self.area_range, self.area_results = None, []
        self.area_lo.SetValue("")
        self.area_hi.SetValue("")
        self.replot(keep_view=True, quick=False)
        self.status("Region cleared")

    def refresh_area(self):
        if self.area_range is not None:
            self.replot(keep_view=True, quick=False)
            if self.area_results:
                self._area_status()

    def _prepare_area_action(self):
        if self.area_range is None:
            self.status("Choose Area and drag over a peak first")
            return False
        # Buttons and exports must use the fields as typed, without requiring Enter.
        a, b = _fnum(self.area_lo.GetValue(), None), _fnum(self.area_hi.GetValue(), None)
        if a is None or b is None:
            self.status("Give both ends of the region")
            return False
        a, b = sorted((a, b))
        old = self.area_range
        # The fields show 4 decimals; a dragged region is unchanged when they still show it.
        if abs(a - old[0]) <= 6e-5 and abs(b - old[1]) <= 6e-5:
            self.refresh_area()
        else:
            self.set_area_range(a, b)
            if self.area_range == old:  # rejected; set_area_range said why
                return False
        if not self.area_results:
            self.status("No plotted traces to integrate")
            return False
        return True

    def _draw_area(self, ax, o, r, z):
        p = r.get("peak")
        if p is None:
            return
        pr = o["proc"]
        ts = p["t"] + o["dx"]
        factor = float(pr.get("factor", 1.0) or 1.0)
        ys = p["y"] * factor + o["dy"]
        base = p["baseline"] * factor + o["dy"]
        art = ax.fill_between(ts, base, np.maximum(ys, base), color=o["colour"], alpha=0.28, lw=0, zorder=z + 0.15)
        art._no_scale = art._area_fill = True
        ln, = ax.plot([ts[0], ts[-1]], [base[0], base[-1]], color=o["colour"], lw=0.8, alpha=0.9,
                      zorder=z + 0.2)
        ln._no_scale = True

    @staticmethod
    def _shift_of(o):
        """Time shift of a drawn trace (alignment plus the shift of its file), min."""
        return float(o.get("align", 0.0) or 0.0) + float(o["entry"].get("shift", 0.0) or 0.0)

    def _area_cells(self, out):
        """Columns and cells of the region areas, one row per trace of out. Peak times are retention times of
        each run (before its shift), as the main peak of the table."""
        units = sorted(set(r.get("units", "") for r in self.area_results))
        unit = units[0] if len(units) == 1 else ""
        cols = [("Region peak (min)", 110), (("Region height (%s)" % unit) if unit else "Region height", 120),
                (("Region area (%s·s)" % unit) if unit else "Region area", 130), ("Relative area (%)", 110),
                ("100 %", 50), (self._area_xlabel_used.strip() or "X", 90)]
        rows = []
        for i, o in enumerate(out):
            r = self.area_results[i] if i < len(self.area_results) else None
            if r is None:
                rows.append([""] * len(cols))
                continue
            rows.append(["%.3f" % (r["rt"] - self._shift_of(o)) if r["rt"] is not None else "",
                         K.number(r["height"]) if r["height"] is not None else "",
                         "%.9g" % r["area"] if r["area"] is not None else "no data",
                         "%.3f" % r["relative"] if r["relative"] is not None else
                         ("" if r["area"] is None else "undefined"),
                         "yes" if r["reference"] else "", "%.9g" % r["x"] if r["x"] is not None else ""])
        return cols, rows

    def _area_info(self):
        """What the table says about the region: its ends, the reference, X values not used."""
        if self.area_range is None:
            return []
        a, b = self.area_range
        out = ["Region %.4f to %.4f min, 100 %%: %s" % (
            a, b, next((r["label"] for r in self.area_results if r["reference"]), "none"))]
        if self.area_results and all(r["relative"] is None for r in self.area_results):
            out.append("the reference has no area: choose another run")
        if self._area_x_error:
            out.append("X not used: %s" % self._area_x_error)
        return out

    def _kinetic_problem(self, message):
        self.status(message)
        wx.MessageBox(message, "Kinetic fitting", wx.OK | wx.ICON_INFORMATION,
                      parent=wx.GetTopLevelParent(self))
        return None

    def open_kinetics(self, e=None):
        if self.area_range is None:
            return self._kinetic_problem("Drag over a peak with the Area tool first.")
        if not self.drawn:
            return self._kinetic_problem("Include at least three traces.")
        if not self.area_x.GetValue().strip():
            if len(self.drawn) < 3:
                return self._kinetic_problem("Include at least three traces.")
            order = "\n".join("%d. %s" % (i + 1, o["label"]) for i, o in enumerate(self.drawn))
            with wx.TextEntryDialog(wx.GetTopLevelParent(self),
                    "One X value per trace, in this order:\n" + order, "X values") as dlg:
                if dlg.ShowModal() != wx.ID_OK:
                    self.status("Kinetic fitting cancelled")
                    return None
                text = dlg.GetValue().strip()
            # Validate before changing the user's existing selected region or areas.
            try:
                values = K.parse_x_values(text)
            except ValueError as ex:
                return self._kinetic_problem(str(ex))
            if len(values) != len(self.drawn):
                return self._kinetic_problem("Give one X value per trace (%d values for %d traces)" %
                                             (len(values), len(self.drawn)))
            if len(set(values)) < 2:
                return self._kinetic_problem("Give at least two different X values.")
            self.area_x.ChangeValue(text)
        if not self._prepare_area_action():
            return self._kinetic_problem("No areas in the region.")
        if self._area_x_error:
            return self._kinetic_problem(self._area_x_error)
        relative = self.area_y.GetSelection() == 1
        if any(r.get("peak") is None for r in self.area_results):
            return self._kinetic_problem("Some traces have no data in the region.")
        if relative and any(r["relative"] is None for r in self.area_results):
            return self._kinetic_problem("The 100 % reference has no area.")
        if len(set(r.get("units", "") for r in self.area_results)) > 1:
            return self._kinetic_problem("The traces have different units.")
        if len(self.area_results) < 3 or len(set(r["x"] for r in self.area_results)) < 2:
            return self._kinetic_problem("Give at least three traces and two different X values.")
        try:
            from ms_kinetics_ui import KineticsFrame
            metadata = [("Signal", self.card.subtitle), ("Region start (min)", self.area_range[0]),
                        ("Region end (min)", self.area_range[1]),
                        ("100% reference", next((r["label"] for r in self.area_results if r["reference"]), "none"))]
            rows = [dict(r, file=getattr(self.drawn[r["index"]]["entry"]["doc"], "path", ""))
                    for r in self.area_results]
            fr = KineticsFrame(self, rows, relative, metadata)
        except Exception as ex:
            import traceback
            traceback.print_exc()
            return self._kinetic_problem("Cannot open kinetic fitting: %s" % ex)
        self._area_frames.append(fr)
        fr.Bind(wx.EVT_CLOSE, lambda ev, f=fr: (self._area_frames.remove(f) if f in self._area_frames else None,
                                              ev.Skip()))
        fr.Show()
        fr.Raise()
        return fr

    def plot_area_graph(self, e=None):
        if not self._prepare_area_action():
            return None
        if self._area_x_error:
            self.status("X not used: %s" % self._area_x_error)
            return None
        relative = self.area_y.GetSelection() == 1
        if relative and not any(r["reference"] for r in self.area_results):
            self.status("The 100 % reference has no area")
            return
        pts = sorted([(float(r["x"]), float(r["relative"] if relative else r["area"]), r["label"])
                      for r in self.area_results if r["area"] is not None], key=lambda q: q[0])
        fr = wx.Frame(self.frame, title="Peak area vs X", size=wx.Size(self.FromDIP(760), self.FromDIP(500)))
        card = U.PlotCard(fr, "Peak area vs X", "region %.4f to %.4f min" % self.area_range, mode="zoom")
        x = np.asarray([p[0] for p in pts], float)
        y = np.asarray([p[1] for p in pts], float)
        U.style_axes(card.ax, self.area_xlabel.GetValue().strip() or "X",
                     "Relative area (%)" if relative else "Actual area")
        card.ax.plot(x, y, color=C["accent"], lw=1.2, marker="o", ms=5, zorder=3)
        for xx, yy, label in pts:
            card.ax.annotate(label, (xx, yy), xytext=(0, 6), textcoords="offset points", ha="center",
                             va="bottom", fontsize=7, color=U.INK)
        if len(x):
            pad = max((float(np.max(x)) - float(np.min(x))) * 0.06, 0.5 if len(x) == 1 else 1e-6)
            card.full = (float(np.min(x)) - pad, float(np.max(x)) + pad)
            card.ax.set_xlim(*card.full)
        card.autoscale_y(pad=0.16)
        card.image_name = lambda: self.image_name() + "_area_vs_x"
        card.on_menu = lambda x_, y_: [("Kinetic fitting…", self.open_kinetics),
                                      ("Export the areas to Excel…", self.on_export_areas_xlsx)]
        s = wx.BoxSizer(wx.VERTICAL)
        s.Add(card, 1, wx.EXPAND | wx.ALL, self.FromDIP(10))
        fr.SetSizer(s)
        fr.Bind(wx.EVT_CLOSE, lambda ev, f=fr: (self._area_frames.remove(f) if f in self._area_frames else None,
                                                ev.Skip()))
        self._area_frames.append(fr)
        fr.Show()
        card.draw()
        return fr

    def _fill_table(self, out, notes):
        """The peaks of each trace (main peak, area % at the guide lines) and, with a region set, its areas in
        more columns; the file details at the end."""
        s = self.s
        ents = [{"label": o["label"], "proc": o["proc"], "units": o["units"], "shift": self._shift_of(o)}
                for o in out]
        hdr, rows, _pk = K.peak_rows(ents, s["guides"], s.get("integ_thr", 1.0))
        # the results first (they were out of sight to the right of the file details)
        widths = [160, 105, 100, 115, 68] + [110] * max(0, len(hdr) - 5)
        pcols = [(h, widths[i]) for i, h in enumerate(hdr)]
        acols, arows = self._area_cells(out) if self.area_range is not None else ([], [[] for _ in out])
        fcols = [("Shift (min)", 80), ("File", 150), ("Sample", 130), ("Acquired", 120)]
        cols = pcols + acols + fcols
        full, peaks = [], []
        for o, r, ar, en in zip(out, rows, arows, ents):
            d = o["entry"]["doc"]
            si = self.info(d)
            acq = si.get("acquired")
            acq = acq.strftime("%d.%m.%Y %H:%M") if isinstance(acq, datetime.datetime) else ""
            sh = en["shift"]
            det = ["%+.3f" % sh if abs(sh) > 1e-9 else "0", os.path.basename(d.attrs.get("path") or ""),
                   si.get("sample_name", ""), acq]
            full.append(list(r) + list(ar) + det)
            peaks.append(list(r) + det)
        if [c[0] for c in self.table.cols] != [c[0] for c in cols]:
            self.table.set_columns(cols)
        self.table.title = "Peaks of each trace"
        info = self._area_info() + list(notes)
        self.table.set_rows(full, "; ".join(info) if info else self._table_note())
        # the table shown (and exported); the peaks and the region areas also apart (the report)
        self._table = ([c[0] for c in cols], full)
        self._peak_table = ([c[0] for c in pcols + fcols], peaks)
        self._area_table = (["Trace"] + [c[0] for c in acols], [[o["label"]] + list(ar) for o, ar in
                                                                zip(out, arows)]) if acols else None

    # ------------------------------------------------------------------ mouse
    def _nearest(self, x, y):
        best, bd = None, None
        for o in self.drawn:
            pr = o["proc"]
            xs = pr["t"] + o["dx"]
            if not (xs[0] <= x <= xs[-1]):
                continue
            v = float(np.interp(x, xs, pr["y"])) + o["dy"]
            dd = abs(v - y)
            if not np.isfinite(dd):  # no value there (a gap): it would win every comparison
                continue
            if bd is None or dd < bd:
                best, bd = o, dd
        return best

    def readout(self, x, y):
        o = self._nearest(x, y)
        if o is None:
            return "%.3f min" % x
        pr = o["proc"]
        v = float(np.interp(x, pr["t"] + o["dx"], pr["y"]))
        unit = "%" if o["scaled"] else o["units"]
        return "%s   %.3f min   %s %s" % (o["label"], x - o["dx"], U._g(v), unit)

    def on_pick(self, x, y, alt=None):
        o = self._nearest(x, y)
        if self.mz_on:  # m/z tool: the peak of that trace there (Alt + click: one scan at that time), and its spectra
            if o is None:
                self.status("m/z: click on a trace")
                return
            if alt is None:
                alt = self._alt_down()
            if alt:
                self.pick_scan(o, x)
            else:
                self.pick_peak(o, x)
            return
        if o is not None:
            self.choose(o["i"])

    def _on_press(self, e):
        """Where a left button press on the plot was (and its keys): the trace a drag starts on, Alt + click."""
        try:
            if e.button == 1 and e.inaxes is self.card.ax and e.xdata is not None and e.ydata is not None:
                self._mz_press = (float(e.xdata), float(e.ydata), str(e.key or ""))
        except Exception:
            self._mz_press = None

    def _alt_down(self):
        if self._mz_press is not None and "alt" in self._mz_press[2]:
            return True
        try:
            return bool(wx.GetKeyState(wx.WXK_ALT))
        except Exception:
            return False

    def _trace_over(self, a, b, y):
        """The trace drawn across the plot times a to b nearest to y in the middle of its part there (a drag that
        starts outside the time range of every trace)."""
        best, bd = None, None
        for o in self.drawn:
            xs = o["proc"]["t"] + o["dx"]
            if not len(xs):
                continue
            lo, hi = max(a, float(xs[0])), min(b, float(xs[-1]))
            if hi < lo:
                continue
            v = float(np.interp(0.5 * (lo + hi), xs, o["proc"]["y"])) + o["dy"]
            dd = abs(v - y) if (y is not None and np.isfinite(v)) else 0.0
            if bd is None or dd < bd:
                best, bd = o, dd
        return best

    def on_mz_range(self, a, b, shift):
        """A drag across the plot while m/z is on (Select tool): the time range of the trace it starts on is
        averaged for its spectra; Shift + drag sets the background range of that run."""
        if not self.mz_on:
            self.card.canvas.draw_idle()
            return
        press = self._mz_press
        o = None
        if press is not None and a - 1e-9 <= press[0] <= b + 1e-9:
            o = self._nearest(press[0], press[1])
        if o is None:
            o = self._trace_over(a, b, press[1] if press is not None else None)
        if o is None:
            self.status("m/z: drag across a trace")
            self.card.canvas.draw_idle()
            return
        if shift:
            self.set_mz_bg(o, a, b)
        else:
            self.pick_range(o, a, b)

    def on_tool(self, tool, kind, x, y):
        if tool == "area" and kind == "drag" and x is not None and y is not None:
            self.set_area_range(x, y)
            return
        if kind != "click" or x is None:
            return
        if tool == "align":
            self.set_align(x, y)
        elif tool == "guide":
            self.toggle_guide(x)

    def _run_time(self, x, y=None):
        """The time in the runs of a point of the plot: x less the skew of the trace there (Offset layout:
        the traces behind the front one are drawn moved right). y picks the trace; without one, or away
        from every trace, the front trace (no skew)."""
        dxs = set(float(o["dx"]) for o in self.drawn)
        if len(dxs) == 1:
            return x - dxs.pop()
        o = self._nearest(x, y) if y is not None else None
        return x - (o["dx"] if o is not None else 0.0)

    def set_align(self, x, y=None):
        if x is not None:
            x = self._run_time(x, y)
        self.align.SetValue("" if x is None else "%.3f" % x)
        self.on_change()

    def toggle_guide(self, x):
        g = list(self.s["guides"])
        x0, x1 = self.card.ax.get_xlim()
        near = [v for v in g if abs(v - x) <= (x1 - x0) / 80.0]
        if near:
            g.remove(near[0])
        else:
            g.append(round(float(x), 2))
        self.guides.SetValue(", ".join("%g" % v for v in sorted(g)))
        self.on_change()

    def plot_menu(self, x, y):
        items = []
        if not self.drawn:
            return items
        x0, x1 = self.card.ax.get_xlim()
        g = self.s["guides"]
        near = [v for v in g if abs(v - x) <= (x1 - x0) / 80.0]
        items.append(("Guide line at %.2f min" % x, lambda: self.toggle_guide(x)) if not near else
                     ("Remove the guide line at %.2f min" % near[0], lambda: self.toggle_guide(near[0])))
        if g:
            items.append(("Remove every guide line", lambda: (self.guides.SetValue(""), self.on_change())))
        items.append((None, None))
        tx = self._run_time(x, y)  # the time in the runs there (Offset: less the skew of that trace)
        items.append(("Align every trace on its peak near %.2f min" % tx, lambda: self.set_align(x, y)))
        if self.s["align"] is not None:
            items.append(("Undo the alignment", lambda: self.set_align(None)))
        items.append(("Scale: each trace at 100 %% at its peak near %.2f min" % tx, lambda: self.set_ref(x, y)))
        if self.area_range is not None:
            items += [("Area vs X…", self.plot_area_graph),
                      ("Kinetic fitting…", self.open_kinetics),
                      ("Areas to Excel…", self.on_export_areas_xlsx),
                      ("Clear the region", self.clear_area)]
        if self.s["t0"] is None and self.s["t1"] is None:
            items.append(("Use the time range shown (%.2f to %.2f min)" % (x0, x1), lambda: self.set_window(x0, x1)))
        else:
            items.append(("The whole run", lambda: self.set_window(None, None)))
        o = self._nearest(x, y)
        if self.mz_on:  # the background ranges of the m/z tool (Shift + drag across a trace)
            mine = o["entry"].get("mzbg") if o is not None else None
            n_bg = sum(1 for en in self.entries if en.get("mzbg"))
            if mine or n_bg:
                items.append((None, None))
            if mine:
                items.append(("Remove the background range of %s (%.2f to %.2f min)" % ((o["label"],) + tuple(mine)),
                              lambda en=o["entry"]: self.clear_mz_bg(en)))
            if n_bg > (1 if mine else 0):
                items.append(("Remove every background range", lambda: self.clear_mz_bg(None)))
        if o is not None:
            items.append((None, None))
            if self.area_range is not None:
                items.append(("Use %s as 100 %% area reference" % o["label"], lambda o=o: self.set_area_reference(o)))
            items.append(("Leave out %s" % o["label"], lambda o=o: self.leave_out(o["i"])))
        items += [(None, None), ("Graph properties\u2026", self.on_graph_props),
                  (None, None), ("Export the traces (CSV)…", self.on_export_traces),
                  ("Export the table (CSV)…", self.on_export_table)]
        return items + self.view_menu()

    def on_graph_props(self):
        """Graph properties of this plot (fonts, axes, frame, the line of each trace, title): kept for the
        session, also in copied and saved images and the report."""
        traces = []
        for o in self.drawn:
            k = self._trace_key(o)
            tp = self.props["traces"].get(k, {})
            traces.append((k, o["label"], o["colour"], tp.get("lw"), tp.get("ls") or "-"))
        saved_colours = {id(en["doc"]): en.get("colour") for en in self.entries}
        shown_colours = {t[0]: (t[2] or "").lower() for t in traces}

        def apply(p, colours, common):
            old = self.props
            self.props = p
            for en in self.entries:
                k = id(en["doc"])
                if colours is None:  # Cancel: the colours of the moment the dialog opened
                    if k in saved_colours:
                        en["colour"] = saved_colours[k]
                elif k in colours:
                    # back to the colour it had when the dialog opened (Reset all): also back to automatic
                    en["colour"] = (saved_colours.get(k) if colours[k].lower() == shown_colours.get(k)
                                    else colours[k])
            if common is not None and abs(float(common) - float(self.s["lw"])) > 1e-9:
                self.s["lw"] = float(common)
                self._busy = True
                try:
                    self.lw.SetValue("%g" % self.s["lw"])
                finally:
                    self._busy = False
                self._save()
            self._fill_entry_fields()
            self.card.ylock = None
            self.replot(keep_view=(old.get("x0"), old.get("x1")) == (p.get("x0"), p.get("x1")))

        ax_t = getattr(self, "_auto_titles", ("Time (min)", self.ylabel))
        dlg = PP.GraphPropsDialog(self, self.props, traces, self.s["lw"], apply, ax_t[0], ax_t[1])
        try:
            dlg.ShowModal()
        finally:
            dlg.Destroy()

    def set_ref(self, x, y=None):
        self.ref_t.SetValue("%.3f" % self._run_time(x, y))
        self.scale.SetSelection(_index(K.SCALES, "ref"))
        self.on_change()

    def set_window(self, a, b):
        self.t0.SetValue("" if a is None else "%.3f" % a)
        self.t1.SetValue("" if b is None else "%.3f" % b)
        self.on_change()

    def leave_out(self, i):
        if 0 <= i < len(self.entries):
            self.entries[i]["include"] = False
            self.flist.Check(i, False)
            self.replot(keep_view=True)

    def on_lambda_max(self, e=None):
        """The absorption maximum (above 210 nm) of the tallest peak of the first file ticked: the peak
        of the PDA trace shown (max plot for the other signals), after the first 8 % of the run (the
        injection peak), its UV spectrum minus the spectrum just before the peak."""
        for en in self.entries:
            if not en["include"] or en["doc"] is self.blank_doc:
                continue
            pda = en["doc"].attrs.get("pda_data")
            if pda is None:
                continue
            # the time window is in plotted time: the run's alignment shift as well as its own time shift
            align = next((o["align"] for o in self.drawn if o["entry"] is en), 0.0)
            r = K.lambda_max(pda, self.s, float(align) + float(en.get("shift", 0.0) or 0.0))
            if r is None:
                continue
            lm, ta, bg = r
            self.wl.SetValue("%.0f" % lm)
            if self.s["signal"] != "pda":
                self.signal.SetSelection(_index(K.SIGNALS, "pda"))
            self.on_change()
            self.status("\u03bbmax %.0f nm (%s, %.2f min)" % (lm, self.label_of(en), ta))
            return
        self.status("No file with PDA data ticked")

    # ------------------------------------------------------------------ m/z tool: spectra of a peak
    @staticmethod
    def _pol_name(pol):
        return "ESI+" if pol == "+" else "ESI−"

    def set_mz(self, on):
        """The m/z tool on (the row of the two spectrum tiles below the plot; a click on a trace picks its
        peak) or off (row, picked peak and its mark gone; a click chooses the file again)."""
        self.mz_on = bool(on)
        b = getattr(self, "tools", None) and self.tools.buttons.get("mzspec")
        if b:
            b.set_active(self.mz_on)
        # while it is on, a drag across a trace averages that time range of its run (Shift: background), as on the
        # chromatograms of the Mass spectrometry view; Ctrl + drag zooms; off: a drag integrates all traces
        self.card.mode = "range"
        self.card.on_range = self.on_mz_range if self.mz_on else self.on_compare_range
        if not self.mz_on:
            self.mzpick = None
            for en in self.entries:  # the background ranges go with the tiles
                en.pop("mzbg", None)
            for c in self.mzcols.values():
                c.update(measure=None, m1=None)
        self._show_spec_tools()
        self._mz_layout()
        for p in ("+", "-"):
            self.plot_mz(p)
        self.replot(keep_view=True, quick=True)
        self.status(self.TOOL_HINTS["mzspec"] if self.mz_on else "m/z off")

    def _mz_layout(self):
        rows = self.rows
        wc, wt = rows.weight_of(self.card, 5.0), rows.weight_of(self.table, 1.5)
        if self.mz_on:
            rows.closed.discard(self.mzrow)  # tiles closed before (right click > Close this tile) come back
            self.mzrow.closed.clear()
            self.mzrow.set_panes([self.mzcards["+"], self.mzcards["-"]],
                                 self.mzrow.all_weights if len(self.mzrow.all_weights) == 2 else [1, 1])
            if rows.frac and self.mzrow not in rows.frac:  # tiles sized by the user: the row below at its share
                rows.frac[self.mzrow] = 0.4
            rows.set_panes([self.card, self.mzrow, self.table], [wc, self._mz_weight, wt])
        else:
            self._mz_weight = rows.weight_of(self.mzrow, self._mz_weight)
            rows.set_panes([self.card, self.table], [wc, wt])
            self.mzrow.Hide()
        self.scroller.refit()

    def _pick_alive(self):
        """The picked peak is dropped when its file was closed or left the list."""
        pk = self.mzpick
        if pk is None:
            return
        docs = getattr(self.frame, "docs", None)
        d = pk["doc"]
        if any(en["doc"] is d for en in self.entries) and (docs is None or d in docs) and d.attrs.get("path"):
            return
        self.mzpick = None
        for p in ("+", "-"):
            self.plot_mz(p)

    def _delay(self, d):
        """MS detector delay of a file (its PDA view): MS time = PDA time + delay."""
        try:
            v = float(U._num(d.attrs.get("pda").delay, 0.0) or 0.0)
            return v if np.isfinite(v) else 0.0  # 'nan' or 'inf' typed: as no delay (the field reads them as numbers)
        except Exception:
            return 0.0

    def _binw(self, d):
        """Bin width of the spectra of the file's MS view (m/z)."""
        try:
            ms_tab = d.attrs.get("ms")
            v = float(ms_tab.bin_width()) if ms_tab is not None and hasattr(ms_tab, "bin_width") else 0.05
            return v if v > 0 else 0.05
        except Exception:
            return 0.05

    def _trace_times(self, o):
        """(shift, pda, file delay, delay) of the drawn trace o: plot time - dx - shift = time of its run; MS time =
        that + delay (PDA traces: the MS detector delay of the file)."""
        en = o["entry"]
        sh = float(o.get("align", 0.0) or 0.0) + float(en.get("shift", 0.0) or 0.0)
        pda = self.s["signal"] in ("pda", "max")
        file_delay = self._delay(en["doc"])
        return sh, pda, file_delay, (file_delay if pda else 0.0)

    def _busy_ms(self, o):
        ms_tab = o["entry"]["doc"].attrs.get("ms")
        if getattr(ms_tab, "_averaging", False):
            self.status("The Mass spectrometry view of %s is still averaging; try again" % o["label"])
            return True
        return False

    def _new_pick(self, o, kind, ta, a, b):
        """The pick (kind "peak", "none", "range" or "scan") of trace o at the times ta (top or middle), a to b of
        the processed trace; its spectra computed and drawn, its file chosen (the plot drawn with the mark)."""
        t_start = time.perf_counter()
        en = o["entry"]
        sh, pda, file_delay, delay = self._trace_times(o)
        off = delay - sh  # processed trace time -> MS time of the run
        self.mzpick = {"doc": en["doc"], "entry": en, "label": o["label"], "kind": kind, "pda": pda, "delay": delay,
                       "ms_apex": ta + off, "ms_t0": a + off, "ms_t1": b + off, "trace_apex": ta - sh,
                       "signal": K.signal_text(self.s), "spectra": {}, "masses": [], "file_delay": file_delay}
        self._mz_compute()
        for q in ("+", "-"):
            self.plot_mz(q)
        self.choose(o["i"])  # its file chosen in the list; the plot drawn again with the mark
        pk = self.mzpick
        pk["seconds"] = time.perf_counter() - t_start
        return pk

    def _pick_status(self, pk, msg):
        sug = K.mass_suggestion(pk.get("masses"), more=0, info=pk.get("mass_info"))
        self.status(msg + ("; " + sug if sug else ""))
        U._log("compare m/z: %s, kind %s, %.3f s" % (pk["label"], pk["kind"], pk["seconds"]))

    def pick_peak(self, o, x):
        """The peak of the drawn trace o nearest to the time x clicked (plot time): its top within 0.15 min
        and its start and end from the library, on the trace as processed; on the baseline, the time
        clicked. Times of the run itself (without the skew, alignment and shift of the trace), MS times
        for the spectra (PDA traces: plus the MS detector delay of the file)."""
        pr = o["proc"]
        t, y = pr["t"], pr["y"]
        ok = np.isfinite(t) & np.isfinite(y)
        if not ok.all():  # points without a value (a gap): the peak search runs on the others
            t, y = t[ok], y[ok]
        tp = float(x) - float(o.get("dx", 0.0) or 0.0)  # time on the processed (shifted) trace
        if self._busy_ms(o):
            return
        apex, p = None, None
        try:
            apex = K.find_apex(t, y, tp, 0.15)
            if apex is not None:
                p = LI.peak_at(t, y, apex, search_s=3.0)
            else:
                p = LI.peak_at(t, y, tp, search_s=9.0)
                if p is not None:
                    apex = float(p["rt"])
            if p is not None and not (p["t0"] - 1e-9 <= apex <= p["t1"] + 1e-9):
                p = None
        except Exception as ex:
            U._log("compare m/z: peak at %.3f min: %s" % (tp, ex))
            apex, p = None, None
        if p is not None:
            kind, ta, a, b = "peak", float(apex), float(p["t0"]), float(p["t1"])
        else:
            kind = "none"
            ta = float(apex) if apex is not None else tp
            a = b = ta
        pk = self._new_pick(o, kind, ta, a, b)
        delay = pk["delay"]
        dl = " (PDA time %.2f min + MS delay %.2f min)" % (pk["trace_apex"], delay) if pk["pda"] else ""
        if kind == "peak":
            msg = "m/z: peak of %s at %.2f min%s, MS scans of %.2f to %.2f min" % (
                o["label"], pk["ms_apex"], dl, pk["ms_t0"], pk["ms_t1"])
        else:
            msg = "m/z: no peak of %s at %.2f min%s: the scans around that time" % (o["label"], pk["ms_apex"], dl)
        self._pick_status(pk, msg)

    def _clip_to(self, o, a, b):
        """Plot times a to b held within the time range of the drawn trace o (a drag past the end of the run), as
        times of the processed trace; None when they miss it."""
        t = o["proc"]["t"]
        t = t[np.isfinite(t)]
        if not len(t) or not (np.isfinite(a) and np.isfinite(b)):
            return None
        dx = float(o.get("dx", 0.0) or 0.0)
        lo, hi = max(min(a, b) - dx, float(t[0])), min(max(a, b) - dx, float(t[-1]))
        if hi < lo:
            return None
        return lo, hi

    def pick_range(self, o, a, b):
        """The time range a to b (plot times) of the drawn trace o averaged for both polarities (a drag across
        the trace while m/z is on), minus the background range of that run if one is set (Shift + drag)."""
        r = self._clip_to(o, a, b)
        if r is None:
            self.status("m/z: the range %.2f to %.2f min misses the trace of %s" % (a, b, o["label"]))
            self.card.canvas.draw_idle()
            return
        if self._busy_ms(o):
            self.card.canvas.draw_idle()
            return
        lo, hi = r
        pk = self._new_pick(o, "range", 0.5 * (lo + hi), lo, hi)
        sp = next((pk["spectra"][q] for q in ("+", "-") if pk["spectra"].get(q, {}).get("range")), None)
        dl = " (PDA time + MS delay %.2f min)" % pk["delay"] if pk["pda"] else ""
        msg = "m/z: %s averaged %.2f to %.2f min%s" % (o["label"], pk["ms_t0"], pk["ms_t1"], dl)
        if sp is not None and sp.get("bg") and sp.get("bg_manual"):
            msg += ", background %.2f to %.2f min subtracted" % tuple(sp["bg"])
        self._pick_status(pk, msg)

    def pick_scan(self, o, x):
        """The scan of each polarity nearest to the time x (plot time) of the drawn trace o (Alt + click while m/z
        is on), as a click on a chromatogram of the Mass spectrometry view; minus the background range of the run if
        one is set."""
        r = self._clip_to(o, x, x)
        if r is None:
            self.status("m/z: click on a trace")
            return
        if self._busy_ms(o):
            return
        pk = self._new_pick(o, "scan", r[0], r[0], r[0])
        dl = " (PDA time %.2f min + MS delay %.2f min)" % (pk["trace_apex"], pk["delay"]) if pk["pda"] else ""
        self._pick_status(pk, "m/z: the scans of %s at %.2f min%s" % (o["label"], pk["ms_apex"], dl))

    def set_mz_bg(self, o, a, b):
        """The background range of the run of the drawn trace o (Shift + drag across it while m/z is on, plot times
        a to b; kept as MS times): subtracted from the spectra of that run, in place of the spectrum just before a
        peak clicked; drawn grey on its trace."""
        r = self._clip_to(o, a, b)
        if r is None or r[1] - r[0] <= 0:
            self.status("m/z: the background range %.2f to %.2f min misses the trace of %s" % (a, b, o["label"]))
            self.card.canvas.draw_idle()
            return
        sh, pda, file_delay, delay = self._trace_times(o)
        en = o["entry"]
        en["mzbg"] = (r[0] - sh + delay, r[1] - sh + delay)
        msg = "m/z: background of %s %.2f to %.2f min%s, subtracted from its spectra" % (
            o["label"], en["mzbg"][0], en["mzbg"][1], " (PDA time + MS delay %.2f min)" % delay if pda else "")
        pk = self.mzpick
        if pk is not None and pk["entry"] is en and not self._busy_ms(o):
            self._mz_compute()  # the spectra shown are those of this run: again with this background
            for q in ("+", "-"):
                self.plot_mz(q, keep_view=True)
        self.replot(keep_view=True, quick=True)
        self.status(msg)

    def clear_mz_bg(self, en=None):
        """The background range of one run (or of every run) removed; spectra of that run computed again."""
        for e in (self.entries if en is None else [en]):
            e.pop("mzbg", None)
        pk = self.mzpick
        if pk is not None and (en is None or pk["entry"] is en):
            self._mz_compute()
            for q in ("+", "-"):
                self.plot_mz(q, keep_view=True)
        self.replot(keep_view=True, quick=True)

    def _mz_compute(self):
        """The averaged spectra of both polarities of the picked peak (C++ ms_average), and the neutral
        masses they suggest."""
        pk = self.mzpick
        if not pk:
            return
        d = pk["doc"]
        ms = d.attrs.get("ms_data")
        binw = self._binw(d)
        pk["binw"] = binw
        pk["spectra"] = {}
        for c in self.mzcols.values():  # new spectra: the measurement goes (the labels pinned stay), as in the MS view
            c.update(measure=None, m1=None)
            if c.get("pin_doc") is not d:  # labels pinned on the spectra of another run: they belong to that run
                c["pins"] = []
            c["pin_doc"] = d
        for pol in ("+", "-"):
            try:
                pk["spectra"][pol] = self._mz_spectrum(ms, pol, pk, binw)
            except Exception as ex:
                U._log("compare m/z: %s spectrum: %s" % (pol, ex))
                pk["spectra"][pol] = {"spec": None, "n": 0, "range": None, "bg": None, "note": "", "e": None,
                                      "msg": "The %s spectrum of %s could not be read (%s)" % (
                                          self._pol_name(pol), pk["label"], ex)}
        both = []
        for pol in ("+", "-"):
            s_ = pk["spectra"][pol].get("spec")
            both.append((s_[:, 0], s_[:, 1]) if s_ is not None and len(s_) else None)
        masses = []
        pk["mass_info"] = None
        if both[0] is not None or both[1] is not None:
            try:
                import ms_adducts
                # min_rel 4 %: the [M+HCOO]- or [M+Na]+ ion of a compound is often a few % of the base peak
                if hasattr(ms_adducts, "neutral_masses_info"):
                    info = ms_adducts.neutral_masses_info(both[0], both[1], tol=0.3, min_rel=0.04, top=8) or {}
                    masses = info.get("masses") or []
                    pk["mass_info"] = {"base": info.get("base") or {}, "explained": info.get("explained") or {}}
                else:
                    masses = ms_adducts.neutral_masses(both[0], both[1], tol=0.3, min_rel=0.04) or []
                # (finite masses only, as in the text of the suggestion)
                masses = [m for m in masses if isinstance(m, dict) and isinstance(m.get("mass"), (int, float)) and
                          np.isfinite(m["mass"])]
            except Exception as ex:
                U._log("compare m/z: suggested masses: %s" % ex)
                masses = []
        pk["masses"] = masses

    def _mz_spectrum(self, ms, pol, pk, binw):
        """The spectrum of one polarity: {"spec" (m/z, intensity) or None, "n" scans, "range" (MS times),
        "bg" range subtracted or None, "note", "msg" (why there is none), "e" scan event}."""
        sp = {"spec": None, "n": 0, "range": None, "bg": None, "bg_manual": False, "note": "", "msg": "", "e": None}
        name = pk["label"]
        word = "positive" if pol == "+" else "negative"
        if ms is None:
            sp["msg"] = "%s has no MS data" % name
            return sp
        e = K._event_of(ms, pol)
        if e is None:
            sp["msg"] = "%s has no %s ion scans" % (name, word)
            return sp
        sp["e"] = e
        rt = np.asarray(ms.rt[ms.event_scans(e)], float)
        if not len(rt):
            sp["msg"] = "%s has no %s ion scans" % (name, word)
            return sp
        k = int(np.argmin(np.abs(rt - pk["ms_apex"])))  # the scan nearest to the top
        # a time outside the MS run (a large MS detector delay): no spectrum (the scan nearest to it, minutes away,
        # is not one of this peak)
        gap = max(2.0 * float(np.median(np.diff(rt))) if len(rt) > 1 else 0.0, 0.02)
        if not abs(float(rt[k]) - pk["ms_apex"]) <= gap and not np.any((rt >= pk["ms_t0"] - gap) &
                                                                        (rt <= pk["ms_t1"] + gap)):
            at = ("%.2f to %.2f min" % (pk["ms_t0"], pk["ms_t1"]) if pk["kind"] == "range" else
                  "%.2f min" % pk["ms_apex"])
            sp["msg"] = "%s has no %s ion scans at %s (MS scans from %.2f to %.2f min)" % (
                name, word, at, float(rt[0]), float(rt[-1]))
            return sp
        notes = []
        kind = pk["kind"]
        if kind == "none":  # no peak: the scan at that time and its neighbours
            a, b = float(rt[max(k - 1, 0)]), float(rt[min(k + 1, len(rt) - 1)])
        elif kind == "scan" or (kind == "peak" and self.s.get("spec_avg") == "scan"):
            a = b = float(rt[k])
        else:  # the peak from start to end, or the range dragged
            a, b = pk["ms_t0"], pk["ms_t1"]
            if not np.any((rt >= a) & (rt <= b)):
                a = b = float(rt[k])
                notes.append("no scan within the range: the scan nearest to it" if kind == "range" else
                             "no scan within the peak: the scan nearest to its top")
        bg = None
        manual = False
        own = pk["entry"].get("mzbg")  # the background range of the run (Shift + drag), MS times
        if own and (kind in ("range", "scan") or self.s.get("spec_bg")):
            # it replaces the spectrum just before a peak clicked (the option of the side panel still decides for
            # clicks); ranges and single scans always have it subtracted
            c, dd = own
            if np.any((rt >= c) & (rt <= dd)):
                bg, manual = (c, dd), True
            else:
                notes.append("no scans in the background range %.2f to %.2f min (nothing subtracted)" % (c, dd))
        elif self.s.get("spec_bg") and kind == "peak":
            # the spectrum just before the peak: the scans of a window as long as the peak (at least 0.05 min)
            # that ends where the peak starts
            w = max(pk["ms_t1"] - pk["ms_t0"], 0.05)
            c, dd = max(pk["ms_t0"] - w, float(rt[0])), pk["ms_t0"]
            if dd > c and np.any((rt >= c) & (rt <= dd)):
                bg = (c, dd)
            else:
                notes.append("no scans before the peak (nothing subtracted)")
        spec, n = ms.average(e, a, b, bg=bg, binw=binw)
        sp.update(range=(a, b), bg=bg, bg_manual=manual, n=int(n), note="; ".join(notes))
        if spec is None or not len(spec) or not n:
            sp["msg"] = "%s: no %s ion scans at %.2f to %.2f min" % (name, word, a, b)
            return sp
        sp["spec"] = np.asarray(spec, float)
        return sp

    def _mz_label(self):
        pk = self.mzpick
        en = pk["entry"]
        return self.label_of(en) if any(e is en for e in self.entries) else pk["label"]

    def _mz_texts(self, pol):
        """(head, time, details) of a spectrum tile: 'ESI+ · run_B', '9.35 to 9.43 min (peak 9.39 min, MS delay
        0.05 min)', '28 scans averaged   ·   ...'."""
        pk = self.mzpick
        sp = pk["spectra"][pol]
        head = "%s · %s" % (self._pol_name(pol), self._mz_label())
        dl = ", MS delay %.2f min" % pk["delay"] if pk["pda"] else ""
        a, b = sp["range"]
        kind = pk["kind"]
        own = ", background %.2f to %.2f min" % tuple(sp["bg"]) if sp["bg"] and sp.get("bg_manual") else ""
        if kind == "range":  # a time range dragged across the trace (and the background range of the run)
            when = ("averaged %.2f to %.2f min (range)" % (a, b) if b - a >= 1e-9 else "scan at %.2f min" % a) + own + dl
        elif kind == "scan":  # Alt + click: one scan
            when = "scan at %.2f min" % a + own + dl
        elif kind == "none":
            when = "%.2f to %.2f min (no peak at %.2f min%s)" % (a, b, pk["ms_apex"], dl)
        elif b - a < 1e-9:
            when = "scan at %.2f min (peak %.2f min%s)" % (a, pk["ms_apex"], dl)
        else:
            when = "%.2f to %.2f min (peak %.2f min%s)" % (a, b, pk["ms_apex"], dl)
        sub = ["%d scans averaged" % sp["n"] if sp["n"] > 1 else "1 scan"]
        if sp["bg"] and kind not in ("range", "scan"):
            sub.append(("background %.2f to %.2f min subtracted" if sp.get("bg_manual") else
                        "spectrum of %.2f to %.2f min subtracted") % tuple(sp["bg"]))
        if sp["note"]:
            sub.append(sp["note"])
        s_ = sp["spec"]
        if s_ is not None and len(s_):
            ib = int(np.argmax(s_[:, 1]))
            sub.append("base peak m/z %.1f (%s)" % (s_[ib, 0], U._g(s_[ib, 1])))
        return head, when, "   ·   ".join(sub)

    def _mz_head(self):
        """The line above the spectrum tiles: the suggested neutral mass, worded as a suggestion."""
        pk = self.mzpick
        if pk is None:
            bold, rest = "", "Click a peak of a trace for its mass spectra"
        else:
            sug = K.mass_suggestion(pk.get("masses"), info=pk.get("mass_info"))
            have = [q for q in ("+", "-") if self._spec(q) is not None]
            if sug:
                best = pk["masses"][0]
                bold = sug
                pols = set(i.get("polarity", "+") for i in (best.get("ions") or []) if isinstance(i, dict))
                rest = ("suggested from both spectra" if best.get("both") or pols == {"+", "-"} else
                        "suggested from ESI+ only" if pols == {"+"} else
                        "suggested from ESI− only" if pols == {"-"} else "suggested")
            elif have:
                bold, rest = "", "No neutral mass suggested from the ions of this %s" % (
                    "time range" if pk["kind"] == "range" else "scan" if pk["kind"] == "scan" else "peak")
            else:
                bold, rest = "", "No mass spectra for this %s of %s" % (
                    "time range" if pk["kind"] == "range" else "time" if pk["kind"] == "scan" else "peak",
                    self._mz_label())
        try:
            self.mzrow.set_head(bold, rest)
        except RuntimeError:
            pass

    def _mz_retitle(self):
        """Titles of the spectrum tiles again (the name of the run follows the list: labels typed, files
        opened or closed)."""
        if self.mzpick is None:
            return
        for p in ("+", "-"):
            c = self.mzcards[p]
            if getattr(c, "_mz_parts", None):
                c._mz_parts = self._mz_texts(p)
                self._fit_mz_title(p)
            else:
                c.set_title(self._pol_name(p) + " · " + self._mz_label(), "")

    def _fit_mz_title(self, pol):
        """Title of a spectrum tile: run, polarity and times when they fit, else the times start the
        subtitle (which the tile shortens; the tooltip of the tile has all)."""
        try:
            card = self.mzcards[pol]
            parts = getattr(card, "_mz_parts", None)
            if not parts:
                return
            head, when, sub = parts
            full = head + " · " + when
            w = card.GetClientSize()[0]
            fits = True
            if w >= 50:
                dc = wx.ClientDC(card)
                dc.SetFont(ui_font(10, 700))
                fits = dc.GetTextExtent(full)[0] <= w - card.FromDIP(14 + 12 + 170 + 12)
            title, s2 = (full, sub) if fits else (head, when + ("   ·   " + sub if sub else ""))
            if (title, s2) != (card.title, card.subtitle):
                card.set_title(title, s2)
        except (RuntimeError, AttributeError, KeyError):
            pass  # closed meanwhile

    def plot_mz(self, pol, keep_view=False):
        """A spectrum tile of the m/z tool, drawn as the Mass spectrometry view draws its spectra (unilcms helpers:
        profile or sticks, the strongest ions labelled with their m/z, labels pinned, the measurement)."""
        card = self.mzcards[pol]
        col = self.mzcols[pol]
        ax = card.ax
        view = ax.get_xlim() if (keep_view and card.full) else None
        card.reset()
        card._mz_parts = None
        card._mz_auto = []
        U.style_axes(ax, "m/z", "Intensity")
        pk = self.mzpick
        sp = pk["spectra"].get(pol) if pk else None
        spec = sp.get("spec") if sp else None
        col.update(spec=spec if (spec is not None and len(spec)) else None, e=sp.get("e") if sp else None,
                   range=sp.get("range") if sp else None)
        if col["spec"] is None:
            col.update(measure=None, m1=None)
            if pk is None:  # (how to pick one: the line above the tiles)
                msg = "No peak picked yet"
            else:
                msg = (sp or {}).get("msg") or "No spectrum"
            ax.text(0.5, 0.5, msg, transform=ax.transAxes, ha="center", va="center", fontsize=9, color=U.MUTED,
                    wrap=True)
            card.full = None
            card.set_title(self._pol_name(pol) + (" · " + self._mz_label() if pk else ""), "")
            card.draw()
            self._mz_head()
            return
        d = pk["doc"]
        mz = spec[:, 0]
        ev = {}
        try:
            ev = d.attrs.get("ms_data").events[sp["e"]]
        except Exception:
            ev = {}
        lo = ev.get("mz_low") or float(mz.min())
        hi = ev.get("mz_high") or float(mz.max())
        n = int(self.s.get("spec_labels", 8))
        try:  # a narrow tile: fewer labels (they ran into each other)
            n = min(n, max(2, int(ax.get_window_extent().width * 72.0 / card.fig.dpi / 40.0)))  # 40 pt each
        except Exception:
            pass
        auto = U.draw_spectrum(card, spec, (lo, hi) if hi > lo else (lo - 1.0, hi + 1.0), self._sticks(), view, n,
                               self.MZ_FMT)
        U.draw_spec_marks(ax, col, auto, self.PIN_FMT, lambda m: U.spec_peak_height(spec, m, self.PICK_MIN))
        card._mz_auto = auto
        card.place()  # (the frame as drawn: the labels are checked for it)
        self._unclash_labels(pol)
        col["desc"] = " ".join(self._mz_texts(pol)[1:2])
        card._mz_parts = self._mz_texts(pol)
        self._fit_mz_title(pol)
        card.draw()
        self._mz_head()

    def _ms_to_plot(self, o, pk_delay=None):
        """MS time -> processed trace time of the drawn trace o: a PDA trace minus the delay of the file (also for
        a time picked on an MS trace, whose own delay is 0), plus its alignment and shift."""
        sh = float(o.get("align", 0.0) or 0.0) + float(o["entry"].get("shift", 0.0) or 0.0)
        back = (self._delay(o["entry"]["doc"]) if pk_delay is None else pk_delay) \
            if self.s["signal"] in ("pda", "max") else 0.0
        return sh - back

    def _shade(self, ax, o, a, b, colour, alpha, hatch=None, under=None):
        """A time range a to b (processed trace time) shaded on the drawn trace o only: a band as tall as the trace
        (as the Mass spectrometry view shades its ranges on the chromatogram) and, if under, the area under the
        trace there."""
        pr = o["proc"]
        t, y = pr["t"], pr["y"]
        ok = np.isfinite(t) & np.isfinite(y)
        t, y = t[ok], y[ok]
        if len(t) < 2 or not b >= a:
            return None
        a, b = max(a, float(t[0])), min(b, float(t[-1]))
        if b < a:
            return None
        top = float(np.max(y)) if len(y) else 0.0
        base = min(0.0, float(np.min(y)))
        if not top > base:
            top = base + 1.0
        z = 3.0 + 0.5  # above the traces (stacked: over a neighbour's fill)
        kw = dict(facecolor=colour, alpha=alpha, lw=0, zorder=z - 0.7)
        if hatch:
            kw.update(hatch=hatch, edgecolor="#FFFFFF")
        if b - a < 1e-9:  # one scan: a thin line
            ln = ax.vlines([a + o["dx"]], o["dy"] + base, o["dy"] + top, color=colour, lw=1.2, alpha=min(1.0, 3 * alpha),
                           zorder=z - 0.7)
            return ln
        band = ax.fill_between([a + o["dx"], b + o["dx"]], o["dy"] + base, o["dy"] + top, **kw)
        if under:
            m = (t >= a) & (t <= b)
            tt, yy = t[m], y[m]
            # the ends of the range on the trace itself
            tt = np.concatenate([[a], tt, [b]])
            yy = np.concatenate([[np.interp(a, t, y)], yy, [np.interp(b, t, y)]])
            ax.fill_between(tt + o["dx"], o["dy"], yy + o["dy"], color=colour, alpha=under, lw=0, zorder=z - 0.6)
        return band

    def _unclash_labels(self, pol):
        """The automatic m/z labels of a tile (tallest first) that would run into a taller one, a pinned label or the
        measurement are hidden, for the present size of the tile (a narrow tile); also for copied and saved images
        and the report (PlotCard.relabel_fn), whose size differs."""
        card = self.mzcards[pol]
        try:
            texts = card.ax.texts
            auto = [a for a in (getattr(card, "_mz_auto", None) or []) if a in texts]
            for a in auto:
                a.set_visible(True)
            r = card.fig.canvas.get_renderer()
            mine = set(id(a) for a in auto)
            kept = [t.get_window_extent(r) for t in texts if id(t) not in mine and t.get_visible() and t.get_text()]
            pad = 1.5 / (72.0 / card.fig.dpi)  # 1.5 pt between two labels
            for a in auto:
                bb = a.get_window_extent(r).padded(pad)
                if any(bb.overlaps(k) for k in kept):
                    a.set_visible(False)
                else:
                    kept.append(bb)
        except Exception as ex:
            U._log("compare m/z: labels: %s" % ex)

    def _relabel_tile(self, pol):
        """A spectrum tile resized: its labels checked again for its size."""
        try:
            if self and getattr(self.mzcards[pol], "_mz_auto", None):
                self._unclash_labels(pol)
                self.mzcards[pol].canvas.draw_idle()
        except RuntimeError:
            pass  # closed meanwhile

    def _draw_mz_mark(self, ax, out):
        """The pick on its trace (part of the plot, so also in copied images and the report): a peak clicked as a
        light shading of its range under the trace and a ring at its top; a time range dragged as a band on the
        trace (the area under it shaded more); a single scan (Alt + click) as a ring. The background range of a run
        (Shift + drag) as a grey hatched band on its trace."""
        pk = self.mzpick
        self._mz_ring = None  # (the ring drawn, for the report: mz_report)
        self._mz_markx = None  # (where the pick is marked, plot time)
        self._mz_bgs = []  # (the runs with a background range drawn)
        if not self.mz_on:
            return
        for o in out:
            bgr = o["entry"].get("mzbg")
            if not bgr:
                continue
            try:
                d = self._ms_to_plot(o)
                if self._shade(ax, o, bgr[0] + d, bgr[1] + d, U.BG_COL, 0.3, hatch="///") is not None:
                    self._mz_bgs.append(o["label"])
            except Exception as ex:
                U._log("compare m/z: background mark: %s" % ex)
        if not pk:
            return
        o = next((o for o in out if o["entry"]["doc"] is pk["doc"]), None)
        if o is None:
            return
        try:
            pr = o["proc"]
            t, y = pr["t"], pr["y"]
            if len(t) < 2:
                return
            d = self._ms_to_plot(o, pk.get("file_delay", pk["delay"]))
            a, b, ta = pk["ms_t0"] + d, pk["ms_t1"] + d, pk["ms_apex"] + d
            col = o["colour"]
            z = 3.0 + 0.5  # above the traces (stacked: over a neighbour's fill)
            if pk["kind"] == "range":
                if self._shade(ax, o, a, b, col, 0.13, under=0.28) is not None:
                    self._mz_markx = 0.5 * (max(a, float(t[0])) + min(b, float(t[-1]))) + o["dx"]
                return
            if b > a:
                m = (t >= a) & (t <= b)
                if np.count_nonzero(m) >= 2:
                    ax.fill_between(t[m] + o["dx"], o["dy"], y[m] + o["dy"], color=col, alpha=0.28, lw=0, zorder=z - 0.6)
            if t[0] <= ta <= t[-1]:
                ya = float(np.interp(ta, t, y)) + o["dy"]
                mk, = ax.plot([ta + o["dx"]], [ya], ls="none", marker="o", ms=5.5, mfc="white", mec=_ink(col),
                              mew=1.3, zorder=z + 2)
                mk._no_scale = True
                self._mz_ring = mk
                self._mz_markx = ta + o["dx"]
        except Exception as ex:
            U._log("compare m/z: mark: %s" % ex)

    def on_spec_opts(self, e=None):
        if self._busy:
            return
        old = (self.s.get("spec_avg"), bool(self.s.get("spec_bg")))
        self.s["spec_avg"] = _key(SPEC_AVG, self.spec_avg.GetSelection())
        self.s["spec_bg"] = bool(self.spec_bg.GetValue())
        self.s["spec_labels"] = int(self.spec_nlab.GetValue())
        self._save()
        if self.mzpick is None:
            return
        if (self.s["spec_avg"], self.s["spec_bg"]) != old:
            self._mz_compute()
        for p in ("+", "-"):
            self.plot_mz(p, keep_view=True)

    def _spec(self, pol):
        pk = self.mzpick
        sp = pk["spectra"].get(pol) if pk else None
        return sp if (sp and sp.get("spec") is not None and len(sp["spec"])) else None

    def _nearest_ion(self, pol, x):
        """m/z of the strongest point near x in the spectrum (as the Mass spectrometry view picks a peak)."""
        sp = self._spec(pol)
        if sp is None or x is None:
            return None
        m = U.spec_nearest_peak(sp["spec"], self.mzcards[pol].ax.get_xlim(), x, self.PICK_MIN)
        if m is None or not U.spec_peak_height(sp["spec"], m, 1e-9) > 0:  # nothing there to label
            return None
        return m

    def _sign(self, pol):
        return -1 if pol == "-" else 1

    def spec_tool(self, pol, tool, x, y):
        """The spectrum tools of the Mass spectrometry view on a tile: Mass chrom. (that ion in every run), Label (pin
        or unpin its label), Measure (two peaks: distance, isotope spacing or charges and mass, for this polarity)."""
        col = self.mzcols[pol]
        if self._spec(pol) is None:
            return
        pk = self._nearest_ion(pol, x)
        if pk is None:
            self.status("Click closer to a peak")
            return
        if tool == "xic":
            self.track_mz(pol, pk)
        elif U.spec_click(col, tool, pk, self.PIN_FMT, lambda m1, m2: U.measure_text(m1, m2, self._sign(pol)),
                          self.status):
            self.plot_mz(pol, keep_view=True)

    def _sticks(self):
        """Display of the spectra: chosen here (right click a tile > Display), else as the Mass spectrometry view of
        the file picked."""
        if self._spec_style is not None:
            return self._spec_style == 1
        try:
            return self.mzpick["doc"].attrs.get("ms").spec_style.GetSelection() == 1
        except Exception:
            return False

    def _set_style(self, i):
        self._spec_style = 1 if i == 1 else 0
        for p in ("+", "-"):
            self.plot_mz(p, keep_view=True)

    def spec_menu(self, pol, x, y):
        """Right click a spectrum tile: the menu of the spectra of the Mass spectrometry view, for the Compare view
        (the mass chromatogram of an ion is followed in every run)."""
        sp = self._spec(pol)
        if sp is None:
            return []
        col = self.mzcols[pol]
        items = []
        m = self._nearest_ion(pol, x)
        if m is not None:
            items.append((("Mass chromatogram of m/z " + self.PIN_FMT + " in every run (%s)") % (m, self._pol_name(pol)),
                          lambda: self.track_mz(pol, m)))
            items.append(U.spec_pin_item(col, m, self.PIN_FMT, lambda v: self.spec_tool(pol, "label", v, None)))
            items.append((None, None))
        items += U.spec_clear_items(col, lambda: self.plot_mz(pol, keep_view=True))
        items.append(U.spec_display_item(self._sticks(), self._set_style))
        items.append(("Show in the Mass spectrometry view", lambda: self.show_in_ms(pol)))
        items += [("Export spectrum\u2026", lambda: self.export_spectrum(pol)),
                  ("Copy spectrum data", lambda: U.copy_spectrum(self, sp["spec"])),
                  ("Open this spectrum in the Deconvolute window", lambda: self.open_full(pol))]
        items += [(None, None), ("Full view", self.mzcards[pol].reset_view)]
        return items + self.view_menu()

    def open_full(self, pol):
        """The spectrum in the full Deconvolute window (separate), written next to its file as the Mass
        spectrometry view does (negative ions as JCAMP-DX, so that it starts in negative mode)."""
        sp = self._spec(pol)
        pk = self.mzpick
        if sp is None or pk is None:
            return None
        return U.open_spectrum_in_deconvolute(self, pk["doc"], sp["spec"], self._spec_name(pol), self._sign(pol),
                                              sp["range"])

    def track_mz(self, pol, m):
        """The compared signal becomes the mass chromatogram of this ion (window 0.5) in every run."""
        self._busy = True
        try:
            self.signal.SetSelection(_index(K.SIGNALS, "xic"))
            self.pol.SetSelection(0 if pol == "+" else 1)
            self.mz.SetValue("%g" % round(float(m), 2))
            self.mz_win.SetValue("0.5")
            self.mz_ppm.SetValue(False)
        finally:
            self._busy = False
        self.on_change()
        self.status("Mass chromatogram of m/z %.2f (%s) in every run" % (m, self._pol_name(pol)))

    def show_in_ms(self, pol):
        """The Mass spectrometry view of the picked file, with the time range of this spectrum averaged."""
        sp = self._spec(pol)
        pk = self.mzpick
        if sp is None or pk is None:
            return
        d = pk["doc"]
        fr = self.frame
        ms_tab = d.attrs.get("ms")
        if d not in getattr(fr, "docs", []) or d.attrs.get("ms_data") is None or ms_tab is None:
            self.status("%s is not open any more" % pk["label"])
            return
        a, b = sp["range"]
        try:
            fr.activate(d)
            d.show_page(ms_tab)  # leaves the Compare view
            if b - a < 1e-9:
                if sp["bg"]:
                    ms_tab.bg = sp["bg"]  # (its range shown on the chromatograms)
                    ms_tab.use_bg.SetValue(True)
                ms_tab.chrom_pick(a, None, None)
                if sp["bg"]:
                    # one scan minus the spectrum before the peak: the MS view shows a single scan without a
                    # background, so its spectra become those of the tiles (the same scans and background)
                    for col in ms_tab.cols:
                        s2 = self._spec(ms_tab.pol(col["e"]) or "+")
                        if s2 is None or s2["e"] != col["e"] or not s2["bg"]:
                            continue
                        col["spec"] = s2["spec"].copy()
                        col["desc"] = "scan at %.3f min, spectrum of %.2f to %.2f min subtracted" % (
                            (s2["range"][0],) + tuple(s2["bg"]))
                        ms_tab.plot_spec(col, keep_view=True)
                    ms_tab.spectra_changed()
                return
            ms_tab.avg = (a, b)
            ms_tab.pick_t = None
            if sp["bg"]:
                ms_tab.bg = sp["bg"]
                ms_tab.use_bg.SetValue(True)
            elif ms_tab.use_bg.GetValue() and ms_tab.bg:
                ms_tab.use_bg.SetValue(False)  # the same spectrum as in the Compare view
            ms_tab.fill_range_fields()
            ms_tab.plot_chroms(keep_view=True)
            ms_tab.compute_average()
            try:
                d.link_time("ms", rng=(a, b))  # and the PDA view at that time
            except Exception:
                pass
        except Exception as ex:
            U._report_error(self, "The Mass spectrometry view could not be shown", ex)

    def _spec_name(self, pol):
        pk = self.mzpick
        sp = self._spec(pol)
        if pk is None or sp is None:
            return None
        base = os.path.splitext(os.path.basename((pk["doc"].attrs.get("path") or "run").rstrip("\\/")))[0]
        a, b = sp["range"]
        return U.safe_file_name("%s_%s_%.2f-%.2fmin" % (base, "pos" if pol == "+" else "neg", a, b))

    def export_spectrum(self, pol, path=None):
        sp = self._spec(pol)
        pk = self.mzpick
        if sp is None:
            self.status("No spectrum to export")
            return None
        if not path:
            try:
                path = pk["doc"].ask_save("Export mass spectrum %s" % self._pol_name(pol), self._spec_name(pol) + ".txt",
                                          U.SPEC_WILDCARD + "|CSV with the run and times (*.csv)|*.csv")
            except Exception:
                path = self._ask("Export the spectrum", self._spec_name(pol) + ".csv")
        if not path:
            return None
        a, b = sp["range"]
        if not path.lower().endswith(".csv"):  # text or JCAMP-DX, as the Mass spectrometry view exports
            try:
                U.write_spectrum_file(path, sp["spec"], self._spec_name(pol), pol, (a, b))
                self.status("Spectrum saved: %s" % path)
                return path
            except Exception as ex:
                U._report_error(self, U._write_failed(path), ex)
                return None
        meta = [("File", os.path.basename((pk["doc"].attrs.get("path") or "").rstrip("\\/"))),
                ("Run", self._mz_label()), ("Polarity", self._pol_name(pol)),
                ("Time range (min)", "%.4f to %.4f" % (a, b)), ("Scans averaged", sp["n"])]
        if pk["kind"] in ("peak", "none"):
            meta.append(("Peak top (min)", "%.4f" % pk["ms_apex"]))
        if pk["pda"]:
            meta.append(("MS detector delay (min)", "%g" % pk["delay"]))
        if sp["bg"]:
            meta.append(("Subtracted", "%.4f to %.4f min%s" % (tuple(sp["bg"]) + (
                " (background range)" if sp.get("bg_manual") else "",))))
        meta.append(("Bin width (m/z)", "%g" % pk.get("binw", 0.05)))
        try:
            K.write_table_csv(path, ["m/z", "Intensity"], [["%.4f" % m_, "%.6g" % i_] for m_, i_ in sp["spec"]], meta)
            self.status("Spectrum saved: %s" % path)
            return path
        except Exception as ex:
            U._report_error(self, "The spectrum could not be saved", ex)

    def mz_report(self):
        """The picked peak and its tiles for the Comparison report, or None (m/z off, nothing picked)."""
        pk = self.mzpick if self.mz_on else None
        if not pk:
            return None
        cards = [(p, self.mzcards[p]) for p in ("+", "-") if self._spec(p) is not None]
        if not cards:
            return None
        # the pick marked in the plot as shown (not when its run is unticked or the blank, or the peak lies outside the
        # time window or the view)
        mx = getattr(self, "_mz_markx", None)
        x0, x1 = self.card.ax.get_xlim()
        marked = mx is not None and x0 <= mx <= x1
        return {"cards": cards, "pick": pk, "label": self._mz_label(), "avg": self.s.get("spec_avg", "peak"),
                "suggestion": K.mass_suggestion(pk.get("masses"), info=pk.get("mass_info")), "marked": marked,
                "kind": pk.get("kind", "peak"), "bgs": list(getattr(self, "_mz_bgs", []) or [])}

    # ------------------------------------------------------------------ export
    def image_name(self):
        s = self.s
        tag = ("%gnm" % s["wl"]) if (s["signal"] == "pda" and s["wl"] is not None) else s["signal"]
        return "comparison_%s" % tag

    def _ask_file(self, title, name, wildcard, extension):
        d = self.entries[0]["doc"] if self.entries else None
        folder = os.path.dirname((d.attrs.get("path") or "").rstrip("\\/")) if d is not None else ""
        dlg = wx.FileDialog(self, title, defaultDir=folder, defaultFile=name, wildcard=wildcard,
                            style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        try:
            if dlg.ShowModal() != wx.ID_OK:
                return None
            p = dlg.GetPath()
            return p if p.lower().endswith(extension.lower()) else p + extension
        finally:
            dlg.Destroy()

    def _ask(self, title, name):
        return self._ask_file(title, name, "CSV (*.csv)|*.csv", ".csv")

    def on_export_traces(self, e=None, path=None):
        if not self.drawn:
            self.status("No trace to export")
            return None
        path = path or self._ask("Export the traces", self.image_name() + "_traces.csv")
        if not path:
            return None
        try:
            K.write_traces_csv(path, [{"label": o["label"], "proc": o["proc"], "units": o["units"],
                                       "scaled": o["scaled"]} for o in self.drawn])
            self.status("Traces saved: %s" % path)
            return path
        except Exception as ex:
            U._report_error(self, "The traces could not be saved", ex)

    def on_export_table(self, e=None, path=None):
        if self.area_range is not None and not self._prepare_area_action():
            return None
        tb = getattr(self, "_table", None)
        if not self.drawn or not tb:
            self.status("No table to export")
            return None
        path = path or self._ask("Export the table", self.image_name() + "_peaks.csv")
        if not path:
            return None
        try:
            K.write_table_csv(path, tb[0], tb[1], [("Signal", self.card.subtitle)])
            self.status("Table saved: %s" % path)
            return path
        except Exception as ex:
            U._report_error(self, "The table could not be saved", ex)

    def on_export_areas_xlsx(self, e=None, path=None):
        if not self._prepare_area_action():
            return None
        path = path or self._ask_file("Export the selected peak areas", self.image_name() + "_areas.xlsx",
                                     "Excel workbook (*.xlsx)|*.xlsx", ".xlsx")
        if not path:
            return None
        xlabel = self.area_xlabel.GetValue().strip() or "X"
        units = sorted(set(r.get("units", "") for r in self.area_results))
        unit = units[0] if len(units) == 1 else ""
        header = ["Trace", "Region start (min)", "Region end (min)", "Region peak (min)",
                  ("Height (%s)" % unit) if unit else "Height",
                  ("Actual area (%s·s)" % unit) if unit else "Actual area", "Relative area (%)",
                  "100% reference", xlabel]
        rows = [[r["label"], r["t0"], r["t1"],
                 r["rt"] - self._shift_of(self.drawn[r["index"]]) if r["rt"] is not None else "", r["height"],
                 r["area"] if r["area"] is not None else "no data", r["relative"], "yes" if r["reference"] else "",
                 r["x"]]
                for r in self.area_results]
        meta = [("Signal", self.card.subtitle), ("Displayed region start (min)", self.area_range[0]),
                ("Displayed region end (min)", self.area_range[1]),
                ("Region peak", "retention time in the run, before its shift"),
                ("Integration", "straight baseline between the signal at both ends; positive area only"),
                ("Display normalisation", "not applied to actual areas")]
        if self._area_x_error:
            meta.append(("X not used", self._area_x_error))
        try:
            K.write_area_xlsx(path, header, rows, meta)
            self.status("Peak areas saved to Excel: %s" % path)
            return path
        except Exception as ex:
            U._report_error(self, "The peak areas could not be saved", ex)

"""
LCMS Analysis (MSpektra): LC-MS and HPLC workspace for Shimadzu LabSolutions
.lcd files and the files of other vendors (lcms_sources.py). The shared plotting, tool bar, integration and deconvolution classes
here are also used by HRMS Analysis (hrms.py).

Part of the portable add-ons; the deconvolution engine itself is unchanged.

Mass spectrometry tab
  tiles of each scan event one below the other (chromatograms first, then
  the spectra): TIC or base
  peak chromatogram, mass chromatograms (right click), spectrum at a time or
  averaged over a range with background subtraction, peak integration,
  export (.txt, .jdx, .csv) and deconvolution in the same window.
PDA tab
  wavelength-time map with a time slider (top edge) and a wavelength slider
  (right edge); the chromatogram and the UV spectrum follow the sliders, and
  the mass spectra of the MS data follow the time chosen (MS detector delay).
  More wavelengths, max plot, UV averages with background subtraction,
  overlay of the MS trace with an adjustable delay, peak integration, export.

Data readers: lcms_data.py (MS, uses OpenSZRaw) and lcms_pda.py (PDA) for .lcd
files; lcms_sources.py for Agilent, Waters, Thermo, mzML, mzXML and ANDI files
(data_formats.py tells the formats apart).
"""
import os
import sys
import time
import threading
import datetime
import types
import subprocess

import numpy as np
import wx

import unidec_theme as T
from unidec_theme import C, ui_font
import lcms_data
import lcms_pda
import lcms_integrate as LI
import data_formats

INK = "#000000"  # axes, ticks and labels
TRACE = "#1D4ED8"  # the main trace (spectra, chromatograms, UV)
NEG_COL = "#B3261E"  # negative ions, when both polarities are shown (not the red of saturated marks)
LW_TRACE, LW_SPEC, LW_STICK = 0.7, 0.6, 0.7  # line widths (pt): chromatograms and UV, profile spectra, centroids
MUTED = "#6B7480"
PAL = ["#2A62C4", "#E0602F", "#0BA064", "#8E5BD0", "#C99A1E", "#D2456F", "#5B6B7C", "#9C6B43"]
AVG_COL = "#2F6FCB"
BG_COL = "#8A94A3"
SEARCH_COLS = ["#0F766E", "#B45309", "#7C3AED", "#BE185D", "#4D7C0F", "#0369A1"]  # the masses of a mass search
NOTE_COL = "#B45309"  # the note of a mass search (another mass fits the same m/z)
ABSENT_COL = "#C4C9D0"  # charge states of a found mass that are not present
TITLE = "LCMS Analysis"
APP = "MSpektra"
OUT_SUFFIX = "_analysis"  # per file output folder: <data file name>_analysis


def _log(*a):
    try:
        print("[%s]" % TITLE, *a)
    except Exception:
        pass


def _num(ctrl, default=None):
    try:
        return float(ctrl.GetValue().strip().replace(",", "."))
    except Exception:
        return default


def _num_s(text, default=None):
    try:
        return float(str(text).replace(",", "."))
    except (TypeError, ValueError):
        return default


def _numlist(text):
    out = []
    for part in text.replace(";", ",").replace(" ", ",").split(","):
        part = part.strip()
        if part:
            try:
                out.append(float(part))
            except ValueError:
                pass
    return out


_XIC_RX = None


def parse_xics(text, default_w, default_ppm):
    """'803.2, 1204.8+-0.3, 889.09+-10ppm' -> [{'mz', 'w', 'ppm'}]."""
    import re
    global _XIC_RX
    if _XIC_RX is None:
        _XIC_RX = re.compile(r"^([0-9]*\.?[0-9]+)(?:(?:\u00b1|\+/-|\+-|/)([0-9]*\.?[0-9]+)(ppm|mz|m/z)?)?$", re.I)
    out = []
    for part in re.split(r"[,;]", text or ""):
        part = part.strip().replace(" ", "")
        if not part:
            continue
        m = _XIC_RX.match(part)
        if not m:
            continue
        w = float(m.group(2)) if m.group(2) else default_w
        ppm = (m.group(3) or "").lower() == "ppm" if m.group(2) else default_ppm
        out.append({"mz": float(m.group(1)), "w": w, "ppm": ppm})
    return out


def xic_tol(x):
    return x["mz"] * x["w"] * 1e-6 if x["ppm"] else x["w"]


def xic_key(x):
    return (round(x["mz"], 5), round(x["w"], 6), bool(x["ppm"]))


def xic_window_text(x):
    return ("\u00b1%g ppm" % x["w"]) if x["ppm"] else ("\u00b1%g" % x["w"])


class XicDialog(wx.Dialog):
    """New or edited mass chromatogram: m/z (one or a list), window and unit,
    scan event."""

    def __init__(self, parent, title, mz_text, w, ppm, events, event, allow_list=True):
        wx.Dialog.__init__(self, wx.GetTopLevelParent(parent), title=title)
        self.SetBackgroundColour(wx.Colour(C["panel"]))
        g = wx.FlexGridSizer(0, 2, self.FromDIP(8), self.FromDIP(10))
        self.mz = wx.TextCtrl(self, value=mz_text, size=wx.Size(self.FromDIP(200), -1))
        self.mz.SetToolTip("m/z value" + ("s, separated by commas" if allow_list else ""))
        self.w = wx.TextCtrl(self, value="%g" % w, size=wx.Size(self.FromDIP(80), -1))
        self.unit = wx.Choice(self, choices=["m/z", "ppm"])
        self.unit.SetSelection(1 if ppm else 0)
        wbox = wx.BoxSizer(wx.HORIZONTAL)
        wbox.Add(wx.StaticText(self, label="\u00b1"), 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, self.FromDIP(4))
        wbox.Add(self.w, 0, wx.ALIGN_CENTER_VERTICAL)
        wbox.Add(self.unit, 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT, self.FromDIP(6))
        self.ev = wx.Choice(self, choices=list(events) + (["All scan events"] if len(events) > 1 and allow_list else []))
        self.ev.SetSelection(max(0, min(event, len(events) - 1)))
        if len(events) < 2 or not allow_list:
            self.ev.Enable(False)
        for lab, ctrl in (("m/z", self.mz), ("Window", wbox), ("Scans", self.ev)):
            g.Add(wx.StaticText(self, label=lab), 0, wx.ALIGN_CENTER_VERTICAL)
            g.Add(ctrl, 0, wx.ALIGN_CENTER_VERTICAL)
        vs = wx.BoxSizer(wx.VERTICAL)
        pad = self.FromDIP(14)
        vs.Add(g, 0, wx.ALL, pad)
        vs.Add(self.CreateButtonSizer(wx.OK | wx.CANCEL), 0, wx.EXPAND | wx.ALL, pad)
        self.SetSizerAndFit(vs)
        self.CentreOnParent()
        self.mz.SetFocus()

    def values(self):
        mzs = [v for v in _numlist(self.mz.GetValue()) if v > 0]
        w = _num(self.w, None)
        if w is None or w <= 0:
            w = 0.5 if self.unit.GetSelection() == 0 else 10.0
        ev = self.ev.GetSelection()
        if self.ev.GetCount() and ev == self.ev.GetCount() - 1 and self.ev.GetString(ev) == "All scan events":
            ev = -1
        return {"mz": mzs, "w": w, "ppm": self.unit.GetSelection() == 1, "event": ev}


def search_text(s, n=3):
    """The masses of a search in a few words (undo steps, report): names where given."""
    rows = search_rows(s)
    return ", ".join(nm or _mass_text(M) for M, w, nm in rows[:n]) + (" \u2026" if len(rows) > n else "")


def search_rows(s):
    """The masses of a mass search: [(mass, window, name)] (searches have "masses", "ws", "names")."""
    ms = list(s.get("masses") or [])
    ws = list(s.get("ws") or [s.get("w", 5.0)] * len(ms))
    nm = list(s.get("names") or [""] * len(ms))
    return [(float(m), float(ws[k] if k < len(ws) else 5.0), str(nm[k] if k < len(nm) else "")) for k, m in enumerate(ms)]


def _mass_text(m):
    """18366, 18364.5, 12243.25 (no trailing zeros)."""
    return ("%.2f" % float(m)).rstrip("0").rstrip(".")


class FindMassDialog(wx.Dialog):
    """Find mass: a table of masses (one row each: neutral mass, window in Da, optional name; an empty row at the
    end for a new mass, x removes a row), then the charge carrier, the lowest m/z of the charge states and the
    scan event, shared by the rows. pols: the polarity of each event (negative ions: the carrier is a lost
    proton)."""

    LOST = "Lost H+"
    MAX_ROWS = 12

    def __init__(self, parent, title, v, events, pols, event):
        import ms_find_mass as F
        wx.Dialog.__init__(self, wx.GetTopLevelParent(parent), title=title)
        self.SetBackgroundColour(wx.Colour(C["panel"]))
        self.pols = list(pols)
        self.names = [n for n, m in F.CARRIERS]
        self.rows = []  # [(mass, window, name, remove button)]
        self.table = wx.FlexGridSizer(0, 4, self.FromDIP(4), self.FromDIP(6))
        for lab in ("Mass (Da)", "\u00b1 (Da)", "Name", ""):
            st = wx.StaticText(self, label=lab)
            st.SetForegroundColour(wx.Colour(C["muted"]))
            self.table.Add(st, 0, wx.ALIGN_BOTTOM)
        for r in list(v.get("rows") or [])[:self.MAX_ROWS]:
            self._add_row(*r)
        self._add_row()
        self.ion = wx.Choice(self, choices=self.names + [self.LOST])
        self.ion.SetMinSize(wx.Size(self.FromDIP(90), -1))
        self._ion = self.names.index(v["carrier"]) if v.get("carrier") in self.names else 0
        self.lo = wx.TextCtrl(self, value="%g" % v.get("mzmin", F.MZMIN), size=wx.Size(self.FromDIP(80), -1))
        lbox = wx.BoxSizer(wx.HORIZONTAL)
        lbox.Add(self.lo, 0, wx.ALIGN_CENTER_VERTICAL)
        lbox.Add(wx.StaticText(self, label="m/z and up"), 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT, self.FromDIP(6))
        self.ev = wx.Choice(self, choices=list(events) + (["All scan events"] if len(events) > 1 else []))
        self.ev.SetSelection(len(events) if (event is not None and event < 0 and len(events) > 1) else
                             max(0, min(event or 0, len(events) - 1)))
        if len(events) < 2:
            self.ev.Enable(False)
        g = wx.FlexGridSizer(0, 2, self.FromDIP(8), self.FromDIP(10))
        for lab, ctrl in (("Charge carrier", self.ion), ("Charge states at", lbox), ("Scans", self.ev)):
            g.Add(wx.StaticText(self, label=lab), 0, wx.ALIGN_CENTER_VERTICAL)
            g.Add(ctrl, 0, wx.ALIGN_CENTER_VERTICAL)
        vs = wx.BoxSizer(wx.VERTICAL)
        pad = self.FromDIP(14)
        vs.Add(self.table, 0, wx.LEFT | wx.RIGHT | wx.TOP, pad)
        vs.Add(g, 0, wx.LEFT | wx.RIGHT | wx.TOP, pad)
        vs.Add(self.CreateButtonSizer(wx.OK | wx.CANCEL), 0, wx.EXPAND | wx.ALL, pad)
        self.SetSizerAndFit(vs)
        self.ion.Bind(wx.EVT_CHOICE, self._ion_chosen)
        self.ev.Bind(wx.EVT_CHOICE, lambda e: self._carrier_shown())
        self._carrier_shown()
        self.CentreOnParent()
        self.rows[0][0].SetFocus()

    def _add_row(self, mass="", w="5", name=""):
        m = wx.TextCtrl(self, value=str(mass), size=wx.Size(self.FromDIP(90), -1))
        wt = wx.TextCtrl(self, value=str(w), size=wx.Size(self.FromDIP(50), -1))
        nm = wx.TextCtrl(self, value=str(name), size=wx.Size(self.FromDIP(170), -1))
        x = wx.Button(self, label="\u00d7", size=wx.Size(self.FromDIP(26), self.FromDIP(24)), style=wx.BU_EXACTFIT)
        x.SetToolTip("Remove this row")
        row = (m, wt, nm, x)
        self.rows.append(row)
        for c in row:
            self.table.Add(c, 0, wx.ALIGN_CENTER_VERTICAL)
        m.Bind(wx.EVT_TEXT, lambda e, r=row: (self._typed(r), e.Skip()))
        x.Bind(wx.EVT_BUTTON, lambda e, r=row: self.remove_row(r))
        for c in (m, wt, nm):
            c.Bind(wx.EVT_KEY_DOWN, lambda e, r=row: self._key(e, r))
        return row

    def _key(self, e, row):
        """Ctrl+Delete removes the row."""
        if e.GetKeyCode() == wx.WXK_DELETE and e.ControlDown():
            self.remove_row(row)
        else:
            e.Skip()

    def _typed(self, row):
        """A mass typed in the last row: a new empty row below it."""
        if row is self.rows[-1] and row[0].GetValue().strip() and len(self.rows) < self.MAX_ROWS:
            self._add_row(w=row[1].GetValue().strip() or "5")
            self._refit()

    def remove_row(self, row):
        if row not in self.rows:
            return
        if len(self.rows) == 1:  # (the last row is emptied instead)
            for c in row[:3]:
                c.ChangeValue("")
            return
        self.rows.remove(row)
        for c in row:
            self.table.Detach(c)
            c.Destroy()
        if self.rows[-1][0].GetValue().strip():
            self._add_row()
        self._refit()

    def _refit(self):
        self.Layout()
        self.Fit()

    def _event(self):
        ev = self.ev.GetSelection()
        if self.ev.GetCount() > len(self.pols) and ev == self.ev.GetCount() - 1:
            return -1
        return max(0, ev)

    def _ion_chosen(self, e):
        if self.ion.GetSelection() < len(self.names):
            self._ion = self.ion.GetSelection()

    def _carrier_shown(self):
        """Negative ions only: the carrier is a lost proton (the choice is greyed)."""
        e = self._event()
        neg = 0 <= e < len(self.pols) and self.pols[e] == "-"
        self.ion.SetSelection(len(self.names) if neg else self._ion)
        self.ion.Enable(not neg)

    def table_text(self):
        """The rows as typed (remembered for the next time): [(mass, window, name)], the empty row left out."""
        out = []
        for m, w, n, x in self.rows:
            if m.GetValue().strip() or n.GetValue().strip():
                out.append((m.GetValue().strip(), w.GetValue().strip(), n.GetValue().strip()))
        return out

    def values(self):
        import ms_find_mass as F
        rows = []
        for m, w, n, x in self.rows:
            M, W = _num(m, None), _num(w, None)
            if M is None or M <= 0:
                continue
            if W is None or W <= 0:
                W = 5.0
            if all(abs(M - r["M"]) > 1e-9 or abs(W - r["w"]) > 1e-9 for r in rows):
                rows.append({"M": float(M), "w": float(W), "name": n.GetValue().strip()})
        lo = _num(self.lo, None)
        if lo is None or lo < 0:
            lo = F.MZMIN
        return {"rows": rows, "table": self.table_text(), "carrier": self.names[self._ion], "mzmin": float(lo),
                "event": self._event()}


class Form(wx.Panel):
    """Compact form for the small dialogs: section titles and rows of a
    label and controls, no scrolling (the dialog takes its natural size)."""

    def __init__(self, parent, label_w=128, note_w=330):
        wx.Panel.__init__(self, parent, style=wx.BORDER_NONE)
        self.SetBackgroundColour(wx.Colour(C["panel"]))
        self.label_w, self.note_w = self.FromDIP(label_w), self.FromDIP(note_w)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        outer = wx.BoxSizer(wx.VERTICAL)
        outer.Add(self.sizer, 1, wx.EXPAND | wx.ALL, self.FromDIP(16))
        self.SetSizer(outer)
        self._first = True
        self.gap_row = self.FromDIP(7)
        T.keep_dip(self, "label_w", "note_w", "gap_row")

    def section(self, title, colour=None):
        s = wx.BoxSizer(wx.HORIZONTAL)
        bar = wx.Panel(self, size=wx.Size(self.FromDIP(3), self.FromDIP(15)))
        bar.SetBackgroundColour(wx.Colour(colour or T.GROUP["blue"]))
        st = wx.StaticText(self, label=title)
        st.SetFont(ui_font(9.5, 700))
        st.SetForegroundColour(wx.Colour(C["text"]))
        s.Add(bar, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, self.FromDIP(8))
        s.Add(st, 0, wx.ALIGN_CENTER_VERTICAL)
        self.sizer.Add(s, 0, wx.TOP | wx.BOTTOM, 0 if self._first else self.FromDIP(12))
        self.sizer.AddSpacer(self.FromDIP(4))
        self._first = False
        return st

    def _label(self, text, muted=False):
        st = wx.StaticText(self, label=text)
        st.SetFont(ui_font(9, 400))
        st.SetForegroundColour(wx.Colour(C["muted"] if muted else C["text"]))
        return st

    def text(self, value="", width=72, tooltip=None):
        t = wx.TextCtrl(self, value=str(value), size=wx.Size(self.FromDIP(width), -1), style=wx.TE_PROCESS_ENTER)
        t.SetFont(ui_font(9, 400))
        if tooltip:
            t.SetToolTip(tooltip)
        return t

    def choice(self, items, width=170, tooltip=None):
        c = wx.Choice(self, choices=list(items))
        c.SetFont(ui_font(9, 400))
        c.SetMinSize(wx.Size(self.FromDIP(width), -1))
        if items:
            c.SetSelection(0)
        if tooltip:
            c.SetToolTip(tooltip)
        return c

    def row(self, label, *items):
        s = wx.BoxSizer(wx.HORIZONTAL)
        lab = self._label(label or "")
        lab.SetMinSize(wx.Size(self.label_w, -1))
        s.Add(lab, 0, wx.ALIGN_CENTER_VERTICAL)
        for it in items:
            if isinstance(it, str):
                s.Add(self._label(it, muted=True), 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT | wx.RIGHT,
                      self.FromDIP(5))
            else:
                s.Add(it, 0, wx.ALIGN_CENTER_VERTICAL)
        self.sizer.Add(s, 0, wx.EXPAND | wx.TOP, self.gap_row)
        return s

    def full(self, ctrl, indent=False):
        self.sizer.Add(ctrl, 0, wx.EXPAND | wx.TOP | (wx.LEFT if indent else 0), self.gap_row)
        return ctrl

    def check(self, label, value=False, tooltip=None):
        cb = wx.CheckBox(self, label=label)
        cb.SetFont(ui_font(9, 400))
        cb.SetValue(bool(value))
        if tooltip:
            cb.SetToolTip(tooltip)
        self.sizer.Add(cb, 0, wx.TOP, self.gap_row + self.FromDIP(2))
        return cb

    def note(self, text, width=None):
        st = wx.StaticText(self, label=text)
        st.SetFont(ui_font(8.5, 400))
        st.SetForegroundColour(wx.Colour(C["muted"]))
        st._wrap_w = self.FromDIP(width) if width else self.note_w
        st.Wrap(st._wrap_w)
        self.sizer.Add(st, 0, wx.EXPAND | wx.TOP, self.gap_row)
        return st

    def set_note(self, st, text):
        st.SetLabel(text)
        st.Wrap(getattr(st, "_wrap_w", self.note_w))
        self.Layout()

    def buttons(self, *btns):
        s = wx.BoxSizer(wx.HORIZONTAL)
        for n, b in enumerate(btns):
            s.Add(b, 1, wx.EXPAND | (wx.LEFT if n else 0), self.FromDIP(6))
        self.sizer.Add(s, 0, wx.EXPAND | wx.TOP, self.gap_row + self.FromDIP(3))
        return s


def dialog_buttons(dlg, *btns, left=()):
    """Bottom bar of a small dialog: a line, then buttons (right aligned)."""
    vs = wx.BoxSizer(wx.VERTICAL)
    line = wx.Panel(dlg, size=wx.Size(-1, 1))
    line.SetBackgroundColour(wx.Colour(C["line"]))
    vs.Add(line, 0, wx.EXPAND)
    bs = wx.BoxSizer(wx.HORIZONTAL)
    for b in left:
        bs.Add(b, 0, wx.RIGHT, dlg.FromDIP(8))
    bs.AddStretchSpacer(1)
    for n, b in enumerate(btns):
        bs.Add(b, 0, wx.LEFT if n else 0, dlg.FromDIP(8))
    vs.Add(bs, 0, wx.EXPAND | wx.ALL, dlg.FromDIP(12))
    return vs


class FormDialog(wx.Dialog):
    """Small window of labelled fields (instead of the side panel).
    fields: dicts with key, label, kind ("text", "pair", "choice", "check"),
    value, and optionally unit, choices, width, tip. values() -> dict."""

    def __init__(self, parent, title, fields, note=None, ok="OK", colour=None):
        wx.Dialog.__init__(self, wx.GetTopLevelParent(parent), title=title)
        self.SetBackgroundColour(wx.Colour(C["panel"]))
        fm = self.form = Form(self)
        fm.section(title, colour or T.GROUP["blue"])
        self.ctrls = {}
        for f in fields:
            k, kind, lab = f["key"], f.get("kind", "text"), f.get("label", "")
            if kind == "text":
                c = fm.text(str(f.get("value", "")), f.get("width", 72), tooltip=f.get("tip"))
                fm.row(lab, c, *([f["unit"]] if f.get("unit") else []))
            elif kind == "pair":
                a, b = f.get("value", ("", ""))
                c = (fm.text(str(a), f.get("width", 64), tooltip=f.get("tip")),
                     fm.text(str(b), f.get("width", 64), tooltip=f.get("tip")))
                fm.row(lab, c[0], "to", c[1], *([f["unit"]] if f.get("unit") else []))
            elif kind == "choice":
                c = fm.choice(list(f["choices"]), f.get("width", 170), tooltip=f.get("tip"))
                c.SetSelection(max(0, min(int(f.get("value", 0) or 0), len(f["choices"]) - 1)))
                fm.row(lab, c)
            else:
                c = fm.check(lab, bool(f.get("value")), tooltip=f.get("tip"))
            self.ctrls[k] = (kind, c)
        if note:
            fm.note(note)
        vs = wx.BoxSizer(wx.VERTICAL)
        vs.Add(fm, 1, wx.EXPAND)
        vs.Add(dialog_buttons(self, flat(self, "Cancel", handler=lambda e: self.EndModal(wx.ID_CANCEL)),
                              flat(self, ok, "primary", handler=lambda e: self.EndModal(wx.ID_OK))), 0, wx.EXPAND)
        self.SetSizerAndFit(vs)
        for kind, c in self.ctrls.values():
            for w in (c if isinstance(c, tuple) else (c,)):
                if isinstance(w, wx.TextCtrl):
                    w.Bind(wx.EVT_TEXT_ENTER, lambda e: self.EndModal(wx.ID_OK))
        self.CentreOnParent()

    def values(self):
        out = {}
        for k, (kind, c) in self.ctrls.items():
            if kind == "text":
                out[k] = c.GetValue().strip()
            elif kind == "pair":
                out[k] = (c[0].GetValue().strip(), c[1].GetValue().strip())
            elif kind == "choice":
                out[k] = c.GetSelection()
            else:
                out[k] = c.GetValue()
        return out


def ask_form(parent, title, fields, note=None, ok="OK", colour=None):
    """values() of a FormDialog, or None when cancelled."""
    dlg = FormDialog(parent, title, fields, note, ok, colour)
    try:
        return dlg.values() if dlg.ShowModal() == wx.ID_OK else None
    finally:
        dlg.Destroy()


def _fmt_t(v):
    return "" if v is None else "%.3f" % v


def draw_badge(gc, win, text, colour, x, y, line_h, font=ui_font):
    """Small rounded label (e.g. "Poor result") in a card header: drawn at x,
    centred on the text line that starts at y and is line_h high. Returns
    its width."""
    col = wx.Colour(colour)
    gc.SetFont(font(8, 700), col)
    tw, th = gc.GetTextExtent(text)
    padx, h = win.FromDIP(7), th + win.FromDIP(3)
    by = y + (line_h - h) / 2.0
    gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(col).Width(1.0)))
    gc.SetBrush(wx.Brush(wx.Colour(col.Red(), col.Green(), col.Blue(), 30)))
    gc.DrawRoundedRectangle(x, by, tw + 2 * padx, h, win.FromDIP(3))  # Windows 11 corners, not a pill
    gc.DrawText(text, x + padx, by + (h - th) / 2.0)
    return tw + 2 * padx


_MEM_FONTS = {}
_WIN_FONTS = {}


def kept_font(points, weight=400):
    """ui_font(points, weight), made once (for window DCs: wx fits it to the DPI of the window when a DC selects
    it)."""
    f = _WIN_FONTS.get((points, weight))
    if f is None:
        f = _WIN_FONTS[(points, weight)] = ui_font(points, weight)
    return f


def own_buffered_dc(win):
    """wx.BufferedPaintDC on a buffer of win's own size, kept: the buffer wx.AutoBufferedPaintDC shares between
    windows is as large as the largest one painted and took about 4 ms each time on a large screen."""
    w, h = win.GetClientSize()
    w, h = max(w, 1), max(h, 1)
    bm = getattr(win, "_paint_buf", None)
    if bm is None or bm.GetWidth() != w or bm.GetHeight() != h:
        bm = win._paint_buf = wx.Bitmap(w, h)
    return wx.BufferedPaintDC(win, bm)


_MEM = {"dc": None, "ext": {}}   # a memory DC for measuring text; text sizes measured (font, text)


def mem_dc():
    """A memory DC kept for measuring text."""
    if _MEM["dc"] is None:
        _MEM["dc"] = wx.MemoryDC(wx.Bitmap(1, 1))
    return _MEM["dc"]


def mem_fonts(win, dc):
    """font(points, weight) for drawing or measuring on the memory DC dc as on win: made for the DPI of win (a
    memory DC has the DPI of the system), kept."""
    try:
        sc = win.GetDPI().y / float(dc.GetPPI().y or 96)
    except Exception:
        sc = 1.0

    def font(points, weight=400):
        fk = (points * sc, weight)
        f = _MEM_FONTS.get(fk)
        if f is None:
            f = _MEM_FONTS[fk] = ui_font(points * sc, weight)
        return f
    return font


class _TextFirst(object):
    """As unidec_theme.crisp(dc) (shapes with a graphics context, GDI text), the graphics context made only when a
    shape is drawn (a header with text only does not need one)."""

    def __init__(self, dc):
        self.__dict__.update(dc=dc, gc=None, _font=None, _col=wx.Colour(0, 0, 0))

    def __getattr__(self, name):
        gc = self.__dict__["gc"]
        if gc is None:
            gc = self.__dict__["gc"] = wx.GraphicsContext.Create(self.__dict__["dc"])
        return getattr(gc, name)

    def SetFont(self, font, colour=None):
        self.__dict__.update(_font=font, _col=colour if colour is not None else wx.Colour(0, 0, 0))
        self.dc.SetFont(font)

    def GetTextExtent(self, text):
        w, h = self.dc.GetTextExtent(text)
        return float(w), float(h)

    def DrawText(self, text, x, y, *args):
        if self.gc is not None:
            self.gc.Flush()
        dc = self.dc
        if self._font is not None:
            dc.SetFont(self._font)
        dc.SetBackgroundMode(wx.TRANSPARENT)
        dc.SetTextForeground(self._col)
        dc.DrawText(text, int(round(x)), int(round(y)))

    def finish(self):
        if self.gc is not None:
            self.gc.Flush()


class TextMeasure(object):
    """SetFont and GetTextExtent of a graphics context, for the sizes only: on a memory DC with fonts made for
    the DPI of win, sizes measured before kept."""

    def __init__(self, win):
        self.dc = mem_dc()
        self.font = mem_fonts(win, self.dc)
        self._f = None

    def SetFont(self, font, colour=None):
        self._f = font

    def GetTextExtent(self, text):
        key = (id(self._f), text)
        ext = _MEM["ext"].get(key)
        if ext is None:
            if len(_MEM["ext"]) > 4000:
                _MEM["ext"].clear()
            self.dc.SetFont(self._f)
            w, h = self.dc.GetTextExtent(text)
            ext = _MEM["ext"][key] = (float(w), float(h))
        return ext


def card_frame(win, gc, w, h, shadows=((3, 10), (1, 16))):
    """Shadow, white card and its edge, w x h (plot and table cards)."""
    r = win.FromDIP(12)
    gc.SetPen(wx.TRANSPARENT_PEN)
    for dy, a in shadows:
        gc.SetBrush(wx.Brush(wx.Colour(14, 28, 48, a)))
        gc.DrawRoundedRectangle(1, 1 + dy, w - 2, h - 2 - dy, r)
    gc.SetBrush(wx.Brush(wx.Colour("#FFFFFF")))
    gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(wx.Colour(C["line"])).Width(1.0)))
    gc.DrawRoundedRectangle(0.5, 0.5, w - 1.5, h - 3, r)


def table_frame(win, gc, w, h):
    """The frame of a table card (one shadow)."""
    card_frame(win, gc, w, h, ((2, 12),))


def inner_frame(win, gc, w, h):
    """A plot inside a tile of several plots (PairCard draws the card): white, no edge."""
    gc.SetPen(wx.TRANSPARENT_PEN)
    gc.SetBrush(wx.Brush(wx.Colour("#FFFFFF")))
    gc.DrawRectangle(0, 0, w, h)


class _HeaderOps(object):
    """Stands in for the graphics context while a header is laid out: what would be drawn where (text, lines,
    rounded labels) and how far to the right it reaches; the text measured as on win."""

    def __init__(self, win):
        self.m = TextMeasure(win)
        self.ops, self.right = [], 0.0

    def __getattr__(self, name):   # (anything else: drawn again as a whole)
        raise AttributeError(name)

    def SetFont(self, font, colour=None):
        self.m.SetFont(font)

    def GetTextExtent(self, text):
        return self.m.GetTextExtent(text)

    def DrawText(self, text, x, y, *args):
        self.ops.append(("t", text, int(round(x)), int(round(y))))
        self.right = max(self.right, x + self.m.GetTextExtent(text)[0])

    def StrokeLine(self, x1, y1, x2, y2):
        self.ops.append(("l", x1, y1, x2, y2))
        self.right = max(self.right, x1, x2)

    def DrawRoundedRectangle(self, x, y, w, h, r):
        self.ops.append(("r", x, y, w, h, r))
        self.right = max(self.right, x + w)

    def CreatePen(self, *args):
        return None

    def SetPen(self, *args):
        pass

    def SetBrush(self, *args):
        pass


def paint_card(win, top, frame, header, key):
    """Paints a card (a window with BG_STYLE_PAINT): frame(win, gc, w, h) draws its shadow and edge, header(gc, w,
    font) the text in its first top rows. Below them every row is the same, except the last ones (the round
    corners): the card is put together from a picture of the frame at a small size (kept while the scale and
    colours stay) and its first rows with the header (kept while the header stays), so a card made taller,
    lower or wider is painted again at once: a wider or narrower card keeps the picture of its header while the
    header is drawn the same (the frame between the header and the right edge has the same columns). The
    header is drawn on a bitmap with fonts made for the DPI of the window (a memory DC has the DPI of the
    system)."""
    w, h = win.GetClientSize()
    k = win.FromDIP(12) + 4   # the rows and columns of a round corner (with the shadow)
    if w < 2 * k + 3 or h < top + 1 + k or top < k:   # (too small: drawn as a whole)
        dc = own_buffered_dc(win)
        dc.SetBackground(wx.Brush(wx.Colour(C["bg"])))
        dc.Clear()
        gc = T.crisp(dc)
        frame(win, gc, w, h)
        header(gc, w)
        return
    bh = top + 1 + k   # the rows kept: header, one straight row, bottom
    look = (top, k, frame, win.GetDPIScaleFactor(), C["bg"], C["line"])
    if getattr(win, "_card_ref_key", None) != look:   # the frame, 2k + 2 wide
        rw = 2 * k + 2
        bm = wx.Bitmap(rw, bh)
        mdc = wx.MemoryDC(bm)
        mdc.SetBackground(wx.Brush(wx.Colour(C["bg"])))
        mdc.Clear()
        gc = T.crisp(mdc)
        frame(win, gc, rw, bh)
        del gc
        mdc.SelectObject(wx.NullBitmap)
        win._card_ref, win._card_ref_key = bm, look
    base = (look, key, C["text"], C["muted"], C["accent_text"])
    lay = getattr(win, "_card_lay", None)   # (bitmap, base, what its header draws where)
    keep = None
    if lay is not None and lay[1] == base and lay[0].GetWidth() != w:
        try:   # the header at this width: drawn the same as on the kept picture?
            ops = _HeaderOps(win)
            header(ops, w, ops.m.font)
            xr = int(ops.right + 1) + win.FromDIP(4)
            if ops.ops == lay[2] and k + 1 <= xr <= min(w, lay[0].GetWidth()) - k - 1:
                keep = xr
        except Exception:
            keep = None
    if keep is None and (lay is None or lay[1] != base or lay[0].GetWidth() != w):   # drawn at this width
        bm = wx.Bitmap(w, bh)
        mdc = wx.MemoryDC(bm)
        ref = wx.MemoryDC(win._card_ref)
        mdc.Blit(0, 0, k + 1, bh, ref, 0, 0)
        mdc.StretchBlit(k + 1, 0, w - 2 * k - 1, bh, ref, k + 1, 0, 1, bh)
        mdc.Blit(w - k, 0, k, bh, ref, k + 2, 0)
        ref.SelectObject(wx.NullBitmap)
        gc = _TextFirst(mdc)
        header(gc, w, mem_fonts(win, mdc))
        gc.finish()
        del gc
        mdc.SelectObject(wx.NullBitmap)
        ops = None
        try:
            rec = _HeaderOps(win)
            header(rec, w, rec.m.font)
            ops = rec.ops
        except Exception:
            pass
        lay = win._card_lay = (bm, base, ops)
    from ms_plot_render import copy_rect
    dc = wx.PaintDC(win)
    mdc = wx.MemoryDC(lay[0])
    # (x, width, source dc, source x, source width) across; (y, height, source y, source height) down
    if keep is None:
        cols = [(0, w, mdc, 0, w)]
    else:
        ref = wx.MemoryDC(win._card_ref)
        cols = [(0, keep, mdc, 0, keep), (keep, w - k - keep, ref, k + 1, 1), (w - k, k, ref, k + 2, k)]
    rows = [(0, top, 0, top), (h - k, k, top + 1, k)]
    if h - k > top:
        rows.append((top, h - k - top, top, 1))
    for x, pw, src, sx, sw in cols:
        for y, ph, sy, sh in rows:
            if not copy_rect(dc, x, y, pw, ph, src, sx, sy, sw, sh):
                dc.StretchBlit(x, y, pw, ph, src, sx, sy, sw, sh)
    if keep is not None:
        ref.SelectObject(wx.NullBitmap)
    mdc.SelectObject(wx.NullBitmap)


# ==========================================================================
# plotting
# ==========================================================================
class PlotCard(wx.Panel):
    """White rounded card with a title and one matplotlib axes.

    Mouse: left drag = range (or zoom, see mode), Shift + drag = second range,
    Ctrl + drag = zoom box, click = pick, double click = full view, right
    click = context menu supplied by the owner. The wheel scrolls the column
    of tiles; Ctrl + wheel zooms the time or m/z axis."""

    own_refresh = True   # (paints itself again when resized: place_panes)

    def __init__(self, parent, title, subtitle="", mode="range"):
        wx.Panel.__init__(self, parent, style=wx.BORDER_NONE)
        _mpl_setup()
        from matplotlib.figure import Figure
        from ms_plot_render import canvas_class
        self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
        self.SetBackgroundColour(C["bg"])
        self.title, self.subtitle, self.mode = title, subtitle, mode
        self.fig = Figure(figsize=(4, 2.6), facecolor="white")
        self.canvas = canvas_class()(self, -1, self.fig)
        self.canvas.SetMinSize(wx.Size(20, 20))  # the tile, not the figure, sets the size
        self.ax = self.fig.add_subplot(111)
        self.full = None  # (x0, x1) full view
        self.xlock = None  # (x0, x1) kept by "m/z range...": double click and new spectra come back to it
        self.on_range = None  # callback(x0, x1, shift)
        self.on_pick = None  # callback(x, y)
        self.on_menu = None  # callback(x, y) -> list of (label, func)
        self.on_view = None  # callback() after zoom
        self.readout_fmt = None  # callback(x, y) -> text shown at the top right while hovering
        self.tip_fmt = None  # callback(x, y) -> tooltip over the plot at the pointer ("" for none)
        self._plot_tip = ""
        self.readout = ""
        self.badge = None  # (text, colour, tooltip) of a warning pill after the title, or None
        self.legend_items = []  # [(colour, label)] shown in the title row
        self.extra_scale = []  # [(x, y)] data that are not Line2D (e.g. sticks) for the y autoscale
        self.label_room_pt = 16.0  # free space above the tallest peak (points) for its label
        self.relabel_fn = None  # callback(): labels placed again for the present size (images of another size)
        self._auto = None  # (pad, ylim) of the last automatic y range
        self.image_name = None  # callback() -> suggested file name (no extension) for Save image
        self.tool_getter = None  # callback() -> active tool of the tab
        self.tools_supported = set()  # tools this plot reacts to (zoom and pan always)
        self.on_tool = None  # callback(tool, "click"|"drag", a, b)
        self.ylock = None  # y range fixed by the zoom box or panning
        self.y_fixed = None  # (lo or None, hi or None): ends of the y range set by the owner (graph properties)
        self.live = []  # artists redrawn by blitting (fast updates while dragging)
        self._bgcache = None
        self._blit_ok = True
        self._press = None
        self._drag_art = None
        self.margins = [0.84, 0.52, 0.16, 0.14]  # left, bottom, right, top (inches)
        self.pair = None  # the PairCard holding this plot with others (one tile), or None
        self.hide_x = False  # no x numbers and title (a plot above another of the same tile, which has them)
        # in a tile of several plots (pair): no title row; the y title is side_name (e.g. "(+) TIC") and the
        # power of ten of the axis, in side_colour; the title row's text is a small note in the plot
        self.side_name = ""
        self.side_colour = None
        self._side_txt = None
        self._side_n = 0
        self._ylabel_orig = ""
        self._unit = ""
        self._info_box = None
        self._info_key = None
        self._ro_art = None
        self._export = False  # an image being made: the x axis shown in full
        self._syncing = False
        s = wx.BoxSizer(wx.VERTICAL)
        self._hh = self.FromDIP(30)
        s.AddSpacer(self._hh)
        s.Add(self.canvas, 1, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, self.FromDIP(8))
        self.SetSizer(s)
        self.Bind(wx.EVT_PAINT, self._paint)
        self.Bind(wx.EVT_SIZE, self._size)
        for ev, fn in (("button_press_event", self._down), ("button_release_event", self._up),
                       ("motion_notify_event", self._move), ("scroll_event", self._scroll),
                       ("resize_event", self._resized),
                       ("figure_leave_event", lambda e: (self._set_readout(""), self._set_plot_tip(""))),
                       ("draw_event", self._on_draw), ("key_press_event", self._hotkey)):
            self.canvas.mpl_connect(ev, fn)
        style_axes(self.ax)
        self.place()
        # the wheel scrolls the column of tiles; only Ctrl + wheel reaches
        # matplotlib (zoom). The matplotlib handler is replaced, not stacked,
        # so the order of the handlers does not matter.
        self.canvas.Unbind(wx.EVT_MOUSEWHEEL)
        self.canvas.Bind(wx.EVT_MOUSEWHEEL, self._wheel)
        self.Bind(wx.EVT_MOUSEWHEEL, lambda e: forward_wheel(self, e))

    def _wheel(self, e):
        if e.ControlDown():
            handler = getattr(self.canvas, "_on_mouse_wheel", None)
            if handler is not None:
                handler(e)  # matplotlib scroll_event: zoom
        else:
            forward_wheel(self, e)

    # ----------------------------------------------------------- appearance
    def set_title(self, title, subtitle=None):
        old = (self.title, self.subtitle)
        self.title = title
        if subtitle is not None:
            self.subtitle = subtitle
        self._update_header()
        self.Refresh()
        if self.pair is not None and old != (self.title, self.subtitle):
            self.draw()  # (the note in the plot)

    def set_badge(self, text=None, colour=None, tip=""):
        """Warning pill after the title (None removes it); tip is added to
        the tooltip of the card."""
        badge = (text, colour or T.GROUP["red"], tip or "") if text else None
        if badge != self.badge:
            self.badge = badge
            self._update_header()
            self.Refresh()
            if self.pair is not None:
                self.draw()

    def set_shared(self):
        """In a tile of several plots (PairCard): no title row and no edge around the plot (place: its y title,
        a note of the title row inside the plot, the unit at the end of the x axis)."""
        it = self.GetSizer().GetItem(1)  # the plot
        it.SetFlag(wx.EXPAND | wx.LEFT | wx.RIGHT)
        self._hdr_key = None
        self._update_header()
        self.place()

    def _badge_width(self, gc, font=ui_font):
        if not self.badge:
            return 0
        gc.SetFont(font(8, 700), wx.Colour(self.badge[1]))
        return self.FromDIP(10) + gc.GetTextExtent(self.badge[0])[0] + 2 * self.FromDIP(7)

    def _size(self, e):
        self._update_header()
        self.Layout()
        self.Refresh()
        e.Skip()

    def _paint(self, e):
        if self.pair is not None:  # (no title row)
            paint_card(self, self._hh, inner_frame, lambda gc, w, font=None: None, ("shared",))
            return
        paint_card(self, self._hh, card_frame, self._draw_header, self._header_key())

    def _view_changed(self):
        """The x range was changed by the user: the owner is told, and the other plots of the tile follow."""
        if self.on_view:
            self.on_view()
        if self.pair is not None and not self._syncing:
            self.pair.sync_x(self)

    def _header_key(self):
        """What the header shows (text, legend, badge, read-out) and its scale."""
        return (self.title, self.subtitle, self.badge, tuple(self.legend_items), self.readout,
                self.readout_fmt is None, self.GetDPIScaleFactor())

    def _header_layout(self, gc, w, font=ui_font):
        """Positions of title, legend and subtitle: one line if they fit,
        otherwise the subtitle (and the legend if needed) on a second line.
        Returns (height, title_y, line2_y, legend_on_line2, sub_on_line2)."""
        gc.SetFont(font(10, 700), wx.Colour(C["text"]))
        tw, th = gc.GetTextExtent(self.title or "Xg")
        tw += self._badge_width(gc, font)
        gc.SetFont(font(9, 500), wx.Colour(C["text"]))
        leg = sum(self.FromDIP(29) + gc.GetTextExtent(lab)[0] for col, lab in self.legend_items)
        gc.SetFont(font(9, 400), wx.Colour(C["muted"]))
        sw = gc.GetTextExtent(self.subtitle)[0] if self.subtitle else 0
        rw = self.FromDIP(170) if self.readout_fmt is not None else 0  # room kept for the read-out
        room = w - self.FromDIP(14) - tw - self.FromDIP(12) - rw - self.FromDIP(12)
        top = self.FromDIP(8)
        line = th + self.FromDIP(4)
        if leg + sw <= room:
            return top + line + self.FromDIP(8), top, None, False, False
        if leg <= room:
            return top + 2 * line + self.FromDIP(6), top, top + line, False, True
        return top + 2 * line + self.FromDIP(6), top, top + line, True, True

    def _draw_header(self, gc, w, font=ui_font):
        hh, y, y2, leg2, sub2 = self._header_layout(gc, w, font)
        gc.SetFont(font(10, 700), wx.Colour(C["text"]))
        tw, th = gc.GetTextExtent(self.title)
        x = self.FromDIP(14)
        gc.DrawText(self.title, x, y)
        if self.badge:
            tw += self.FromDIP(10) + draw_badge(gc, self, self.badge[0], self.badge[1], x + tw + self.FromDIP(10), y,
                                                th, font)
        rw = 0
        if self.readout:
            gc.SetFont(font(9, 600), wx.Colour(C["accent_text"]))
            rw, rh = gc.GetTextExtent(self.readout)
            gc.DrawText(self.readout, w - rw - self.FromDIP(16), y + (th - rh) / 2.0)
            rw += self.FromDIP(24)
        limit = w - rw - self.FromDIP(12)
        xx = x + tw + self.FromDIP(12)
        ly = y
        if leg2:
            xx, ly = x, y2
            limit = w - self.FromDIP(12)
        # legend: coloured dashes and names
        gc.SetFont(font(9, 500), wx.Colour(C["text"]))
        for col, lab in self.legend_items:
            lw_, lh_ = gc.GetTextExtent(lab)
            need = self.FromDIP(18) + lw_ + self.FromDIP(12)
            if xx + need > limit:
                break
            cy = ly + th / 2.0
            gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(wx.Colour(col)).Width(2.5 * self.FromDIP(100) / 100.0)
                                   .Cap(wx.CAP_ROUND)))
            gc.StrokeLine(xx, cy, xx + self.FromDIP(12), cy)
            gc.DrawText(lab, xx + self.FromDIP(17), ly + (th - lh_) / 2.0)
            xx += need
        if self.subtitle:
            gc.SetFont(font(9, 400), wx.Colour(C["muted"]))
            sub = self.subtitle
            sy = ly
            if sub2:
                sy = y2
                if not leg2:
                    xx = x
                limit = w - self.FromDIP(12)
            room = limit - xx
            if gc.GetTextExtent(sub)[0] > room:
                while sub and gc.GetTextExtent(sub + "\u2026")[0] > room:
                    sub = sub[:-1]
                sub = (sub.rstrip() + "\u2026") if sub else ""
            sw, sh = gc.GetTextExtent(sub)
            gc.DrawText(sub, xx, sy + (th - sh) / 2.0)

    def _update_header(self):
        """Header height follows its text (one or two lines)."""
        try:
            if self.pair is not None:  # no title row: a thin strip above the plot
                hh = self.FromDIP(1)
                if hh != self._hh:
                    self._hh = hh
                    self.GetSizer().GetItem(0).AssignSpacer(wx.Size(0, hh))
                    self.Layout()
                return
            w = self.GetClientSize()[0]
            if w < 50:
                return
            key = (w, self._header_key())
            if key == getattr(self, "_hdr_key", None):
                return   # (the same width and text: the same header; e.g. a tile made taller)
            self._hdr_key = key
            m = TextMeasure(self)
            hh = int(round(self._header_layout(m, w, m.font)[0]))
            if abs(hh - self._hh) > 1:
                self._hh = hh
                self.GetSizer().GetItem(0).AssignSpacer(wx.Size(0, hh))
                self.Layout()
            tip = self.title + ((": " + self.subtitle) if self.subtitle else "")
            if self.badge and self.badge[2]:
                tip += "\n" + self.badge[2]
            if tip != getattr(self, "_tip", None):
                self._tip = tip
                self.SetToolTip(tip)
        except RuntimeError:
            pass

    def set_legend(self, items):
        items = list(items)
        if items != self.legend_items:
            self.legend_items = items
            self._update_header()
            self.RefreshRect(wx.Rect(0, 0, self.GetClientSize()[0], self._hh))
            if self.pair is not None:
                self.draw()

    def _set_readout(self, text):
        if text != self.readout:
            self.readout = text
            if self.pair is not None:
                self._shared_readout(text)
                return
            w = self.GetClientSize()[0]
            self.RefreshRect(wx.Rect(0, 0, w, self._hh))

    def _shared_readout(self, text):
        """The read-out of a plot without a title row: at the top left inside the plot, drawn over the stored
        picture of the plot (not in images)."""
        try:
            ax = self.ax
            art = self._ro_art
            if art is None or art not in ax.texts:
                if not text:
                    return
                art = ax.text(0.006, 0.975, "", transform=ax.transAxes, ha="left", va="top", fontsize=8,
                              color=C["accent_text"], fontweight="bold", zorder=12,
                              bbox=dict(boxstyle="square,pad=0.25", fc="white", ec="none", alpha=0.85))
                art._ui_only = True
                self._ro_art = art
                self.add_live(art)
            art.set_text(text)
            art.set_visible(bool(text))
            self.blit_live()
        except RuntimeError:
            pass

    # ------------------------------------------- a plot in a tile of several plots
    SHARED_L = 0.78  # left margin (inches): the y title and the numbers
    SHARED_DX = 0.64  # the middle of the y title from the frame (inches)
    SHARED_B = (0.06, 0.27)  # bottom margin without and with the x numbers
    SHARED_T = 0.05  # top margin

    def _side_base(self):
        name = self.side_name or self.title or ""
        if "(%)" in (self._ylabel_orig or ""):
            name += "  (%)"
        return name

    def _side_power(self, n):
        """The y title with the power of ten of the axis (SideFormatter, at every drawing)."""
        try:
            self._side_n = n = int(n or 0)
            txt = self._side_base() + (("  ×10" + ("%d" % n).translate(_SUP)) if n else "")
            lab = self.ax.yaxis.label
            if lab.get_text() != txt:
                lab.set_text(txt)
            self._side_txt = txt
        except Exception:
            pass

    def _info_items(self):
        """[(colour, text)] of the title row shown as a note in the plot: the legend (except the trace named
        in the y title), the warning and the subtitle."""
        side = " " + (self.side_name or "") + " "
        items = [(col, lab) for col, lab in self.legend_items if lab and (" " + lab + " ") not in side]
        if self.badge:
            items.append((self.badge[1], self.badge[0]))
        if self.subtitle:
            items.append((C["muted"], " ".join(self.subtitle.split())))
        return items

    def _unit_shown(self, on):
        u = getattr(self.ax, "_unit_art", None)
        if u is not None:
            u.set_visible(bool(on) and u in self.ax.texts)

    def _shared_look(self, show_x, w_in):
        """y title (name, power of ten), unit at the right end of the x axis, more x numbers, the note at the
        top right: made again after every new drawing of the plot (style_axes, cla)."""
        from matplotlib.ticker import AutoMinorLocator
        from matplotlib.transforms import ScaledTranslation
        ax = self.ax
        # y title
        lab = ax.yaxis.label
        if lab.get_text() != self._side_txt:  # (a new drawing: the axis title of style_axes)
            self._ylabel_orig = lab.get_text()
            self._side_txt = None
        if not getattr(ax.yaxis.get_major_formatter(), "_side", False):
            ax.yaxis.set_major_formatter(_side_formatter(self))
        ax.yaxis.get_offset_text().set_visible(False)
        self._side_power(self._side_n)
        lab.set_fontsize(9.5)
        lab.set_fontweight("bold")
        lab.set_color(self.side_colour or INK)
        lab.set_verticalalignment("center")
        lab.set_visible(True)
        ax.yaxis.set_label_coords(0.0, 0.5, transform=ax.transAxes + ScaledTranslation(
            -self.SHARED_DX, 0.0, self.fig.dpi_scale_trans))
        # x axis: more numbers, the unit on their line at the right end
        loc = ax.xaxis.get_major_locator()
        if hasattr(loc, "dense"):
            loc.dense = True
        ml = ax.xaxis.get_minor_locator()
        if isinstance(ml, AutoMinorLocator) and getattr(ml, "ndivs", None) != 2:
            ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        if type(ax.xaxis.get_major_formatter()).__name__ == "ScalarFormatter":
            ax.xaxis.set_major_formatter(_unit_formatter())
        xl = ax.xaxis.label.get_text()
        if xl:
            self._unit = xl
            ax.xaxis.label.set_text("")
        u = getattr(ax, "_unit_art", None)
        if (u is None or u not in ax.texts) and self._unit:
            u = ax.text(1.0, 0.0, self._unit, ha="right", va="top", fontsize=9, color=INK, clip_on=False,
                        style="italic" if self._unit == "m/z" else "normal",
                        transform=ax.transAxes + ScaledTranslation(0.0, -6.5 / 72.0, self.fig.dpi_scale_trans))
            u._unit_mark = True
            ax._unit_art = u
        elif u is not None and self._unit and u.get_text() != self._unit:
            u.set_text(self._unit)
            u.set_fontstyle("italic" if self._unit == "m/z" else "normal")
        self._unit_shown(show_x)
        if self._export:
            return  # (an image being made: the note stays as it is, hidden)
        # the note: what the title row said
        items = self._info_items()
        try:
            i = self.pair.cards.index(self)
            if i > 0 and self.pair.cards[i - 1]._info_items() == items:
                items = []  # (the same as on the plot above)
        except (ValueError, AttributeError):
            pass
        room = max(40.0, 0.86 * (w_in - self.SHARED_L - self.margins[2]) * 72.0)  # (points)
        drop = bool(getattr(self, "top_tags", False))  # charge tags along the top: the note goes below them
        key = (tuple(items), round(room), drop)
        box = self._info_box
        if box is not None and (key != self._info_key or box not in ax.get_children()):
            try:
                box.remove()
            except Exception:
                pass
            box = self._info_box = None
        self._info_key = key
        if box is None and items:
            from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, TextArea
            parts, used = [], 0.0
            for j, (col, txt) in enumerate(items):
                sep = 5.0 + _text_w_pt("·", 8) if j else 0.0
                tw = _text_w_pt(txt, 8)
                if used + sep + tw > room:  # (a narrow plot: the text shortened)
                    while txt and used + sep + _text_w_pt(txt + "…", 8) > room:
                        txt = txt[:-1]
                    if not txt.strip():
                        break
                    txt = txt.rstrip() + "…"
                    tw = room
                if j:
                    parts.append(TextArea("·", textprops=dict(color=C["muted"], fontsize=8)))
                parts.append(TextArea(txt, textprops=dict(color=col, fontsize=8)))
                used += sep + tw
                if used >= room:
                    break
            if parts:
                kw = {}
                if drop:
                    h_pt = max(1.0, ax.get_position().height * self.fig.get_figheight() * 72.0)
                    kw = dict(bbox_to_anchor=(1.0, 1.0 - min(0.4, 13.0 / h_pt)), bbox_transform=ax.transAxes)
                box = AnchoredOffsetbox(loc="upper right", child=HPacker(children=parts, align="baseline", pad=0,
                                                                          sep=5), pad=0.2, borderpad=0.3,
                                        frameon=True, **kw)
                box.patch.set_facecolor("white")
                box.patch.set_edgecolor("none")
                box.patch.set_alpha(0.78)
                box.set_zorder(8)
                box._ui_only = True  # (as the title row: not in images)
                ax.add_artist(box)
                self._info_box = box

    def _side_title_at(self, dx=None):
        """In an image: the y title just left of the y numbers, dx inches from the frame (measured when None;
        the plots of one image use the largest, so their titles line up). Returns dx."""
        from matplotlib.transforms import ScaledTranslation
        try:
            if dx is None:
                r = self.fig.canvas.get_renderer()
                ticks = self.ax.yaxis._update_ticks()
                ws = [t.label1.get_window_extent(r).width for t in ticks
                      if t.label1.get_visible() and t.label1.get_text()]
                lab = self.ax.yaxis.label
                tick = ticks[0].get_tick_padding() + ticks[0].get_pad() if ticks else 6.5   # (points)
                dx = (max(ws or [0.0]) + 0.5 * lab.get_window_extent(r).width) / self.fig.dpi + \
                    (tick + 0.45 * lab.get_fontsize()) / 72.0
            self.ax.yaxis.set_label_coords(0.0, 0.5, transform=self.ax.transAxes + ScaledTranslation(
                -dx, 0.0, self.fig.dpi_scale_trans))
            return dx
        except Exception as ex:
            _log("y title of an image: %s" % ex)
            return None

    def _place_shared(self, w, h):
        show_x = self._export or not self.hide_x
        l, r, t = self.SHARED_L, self.margins[2], self.SHARED_T
        b = self.SHARED_B[1] if show_x else self.SHARED_B[0]
        self.ax.tick_params(axis="x", labelbottom=show_x)
        self.ax.xaxis.label.set_visible(False)
        self._shared_look(show_x, w)
        self.ax.set_position([l / w, b / h, max(0.05, 1 - (l + r) / w), max(0.05, 1 - (b + t) / h)])
        self._fit_notes(w * (1 - (l + r) / w) * self.fig.dpi)

    def _set_plot_tip(self, tip):
        if tip != self._plot_tip:
            self._plot_tip = tip
            try:
                if tip:
                    self.canvas.SetToolTip(tip)
                else:
                    self.canvas.UnsetToolTip()
            except Exception:
                pass

    # ------------------------------------------------------ live layer
    def reset(self):
        """Clear the axes (and the live layer) before a full replot."""
        self.ax.cla()
        self.live = []
        self.extra_scale = []
        self._bgcache = None

    def add_live(self, art):
        if art is not None:
            art.set_animated(self._blit_ok)
            self.live.append(art)
        return art

    def remove_live(self, arts):
        for a in arts:
            try:
                self.live.remove(a)
            except ValueError:
                pass
            try:
                a.remove()
            except Exception:
                pass

    def _on_draw(self, e):
        if not self.live or not self._blit_ok:
            self._bgcache = None  # an older stored image must not come back
            return
        try:
            self._bgcache = self.canvas.copy_from_bbox(self.fig.bbox)
        except Exception as ex:
            _log("blitting off (%s)" % ex)
            self._no_blit()
            return
        # the wx bitmap is updated after this event: draw the live layer afterwards
        wx.CallAfter(self._safe_blit)

    def _safe_blit(self):
        try:
            self.blit_live()
        except RuntimeError:
            pass  # window closed

    def _no_blit(self):
        self._blit_ok = False
        self._bgcache = None
        for a in self.live:
            a.set_animated(False)

    def blit_live(self):
        if not self._blit_ok or self._bgcache is None:
            self.canvas.draw_idle()
            return
        try:
            self.canvas.restore_region(self._bgcache)
            for a in self.live:
                if a.get_visible() and a.axes is not None:
                    self.ax.draw_artist(a)
            self.canvas.blit(self.fig.bbox)
        except Exception as ex:
            _log("blitting off (%s)" % ex)
            self._no_blit()
            self.canvas.draw_idle()

    # ------------------------------------------------------ image export
    def _render(self, fmt="png", dpi=None, target=None, measure=False, align=None, shared=False, height_in=None):
        """The plot as an image: at the journal size of the image settings
        (width x height in cm, text size in points), or as on screen. For the image of a whole tile
        (PairCard): shared, a plot above another has no x numbers; measure returns its margins (left, right,
        bottom, inches; None as on screen) and the place of its y title instead; align (left, right, y title)
        sets them; height_in its height."""
        import io
        from matplotlib.text import Text
        cfg = export_settings()
        dpi = dpi or cfg["dpi"]
        self._export = not (shared and self.hide_x)
        self.ax.tick_params(axis="x", labelbottom=self._export or not self.hide_x)
        self._unit_shown(self._export or not self.hide_x)  # (a plot of a tile of several: the unit of x)
        ui = [a for a in self.fig.findobj() if getattr(a, "_ui_only", False) and a.get_visible()]
        for a in ui:
            a.set_visible(False)
        for a in self.live:
            a.set_animated(False)
        # saved images always carry the axis titles, even from a short tile
        self.ax.xaxis.label.set_visible(not (shared and self.hide_x))
        self.ax.yaxis.label.set_visible(True)
        out = target if target is not None else io.BytesIO()
        old_size = self.fig.get_size_inches().copy()
        old_pos = [(a, a.get_position()) for a in self.fig.axes]
        old_ylim, old_auto = self.ax.get_ylim(), self._auto
        fonts = []
        relabeled = False
        titles = None
        try:
            if cfg["mode"] == "screen":
                if measure:
                    return None
                self.fig.savefig(out, format=fmt, dpi=dpi, facecolor="white", bbox_inches="tight", pad_inches=0.06)
            else:
                k = float(cfg["text"]) / 10.0  # on screen the axis titles are 10 pt
                for t in self.fig.findobj(Text):
                    fonts.append((t, t.get_fontsize()))
                    t.set_fontsize(max(6.5, t.get_fontsize() * k))
                self.fig.set_size_inches(cfg["width_cm"] / 2.54, height_in or cfg["height_cm"] / 2.54, forward=False)
                titles = _titles_state(self.fig)
                fit_axes(self.fig)
                fit_titles(self.fig)
                fit_axes(self.fig)
                # an automatic y range is made again for the image, whose
                # proportions differ from the tile (no empty band at the top)
                a = self._auto
                if a is not None and self.ylock is None and np.allclose(self.ax.get_ylim(), a[1]):
                    self.autoscale_y(a[0])
                    fit_axes(self.fig)
                side_dx = None
                if self.pair is not None:  # the y title next to the numbers (in an image of the tile: lined up)
                    side_dx = self._side_title_at(align[2] if align is not None and len(align) > 2 else None)
                    fit_axes(self.fig)
                if self.relabel_fn is not None:  # labels placed for the image (they overlapped in it)
                    relabeled = True
                    for _ in range(2):
                        h0 = self.ax.get_position().height
                        self.relabel_fn()
                        for t in self.fig.findobj(Text):
                            if getattr(t, "_dec_label", False):
                                fonts.append((t, t.get_fontsize()))
                                t.set_fontsize(max(6.5, t.get_fontsize() * k))
                        fit_axes(self.fig)
                        # labels that changed the height of the plot (a title on more lines than on screen)
                        # are placed once more for it: names placed for a taller plot ran into each other
                        if abs(self.ax.get_position().height - h0) <= 0.03 * h0:
                            break
                W = self.fig.get_size_inches()[0]
                pos = self.ax.get_position()
                if measure:
                    return pos.x0 * W, (1.0 - pos.x1) * W, pos.y0 * self.fig.get_size_inches()[1], side_dx
                if align is not None:  # the plots of a tile one below the other: their frames line up
                    self.ax.set_position([align[0] / W, pos.y0, max(0.05, 1.0 - (align[0] + align[1]) / W),
                                          pos.height])
                self.fig.savefig(out, format=fmt, dpi=dpi, facecolor="white")
        finally:
            self._export = False
            if titles is not None:
                _restore_titles(titles)  # their texts (sizes: next line, as before the image)
            for t, fs in fonts:
                t.set_fontsize(fs)
            self.fig.set_size_inches(old_size, forward=False)
            for a, pos in old_pos:
                a.set_position(pos)
            self.ax.set_ylim(*old_ylim)
            self._auto = old_auto
            for a in ui:
                a.set_visible(True)
            for a in self.live:
                a.set_animated(self._blit_ok)
            if relabeled:
                try:
                    self.relabel_fn()  # back to the labels of the tile
                except Exception as ex:
                    _log("labels after the image: %s" % ex)
            self.place()
            self.canvas.draw_idle()
        return out

    def shown(self):
        """Is this tile on screen (not closed, not in a hidden row)?"""
        try:
            return bool(self.IsShownOnScreen()) or (self.IsShown() and self.GetParent().IsShown())
        except RuntimeError:
            return False

    def snapshot(self, path, w_cm, h_cm, text_scale=0.75, dpi=300, relabel=None, hide=None, legend=True,
                 legend_keep=None, min_scale=0.6, h_range=None):
        """The plot as it is on screen (same view, traces, labels, markers, overlays and shaded
        ranges) as an image for a report figure of w_cm x h_cm whose text is text_scale x the
        text on screen. The figure is laid out at (w_cm, h_cm) / text_scale, where every text keeps
        its size on screen, and saved at dpi x text_scale: placed at w_cm in the report it has the
        report's resolution and type size. relabel(): places the labels again for this size (the
        labels of a deconvolution result depend on the room in points), and once more for the tile
        afterwards. hide(artist) -> True: left out of the image (e.g. shaded ranges not wanted).
        The legend of the title row (trace names) is drawn into the plot (legend_keep(label) -> False:
        an entry left out).

        A tile wider on screen than the report figure is drawn with smaller text (down to min_scale
        x the screen text), so that the labels placed for its size are those of the screen (a
        zoomed mass spectrum keeps its isotope labels). h_range (lo, hi): the height follows the
        proportions of the tile on screen within these limits (cm) instead of h_cm."""
        from matplotlib.lines import Line2D
        if relabel is None:
            relabel = self.relabel_fn
        k = max(0.3, float(text_scale))
        try:
            ws, hs = [float(v) * 2.54 for v in self.fig.get_size_inches()]  # the tile on screen (cm)
            if ws > 1.0:
                k = max(min(k, float(min_scale)), min(k, w_cm / ws))
                if h_range is not None and hs > 0.5:
                    h_cm = max(float(h_range[0]), min(float(h_range[1]), w_cm * hs / ws))
        except Exception:
            pass
        ui = [a for a in self.fig.findobj() if a.get_visible() and (getattr(a, "_ui_only", False) or
                                                                    (hide is not None and hide(a)))]
        for a in ui:
            a.set_visible(False)
        for a in self.live:
            a.set_animated(False)
        xl_vis, yl_vis = self.ax.xaxis.label.get_visible(), self.ax.yaxis.label.get_visible()
        self._export = True
        self.ax.tick_params(axis="x", labelbottom=True)
        self._unit_shown(True)
        self.ax.xaxis.label.set_visible(True)
        self.ax.yaxis.label.set_visible(True)
        old_size = self.fig.get_size_inches().copy()
        old_pos = [(a, a.get_position()) for a in self.fig.axes]
        old_ylim, old_auto = self.ax.get_ylim(), self._auto
        leg = None
        titles = None
        try:
            items = [it for it in (self.legend_items or []) if it and it[1] and
                     (legend_keep is None or legend_keep(it[1]))]
            if legend and items:
                hs = [Line2D([], [], color=c, lw=1.4) for c, _ in items]
                leg = self.ax.legend(hs, [t for _, t in items], loc="upper right", fontsize=8, frameon=False,
                                     handlelength=1.4, borderaxespad=0.3, labelspacing=0.25)
            self.fig.set_size_inches(w_cm / 2.54 / k, h_cm / 2.54 / k, forward=False)
            self.place()
            self.ax.xaxis.label.set_visible(True)
            self.ax.yaxis.label.set_visible(True)
            titles = _titles_state(self.fig)
            if self.pair is not None:
                self._side_title_at()
            fit_axes(self.fig)
            fit_titles(self.fig)
            fit_axes(self.fig)
            a = self._auto
            if a is not None and self.ylock is None and np.allclose(self.ax.get_ylim(), a[1]):
                self.autoscale_y(a[0])  # the room above the peaks for their labels is in points
                fit_axes(self.fig)
            if relabel is not None:
                relabel()
                fit_axes(self.fig)
            self.fig.savefig(path, dpi=dpi * k, facecolor="white")
        finally:
            if leg is not None:
                try:
                    leg.remove()
                except Exception:
                    pass
            if titles is not None:
                _restore_titles(titles)
            self.fig.set_size_inches(old_size, forward=False)
            for a, pos in old_pos:
                a.set_position(pos)
            self.ax.set_ylim(*old_ylim)
            self._auto = old_auto
            for a in ui:
                a.set_visible(True)
            for a in self.live:
                a.set_animated(self._blit_ok)
            self._export = False
            self.ax.xaxis.label.set_visible(xl_vis)
            self.ax.yaxis.label.set_visible(yl_vis)
            if relabel is not None:
                try:
                    relabel()
                except Exception as ex:
                    _log("labels after the report image: %s" % ex)
            self.place()
            self.canvas.draw_idle()
        return path

    def copy_image(self):
        """PNG with its resolution (Word and PowerPoint paste it at the set
        size in cm) plus a bitmap (DIB) for other programs. On Windows the
        clipboard is written directly (Win32), which every program reads;
        elsewhere through wx."""
        try:
            cfg = export_settings()
            data = self._render("png").getvalue()
            if not data.startswith(b"\x89PNG"):
                raise RuntimeError("the image could not be drawn")
            how = copy_png_to_clipboard(self, data)
            if cfg["mode"] == "screen":
                msg = "Image copied (%d dpi, as on screen): paste it with Ctrl+V" % cfg["dpi"]
            else:
                msg = "Image copied: %g x %g cm, text %g pt, %d dpi. Paste it with Ctrl+V" % (
                    cfg["width_cm"], cfg["height_cm"], cfg["text"], cfg["dpi"])
            _status(self, msg)
            toast(self, msg)
            _log("copy image: %d bytes, %s" % (len(data), how))
        except Exception as ex:
            _report_error(self, "The image could not be copied", ex)

    def save_image(self):
        try:
            frame = wx.GetTopLevelParent(self)
            name = None
            try:
                name = self.image_name() if self.image_name else None
            except Exception:
                name = None
            if not name:
                base = frame.base_name() if hasattr(frame, "base_name") else "plot"
                name = base + "_" + "".join(ch if ch.isalnum() else "_" for ch in self.title.lower()).strip("_")
            name = safe_file_name(name)
            wild = ("PNG image, 300 dpi (*.png)|*.png|TIFF image, 300 dpi (*.tif)|*.tif|"
                    "PDF, vector (*.pdf)|*.pdf|SVG, vector (*.svg)|*.svg")
            if hasattr(frame, "ask_save"):
                path = frame.ask_save("Save image", name + ".png", wild)
            else:
                path = ask_save_file(frame, "Save image", name + ".png", wild)
            if not path:
                return
            fmt = os.path.splitext(path)[1].lower().lstrip(".")
            fmt = {"tif": "tiff", "jpg": "jpeg"}.get(fmt, fmt)
            if fmt not in ("png", "tiff", "pdf", "svg", "jpeg"):
                fmt = "png"
                path += ".png"
            data = self._render(fmt).getvalue()
            with open(path, "wb") as fh:
                fh.write(data)
            if not os.path.isfile(path) or os.path.getsize(path) == 0:
                raise RuntimeError("nothing was written to " + path)
            _status(self, "Saved " + path)
            toast(self, "Saved %s\nin %s" % (os.path.basename(path), os.path.dirname(path)))
            _log("save image: %s (%d bytes)" % (path, len(data)))
        except Exception as ex:
            _report_error(self, "The image could not be saved", ex)

    def place(self):
        """Axes position from margins in inches; short tiles get a compact
        layout (no axis titles, smaller margins) so the plot keeps room."""
        try:
            w, h = self.fig.get_size_inches()
            if w <= 0.3 or h <= 0.3:
                return
            if self.pair is not None:
                self._place_shared(w, h)
                return
            l, b, r, t = self.margins
            if self.ax.get_title(loc="right") or self.ax.get_title(loc="left") or self.ax.get_title():
                t = max(t, 0.30)  # room for a caption above the plot
            compact = h < 1.45 and self.mode != "map"
            if compact:
                l, b, t = min(l, 0.68), 0.34, min(t, 0.30)
            self.ax.xaxis.label.set_visible(not compact)
            self.ax.yaxis.label.set_visible(not compact)
            if self.hide_x:
                # the plot below in the tile has the x numbers and title (an image of this plot alone has them)
                self.ax.tick_params(axis="x", labelbottom=self._export)
                if not self._export:
                    b = 0.10
                    self.ax.xaxis.label.set_visible(False)
            self.ax.set_position([l / w, b / h, max(0.05, 1 - (l + r) / w), max(0.05, 1 - (b + t) / h)])
            self._fit_notes(w * (1 - (l + r) / w) * self.fig.dpi)
        except Exception:
            pass

    def _fit_notes(self, width_px):
        """A note in the middle of the plot (an empty plot's hint) on two lines when one is wider than the plot,
        so it does not run over the axis numbers of a narrow tile."""
        notes = [tx for tx in self.ax.texts if tx.get_transform() is self.ax.transAxes and
                 tuple(tx.get_position()) == (0.5, 0.5) and tx.get_ha() == "center"]
        if not notes:
            return
        r = self.canvas.get_renderer()
        for tx in notes:
            cur = tx.get_text()
            if cur != getattr(tx, "_two_lines", None):   # (a text set since: that one)
                tx._one_line = cur
            one = tx._one_line
            tx._two_lines = None
            tx.set_text(one)
            words = one.split()
            if "\n" in one or len(words) < 2 or tx.get_window_extent(r).width <= 0.94 * width_px:
                continue
            k = min(range(1, len(words)), key=lambda i: abs(len(" ".join(words[:i])) - len(" ".join(words[i:]))))
            tx._two_lines = " ".join(words[:k]) + "\n" + " ".join(words[k:])
            tx.set_text(tx._two_lines)

    def draw(self):
        self.place()
        self.canvas.draw_idle()

    def _resized(self, e=None):
        """New size of the plot: the margins, and an automatic y range made
        again, since the room kept above the tallest peak for its label is
        in points (a row drawn before its final size had labels across the
        top of the frame)."""
        self.place()
        a = self._auto
        if a is None or self.ylock is not None or self.mode == "map":
            return
        try:
            if np.allclose(self.ax.get_ylim(), a[1]):
                self.autoscale_y(a[0])
        except Exception:
            pass

    # --------------------------------------------------------------- mouse
    TOOL_CURSORS = {"zoom": wx.CURSOR_MAGNIFIER, "pan": wx.CURSOR_SIZING, "background": wx.CURSOR_SIZEWE,
                    "drag_peak": wx.CURSOR_SIZEWE, "click_peak": wx.CURSOR_HAND, "split": wx.CURSOR_CROSS,
                    "delete": wx.CURSOR_HAND, "xic": wx.CURSOR_HAND, "measure": wx.CURSOR_CROSS,
                    "calibrate": wx.CURSOR_SIZEWE, "area": wx.CURSOR_SIZEWE,
                    "label": wx.CURSOR_HAND}
    SPAN_TOOLS = ("background", "drag_peak", "calibrate", "area")
    CLICK_TOOLS = ("click_peak", "split", "delete", "xic", "measure", "label")

    def _tool(self):
        t = self.tool_getter() if self.tool_getter else "select"
        return t if (t in ("zoom", "pan") or t in self.tools_supported) else "select"

    def _clamped(self, e):
        """Data coordinates of the mouse, held at the edge of the plot when the
        pointer is outside it (so dragging past the axes keeps working)."""
        if e.x is None or e.y is None:
            return None, None
        bb = self.ax.bbox
        px = min(max(e.x, bb.x0), bb.x1)
        py = min(max(e.y, bb.y0), bb.y1)
        try:
            x, y = self.ax.transData.inverted().transform((px, py))
        except Exception:
            return None, None
        return float(x), float(y)

    def _clear_drag_art(self):
        if self._drag_art is not None:
            art, self._drag_art = self._drag_art, None
            if art in self.live:
                self.remove_live([art])
            else:
                try:
                    art.remove()
                except Exception:
                    pass

    def _down(self, e):
        import shortcuts
        shortcuts.plot_clicked(self)  # the plot takes the keyboard focus (tool letters, Ctrl+C, Ctrl+E)
        if e.inaxes is not self.ax or e.xdata is None:
            return
        if e.button == 1 and e.dblclick:
            self.reset_view()
            return
        if e.button == 3:
            if self.on_menu:
                items = self.on_menu(e.xdata, e.ydata) or []
                if items:
                    wx.CallAfter(self._popup, items)
            return
        if e.button != 1:
            return
        tool = self._tool()
        ctrl = bool(e.key and ("control" in e.key or "ctrl" in e.key)) or wx.GetKeyState(wx.WXK_CONTROL)
        if ctrl and self.mode != "map":
            self._press = ("zoom", e.xdata, e.ydata, e.x, e.y)
            return
        if tool == "select":
            shift = bool(e.key and "shift" in e.key) or wx.GetKeyState(wx.WXK_SHIFT)
            self._press = ("select", e.xdata, e.x, shift)
        elif tool == "zoom":
            self._press = ("zoom", e.xdata, e.ydata, e.x, e.y)
        elif tool == "pan":
            self._press = ("pan", e.x, e.y, self.ax.get_xlim(), self.ax.get_ylim())
        elif tool in self.SPAN_TOOLS:
            self._press = ("span", tool, e.xdata, e.x)
        else:
            self._press = ("click", tool, e.xdata, e.ydata, e.x, e.y)

    def _hotkey(self, e):
        """Ctrl+C copies the plot as an image, Ctrl+S saves it (after a click
        on the plot)."""
        k = (e.key or "").lower()
        if k in ("ctrl+c", "control+c"):
            wx.CallAfter(self.copy_image)
        elif k in ("ctrl+s", "control+s"):
            wx.CallAfter(self.save_image)

    def _popup(self, items):
        items = [it for it in items if not (it[0] and it[0].startswith("Fit every tile"))]
        top = wx.GetTopLevelParent(self)
        if hasattr(top, "make_report"):
            items = list(items) + [(None, None), ("Create a report\u2026 (Ctrl+R)", lambda: top.make_report())]
        import undo
        items = list(items) + [(None, None)] + undo.menu_items(top)  # Undo: <step> (Ctrl+Z), Redo (Ctrl+Y)
        items = list(items) + [(None, None), ("Copy image (Ctrl+C)", self.copy_image),
                               ("Save image as\u2026 (Ctrl+S)", self.save_image)]
        if self.pair is not None:  # (the entries above: this plot only)
            items += [("Copy image of both plots", self.pair.copy_image),
                      ("Save image of both plots as\u2026", self.pair.save_image)]
        items += [("Image size for copy and save\u2026", lambda: ask_export_settings(self))]
        if _stack_of(self) is not None:
            items += tile_menu(self)
        menu = build_menu(self, items)
        popup_menu(self.canvas, menu)
        menu.Destroy()

    def _move(self, e):
        p = self._press
        if p is None:
            tool = self._tool()
            inside = e.inaxes is self.ax
            self.canvas.SetCursor(wx.Cursor(self.TOOL_CURSORS.get(tool, wx.CURSOR_ARROW) if inside
                                            else wx.CURSOR_ARROW))
            if self.readout_fmt is not None:
                ok = inside and e.xdata is not None and self.full is not None
                self._set_readout(self.readout_fmt(e.xdata, e.ydata) if ok else "")
            if self.tip_fmt is not None:
                try:
                    tip = self.tip_fmt(e.xdata, e.ydata) if inside and e.xdata is not None else ""
                except Exception:
                    tip = ""
                self._set_plot_tip(tip or "")
            return
        kind = p[0]
        if kind == "pan":
            bb = self.ax.bbox
            (x0, x1), (y0, y1) = p[3], p[4]
            dx = (e.x - p[1]) / max(bb.width, 1) * (x1 - x0)
            dy = (e.y - p[2]) / max(bb.height, 1) * (y1 - y0)
            self.ax.set_xlim(x0 - dx, x1 - dx)
            self.ax.set_ylim(y0 - dy, y1 - dy)
            self.ylock = self.ax.get_ylim()
            self.canvas.draw_idle()
            return
        if kind == "click":
            return
        ex, ey = self._clamped(e)
        if ex is None:
            return
        self._clear_drag_art()
        if kind == "zoom":
            from matplotlib.patches import Rectangle
            x0, y0 = p[1], p[2]
            self._drag_art = self.ax.add_patch(Rectangle((min(x0, ex), min(y0, ey)), abs(ex - x0),
                                                         abs(ey - y0), fill=True, fc=(0.16, 0.38, 0.77, 0.08),
                                                         ec=C["accent"], lw=0.9, ls="--"))
        else:
            x0 = p[1] if kind == "select" else p[2]
            if kind == "select":
                col, alpha = (BG_COL if p[3] else AVG_COL), 0.18
            else:
                col, alpha = {"background": (BG_COL, 0.2), "calibrate": ("#0BA064", 0.18)}.get(p[1], ("#D2456F", 0.16))
            self._drag_art = self.ax.axvspan(min(x0, ex), max(x0, ex), color=col, alpha=alpha, lw=0)
        if self._blit_ok:
            # only the rectangle is drawn again, over the stored image of the plot
            self.add_live(self._drag_art)
            self.blit_live()
        else:
            self.canvas.draw_idle()

    def _up(self, e):
        p = self._press
        if p is None:
            return
        self._press = None
        self._clear_drag_art()
        kind = p[0]
        inside = e.xdata is not None and e.inaxes is self.ax
        if kind in ("zoom", "select", "span") and not inside:
            # released outside the plot: use the point held at its edge
            cx, cy = self._clamped(e)
            if cx is not None:
                e.xdata, e.ydata, inside = cx, cy, True
        if kind == "pan":
            self.canvas.draw_idle()
            self._view_changed()
            return
        if kind == "zoom":
            px, py = p[3], p[4]
            if inside and abs(e.x - px) > 5 and abs(e.y - py) > 5:
                self.ax.set_xlim(min(p[1], e.xdata), max(p[1], e.xdata))
                self.ax.set_ylim(min(p[2], e.ydata), max(p[2], e.ydata))
                self.ylock = self.ax.get_ylim()
                self._view_changed()
            self.canvas.draw_idle()
            return
        if kind == "click":
            self.canvas.draw_idle()
            if inside and abs(e.x - p[4]) < 5 and abs(e.y - p[5]) < 5 and self.on_tool:
                self.on_tool(p[1], "click", e.xdata, e.ydata)
            return
        if kind == "span":
            tool, x0, px0 = p[1], p[2], p[3]
            self.canvas.draw_idle()
            if not inside or abs(e.x - px0) < 4:
                return
            a, b = min(x0, e.xdata), max(x0, e.xdata)
            if tool == "background":
                if self.on_range:
                    self.on_range(a, b, True)
            elif self.on_tool:
                self.on_tool(tool, "drag", a, b)
            return
        # select
        x0, px0, shift = p[1], p[2], p[3]
        x1 = e.xdata if inside else None
        if x1 is None or abs(e.x - px0) < 4:
            self.canvas.draw_idle()
            if abs(e.x - px0) < 4 and self.on_pick and x1 is not None:
                self.on_pick(x1, e.ydata)
            return
        a, b = min(x0, x1), max(x0, x1)
        if self.mode == "zoom" and not shift:
            self.ax.set_xlim(a, b)
            self.ylock = None
            self.autoscale_y()
            self.canvas.draw_idle()
            self._view_changed()
        elif self.on_range:
            self.on_range(a, b, shift)
        else:
            self.canvas.draw_idle()

    def _scroll(self, e):
        if e.inaxes is not self.ax or e.xdata is None:
            return
        f = 0.8 if e.button == "up" else 1.25
        x0, x1 = self.ax.get_xlim()
        a = e.xdata - (e.xdata - x0) * f
        b = e.xdata + (x1 - e.xdata) * f
        if self.full:
            a, b = max(a, self.full[0]), min(b, self.full[1])
        if b > a:
            self.ax.set_xlim(a, b)
            self.ylock = None
            self.autoscale_y()
            self.canvas.draw_idle()
            self._view_changed()

    def reset_view(self):
        self.ylock = None
        if self.xlock or self.full:
            self.ax.set_xlim(*(self.xlock or self.full))
        self.autoscale_y()
        self.canvas.draw_idle()
        self._view_changed()

    def autoscale_y(self, pad=0.08):
        """y range from the visible part of the lines (unless the y range was
        set with the zoom box or by panning)."""
        base_pad = pad
        if self.ylock is not None:
            self.ax.set_ylim(*self.ylock)
            return
        x0, x1 = self.ax.get_xlim()
        lo, hi = None, None
        sources = [(ln.full_data() if hasattr(ln, "full_data") else (ln.get_xdata(), ln.get_ydata()))
                   for ln in self.ax.get_lines()
                   if ln.get_visible() and not getattr(ln, "_no_scale", False)] + list(self.extra_scale)
        for x, y in sources:
            if len(x) == 0:
                continue
            x = np.asarray(x, float)
            y = np.asarray(y, float)
            m = (x >= min(x0, x1)) & (x <= max(x0, x1))
            if not np.any(m):
                continue
            ym = y[m]
            ym = ym[np.isfinite(ym)]
            if not len(ym):
                continue
            lo = ym.min() if lo is None else min(lo, ym.min())
            hi = ym.max() if hi is None else max(hi, ym.max())
        if lo is None:
            fx = self.y_fixed
            if fx is not None and fx[0] is not None and fx[1] is not None:  # nothing in view: the range set
                self.ax.set_ylim(*fx)
            return
        if hi <= lo:
            hi = lo + 1
        # room above the tallest peak for its label, also on short tiles
        try:
            h_pt = self.ax.get_position().height * self.fig.get_size_inches()[1] * 72.0
            f = min(0.4, self.label_room_pt / max(h_pt, 1.0))
            pad = max(pad, f / (1.0 - f))
        except Exception:
            pass
        ymin = 0.0 if lo >= 0 else lo - 0.03 * (hi - lo)
        self.ax.set_ylim(ymin, hi + pad * (hi - ymin))
        self._auto = (base_pad, tuple(self.ax.get_ylim()))
        fx = self.y_fixed
        if fx is not None and (fx[0] is not None or fx[1] is not None):  # set by the owner: kept (also in images)
            a, b = self.ax.get_ylim()
            lo, hi = (fx[0] if fx[0] is not None else a), (fx[1] if fx[1] is not None else b)
            if hi <= lo:  # one end set beyond the other end of the automatic range (the axis turned upside down)
                if fx[0] is not None:
                    hi = lo + (b - a)
                else:
                    lo = hi - (b - a)
            self.ax.set_ylim(lo, hi)
            self._auto = None


class PairCard(wx.Panel):
    """One tile with plots one below the other that share the x axis (the positive and the negative ions of a
    file): each plot is a PlotCard of its own (mouse, tools, menus, images of that plot alone), drawn without its
    own card and title row (its name and the power of ten of its y axis are its y title, the text of the title
    row a note inside the plot); the x numbers and the unit only under the lowest. Zooming one plot zooms the
    others."""

    own_refresh = True   # (paints itself again when resized: place_panes)

    def __init__(self, parent, title):
        wx.Panel.__init__(self, parent, style=wx.BORDER_NONE)
        self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
        self.SetBackgroundColour(C["bg"])
        self.title = title  # (Closed tiles, file names of images)
        self.cards = []
        self.Bind(wx.EVT_PAINT, self._paint)
        self.Bind(wx.EVT_SIZE, lambda e: (self.relayout(), e.Skip()))
        self.Bind(wx.EVT_MOUSEWHEEL, lambda e: forward_wheel(self, e))

    def set_cards(self, cards):
        self.cards = list(cards)
        for i, c in enumerate(self.cards):
            c.pair = self
            c.hide_x = i < len(self.cards) - 1
            c.set_shared()
        self.relayout()

    def relayout(self):
        w, h = self.GetClientSize()
        n = len(self.cards)
        if not n or w < 10 or h < 10:
            return
        sx, top, bot = self.FromDIP(4), self.FromDIP(3), self.FromDIP(7)
        inner = max(n, h - top - bot)
        # the x numbers under the lowest plot (its bottom margin minus that of the others): frames of one height
        try:
            dpi = float(self.cards[-1].fig.dpi)
        except Exception:
            dpi = float(self.FromDIP(100))
        b0, b1 = PlotCard.SHARED_B
        extra = min(int(round((b1 - b0) * dpi)), inner // 3)
        each = max(1, (inner - extra) // n)
        rects, y = [], top
        for i, c in enumerate(self.cards):
            ch = each if i < n - 1 else inner - each * (n - 1)
            rects.append((c, (sx, y, max(1, w - 2 * sx), max(1, ch))))
            y += ch
        place_panes(self, rects)

    def _paint(self, e):
        paint_card(self, self.FromDIP(12) + 4, card_frame, lambda gc, w, font=None: None, ("tile",))

    def sync_x(self, src):
        """The other plots at the x range of src (and their y ranges made again)."""
        xl = src.ax.get_xlim()
        for c in self.cards:
            if c is src or not c.full:
                continue
            c._syncing = True
            try:
                c.ax.set_xlim(*xl)
                c.ylock = None
                c.autoscale_y()
                if c.on_view:
                    c.on_view()
                c.draw()
            except RuntimeError:
                pass
            finally:
                c._syncing = False

    # ------------------------------------------------------ image of the whole tile
    def _image(self, fmt):
        """The plots one below the other as one image (PNG, TIFF or SVG), each at the size of the image
        settings (the upper ones without x numbers), their frames lined up."""
        import io
        cards = [c for c in self.cards if c.full is not None] or list(self.cards)
        meas = [c._render("png", measure=True, shared=True) for c in cards]
        align = None
        if all(m is not None for m in meas):
            dx = [m[3] for m in meas if m[3] is not None]
            align = (max(m[0] for m in meas), max(m[1] for m in meas), max(dx) if dx else None)
        cfg = export_settings()
        h_in = float(cfg["height_cm"]) / 2.54
        parts = []
        for i, c in enumerate(cards):
            last = i == len(cards) - 1
            hi = None if (last or meas[-1] is None) else max(0.6, h_in - meas[-1][2] + 0.08)
            parts.append(c._render("svg" if fmt == "svg" else "png", align=align, shared=True,
                                   height_in=hi).getvalue())
        if fmt == "svg":
            return _stack_svg(parts)
        from PIL import Image
        ims = [Image.open(io.BytesIO(b)).convert("RGB") for b in parts]
        W = max(im.width for im in ims)
        out = Image.new("RGB", (W, sum(im.height for im in ims)), "white")
        y = 0
        for im in ims:
            out.paste(im, (0, y))
            y += im.height
        buf = io.BytesIO()
        dpi = int(cfg["dpi"])
        out.save(buf, format="TIFF" if fmt in ("tif", "tiff") else "PNG", dpi=(dpi, dpi))
        return buf.getvalue()

    def _name(self):
        name = None
        try:
            name = self.cards[0].image_name() if self.cards and self.cards[0].image_name else None
        except Exception:
            name = None
        if name and "_pos" in name:
            return name.replace("_pos", "_pos_neg", 1)
        frame = wx.GetTopLevelParent(self)
        base = frame.base_name() if hasattr(frame, "base_name") else "plot"
        return (name or base) + "_" + "".join(ch if ch.isalnum() else "_" for ch in self.title.lower()).strip("_")

    def copy_image(self):
        try:
            data = self._image("png")
            how = copy_png_to_clipboard(self, data)
            msg = "Image of both plots copied: paste it with Ctrl+V"
            _status(self, msg)
            toast(self, msg)
            _log("copy image of a tile: %d bytes, %s" % (len(data), how))
        except Exception as ex:
            _report_error(self, "The image could not be copied", ex)

    def save_image(self):
        try:
            frame = wx.GetTopLevelParent(self)
            name = safe_file_name(self._name())
            wild = "PNG image (*.png)|*.png|TIFF image (*.tif)|*.tif|SVG, vector (*.svg)|*.svg"
            if hasattr(frame, "ask_save"):
                path = frame.ask_save("Save image of both plots", name + ".png", wild)
            else:
                path = ask_save_file(frame, "Save image of both plots", name + ".png", wild)
            if not path:
                return
            fmt = os.path.splitext(path)[1].lower().lstrip(".")
            if fmt not in ("png", "tif", "tiff", "svg"):
                fmt = "png"
                path += ".png"
            data = self._image(fmt)
            with open(path, "wb") as fh:
                fh.write(data)
            _status(self, "Saved " + path)
            toast(self, "Saved %s\nin %s" % (os.path.basename(path), os.path.dirname(path)))
        except Exception as ex:
            _report_error(self, "The image could not be saved", ex)


def _stack_svg(parts):
    """SVG images (bytes) one below the other in one SVG."""
    import re
    inner, W, H = [], 0.0, 0.0
    for b in parts:
        t = b.decode("utf-8")
        t = re.sub(r"<\?xml[^>]*\?>", "", t)
        t = re.sub(r"<!DOCTYPE[^>]*>", "", t, flags=re.S)
        m = re.search(r"<svg\b[^>]*>", t)
        head = m.group(0)
        w = float(re.search(r'width="([\d.]+)', head).group(1))
        h = float(re.search(r'height="([\d.]+)', head).group(1))
        inner.append(t[:m.start()] + head.replace("<svg", '<svg x="0" y="%g"' % H, 1) + t[m.end():])
        W, H = max(W, w), H + h
    out = ('<?xml version="1.0" encoding="utf-8" standalone="no"?>\n'
           '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" version="1.1" '
           'width="%gpt" height="%gpt" viewBox="0 0 %g %g">\n' % (W, H, W, H)) + "\n".join(inner) + "\n</svg>\n"
    return out.encode("utf-8")


_MPL_READY = []


def _mpl_setup():
    """Publication fonts for the plots: Arial (as journals ask), embedded as
    real text in PDF and SVG files; crisper text on screen."""
    if _MPL_READY:
        return
    _MPL_READY.append(True)
    try:
        import matplotlib
        from matplotlib import font_manager
        names = {f.name for f in font_manager.fontManager.ttflist}
        fam = [n for n in ("Arial", "Helvetica", "Liberation Sans", "DejaVu Sans") if n in names]
        if fam:
            matplotlib.rcParams["font.family"] = "sans-serif"
            matplotlib.rcParams["font.sans-serif"] = fam + list(matplotlib.rcParams["font.sans-serif"])
        matplotlib.rcParams["pdf.fonttype"] = 42
        matplotlib.rcParams["ps.fonttype"] = 42
        matplotlib.rcParams["svg.fonttype"] = "none"
        matplotlib.rcParams["text.hinting_factor"] = 1
        matplotlib.rcParams["text.hinting"] = "default"
        matplotlib.rcParams["mathtext.default"] = "regular"
    except Exception as ex:
        _log("plot fonts: %s" % ex)


def _status(win, text):
    try:
        wx.GetTopLevelParent(win).SetStatusText(text, 0)
    except Exception:
        pass


EXPORT_DEFAULTS = {"mode": "journal", "width_cm": 8.5, "height_cm": 6.0, "text": 9.0, "dpi": 300}
_SET = {"mtime": None, "data": {}, "gen": None}


_ENV = {"cls": None}


def _envelope_class():
    """Line2D of a dense trace drawn from its lowest and highest point in
    each pixel column of the current view (at most two points per column):
    it looks the same as the full trace but draws many times faster, while
    zooming, panning and dragging. Saved images get the envelope at their
    own resolution; the data themselves are never changed."""
    if _ENV["cls"] is not None:
        return _ENV["cls"]
    from matplotlib.lines import Line2D
    from matplotlib.path import Path

    class EnvelopeLine(Line2D):
        def __init__(self, x, y, fill=None, fill_base=0.0, **kw):
            x, y = np.asarray(x, float), np.asarray(y, float)
            self._fill_base = float(fill_base)  # the area reaches down to this y (a stacked trace: its offset)
            bad = ~np.isfinite(x)
            if bad.any():
                # gaps (a fit is NaN between its m/z windows): each part on its own;
                # searchsorted on an array with NaN in it found wrong ends, and a
                # 45000 point fit was drawn from 4 points (or none) when not zoomed
                cut = np.flatnonzero(bad)
                starts, ends = np.r_[0, cut + 1], np.r_[cut, len(x)]
                self._segs = [(s, e) for s, e in zip(starts, ends) if e > s]
                unsorted = [(s, e) for s, e in self._segs if e - s > 1 and np.any(np.diff(x[s:e]) < 0)]
                if unsorted:
                    x, y = x.copy(), y.copy()  # the caller's arrays stay as they are
                    for s, e in unsorted:
                        o = s + np.argsort(x[s:e], kind="stable")
                        x[s:e], y[s:e] = x[o], y[o]
            else:
                self._segs = None
                if len(x) > 1 and np.any(np.diff(x) < 0):
                    o = np.argsort(x, kind="stable")
                    x, y = x[o], y[o]
            self._fx, self._fy = x, y
            self._fill = fill  # (colour, alpha): area under the trace, or None
            Line2D.__init__(self, x, y, **kw)

        def full_data(self):
            return self._fx, self._fy

        def _view_data(self):
            x, y = self._fx, self._fy
            ax = self.axes
            if ax is None or len(x) < 3:
                return x, y
            a, b = sorted(ax.get_xlim())
            npx = max(50, int(ax.bbox.width))
            if self._segs is None:
                return self._reduce(x, y, a, b, 2 * npx)
            outx, outy = [], []
            for s, e in self._segs:
                xs = x[s:e]
                if xs[-1] < a or xs[0] > b:
                    continue
                share = (min(xs[-1], b) - max(xs[0], a)) / ((b - a) or 1.0)
                rx, ry = self._reduce(xs, y[s:e], a, b, max(8, int(2 * npx * share) + 2))
                if outx:
                    outx.append([np.nan])
                    outy.append([np.nan])
                outx.append(rx)
                outy.append(ry)
            if not outx:
                return x[:0], y[:0]
            return np.concatenate(outx), np.concatenate(outy)

        @staticmethod
        def _reduce(x, y, a, b, nb):
            """Points of a rising, gapless x for the view a to b in nb half
            pixel columns."""
            i0 = max(0, int(np.searchsorted(x, a)) - 1)
            i1 = min(len(x), int(np.searchsorted(x, b)) + 1)
            xs, ys = x[i0:i1], y[i0:i1]
            if len(xs) <= 3 * nb:
                return xs, ys
            # lowest and highest point of each column, at their own x and in their order, in O(n):
            # x rises, so each column is a run of points; the same points as the lexsort version
            n = len(xs)
            starts = np.unique(np.searchsorted(xs, np.linspace(xs[0], xs[-1], nb + 1)[:-1], side="left"))
            starts = starts[starts < n]
            run = np.repeat(np.arange(len(starts)), np.diff(np.r_[starts, n]))
            i_lo = np.flatnonzero(ys == np.minimum.reduceat(ys, starts)[run])
            i_lo = i_lo[np.r_[True, run[i_lo][1:] != run[i_lo][:-1]]]
            i_hi = np.flatnonzero(ys == np.maximum.reduceat(ys, starts)[run])
            i_hi = i_hi[np.r_[run[i_hi][1:] != run[i_hi][:-1], True]]
            keep = np.zeros(n, bool)
            keep[i_lo] = True
            keep[i_hi] = True
            keep[0] = keep[-1] = True
            pick = np.flatnonzero(keep)
            return xs[pick], ys[pick]

        def draw(self, renderer):
            if not self.get_visible():
                return
            try:
                # pixel images (screen, PNG, TIFF) get the envelope; vector
                # files (PDF, SVG) keep every point, as they may be zoomed
                raster = type(renderer).__name__ == "RendererAgg"
                xs, ys = self._view_data() if raster else (self._fx, self._fy)
                Line2D.set_data(self, xs, ys)
                if self._fill and len(xs) > 1:
                    colour, alpha = self._fill
                    # closed polygon: matplotlib ignores the last vertex of a closed Path (it
                    # becomes the CLOSEPOLY code), so the corner (xs[0], 0) is given twice;
                    # without it the area was closed by a diagonal from the bottom right to
                    # the first point (a sloped band when the view starts on a peak)
                    ok = np.isfinite(xs) & np.isfinite(ys)
                    runs = np.split(np.arange(len(xs)), np.flatnonzero(~ok)) if not ok.all() else [np.arange(len(xs))]
                    gc = renderer.new_gc()
                    try:
                        from matplotlib.colors import to_rgba
                        gc.set_alpha(alpha)
                        gc.set_linewidth(0)
                        gc.set_clip_rectangle(self.axes.bbox)
                        rgb = to_rgba(colour)
                        for idx in runs:  # each part between gaps has its own area
                            idx = idx[ok[idx]]
                            if len(idx) < 2:
                                continue
                            px, py = xs[idx], ys[idx]
                            fb = self._fill_base
                            verts = np.concatenate([np.column_stack([px, py]),
                                                    [[px[-1], fb], [px[0], fb], [px[0], fb]]])
                            renderer.draw_path(gc, Path(verts, closed=True), self.get_transform(),
                                               (rgb[0], rgb[1], rgb[2], alpha))
                    finally:
                        gc.restore()
            except Exception as ex:
                _log("envelope: %s" % ex)
            Line2D.draw(self, renderer)

    _ENV["cls"] = EnvelopeLine
    return EnvelopeLine


def plot_trace(ax, x, y, fill=None, fill_base=0.0, **kw):
    """A dense trace (spectrum, mass spectrum) as an envelope line; fill =
    (colour, alpha) shades the area under it, down to fill_base."""
    ln = _envelope_class()(x, y, fill=fill, fill_base=fill_base, **kw)
    ax.add_line(ln)
    ax.autoscale_view()
    return ln


def settings():
    """The settings file, read again only when it changed (layout code asks
    for settings on every mouse move)."""
    try:
        path = T._ST.get("settings_path")
        m = os.path.getmtime(path) if path and os.path.isfile(path) else None
    except Exception:
        m = None
    gen = T._ST.get("gen", 0)  # (a save of this program within one tick of the file time)
    if m != _SET["mtime"] or m is None or gen != _SET["gen"]:
        _SET["data"] = T._load()
        _SET["mtime"] = m
        _SET["gen"] = gen
    return _SET["data"]


def _report_error(win, what, ex):
    """Message with the reason, and the full traceback in the log file."""
    import traceback
    tb = traceback.format_exc()
    try:
        print("%s:\n%s" % (what, tb))
        sys.stdout.flush()
    except Exception:
        pass
    try:
        wx.MessageBox("%s:\n%s\n\nDetails are in the log file." % (what, ex or ex.__class__.__name__),
                      APP, wx.ICON_ERROR, wx.GetTopLevelParent(win) if win else None)
    except Exception:
        pass


def safe_file_name(name):
    bad = '<>:"/\\|?*'
    out = "".join("_" if (c in bad or ord(c) < 32) else c for c in str(name)).strip(" .")
    return out[:150] or "plot"


def ask_save_file(parent, title, name, wildcard):
    dlg = wx.FileDialog(parent, title, defaultFile=name, wildcard=wildcard,
                        style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
    try:
        if dlg.ShowModal() != wx.ID_OK:
            return None
        path = dlg.GetPath()
        ext = wildcard.split("|")[1 + 2 * dlg.GetFilterIndex()].replace("*", "").split(";")[0]
        if ext and not path.lower().endswith(ext.lower()):
            path += ext
        return path
    finally:
        dlg.Destroy()


def _win_clipboard_png(hwnd, png):
    """Put a PNG (registered format "PNG", read by Word, PowerPoint, Outlook)
    and the same picture as a device independent bitmap (CF_DIB, read by
    everything else) on the Windows clipboard with the Win32 API."""
    import io
    import ctypes
    from ctypes import wintypes
    from PIL import Image
    u32 = ctypes.WinDLL("user32", use_last_error=True)
    k32 = ctypes.WinDLL("kernel32", use_last_error=True)
    u32.OpenClipboard.argtypes = [wintypes.HWND]
    u32.OpenClipboard.restype = wintypes.BOOL
    u32.EmptyClipboard.argtypes = []
    u32.EmptyClipboard.restype = wintypes.BOOL
    u32.CloseClipboard.argtypes = []
    u32.CloseClipboard.restype = wintypes.BOOL
    u32.SetClipboardData.argtypes = [wintypes.UINT, ctypes.c_void_p]
    u32.SetClipboardData.restype = ctypes.c_void_p
    u32.RegisterClipboardFormatW.argtypes = [wintypes.LPCWSTR]
    u32.RegisterClipboardFormatW.restype = wintypes.UINT
    u32.IsClipboardFormatAvailable.argtypes = [wintypes.UINT]
    u32.IsClipboardFormatAvailable.restype = wintypes.BOOL
    k32.GlobalAlloc.argtypes = [wintypes.UINT, ctypes.c_size_t]
    k32.GlobalAlloc.restype = ctypes.c_void_p
    k32.GlobalLock.argtypes = [ctypes.c_void_p]
    k32.GlobalLock.restype = ctypes.c_void_p
    k32.GlobalUnlock.argtypes = [ctypes.c_void_p]
    k32.GlobalUnlock.restype = wintypes.BOOL
    k32.GlobalFree.argtypes = [ctypes.c_void_p]
    k32.GlobalFree.restype = ctypes.c_void_p
    GMEM_MOVEABLE, CF_DIB = 0x0002, 8

    im = Image.open(io.BytesIO(png))
    dpi = im.info.get("dpi") or (300, 300)
    bmp = io.BytesIO()
    im.convert("RGB").save(bmp, "BMP", dpi=(round(dpi[0]), round(dpi[1])))
    dib = bmp.getvalue()[14:]  # without the 14 byte file header
    fmt_png = u32.RegisterClipboardFormatW("PNG")

    def glob(data):
        h = k32.GlobalAlloc(GMEM_MOVEABLE, len(data))
        if not h:
            raise OSError("no memory for the clipboard")
        p = k32.GlobalLock(h)
        if not p:
            k32.GlobalFree(h)
            raise OSError("no memory for the clipboard")
        ctypes.memmove(p, data, len(data))
        k32.GlobalUnlock(h)
        return h

    opened = False
    for _ in range(20):  # another program may hold the clipboard for a moment
        if u32.OpenClipboard(hwnd):
            opened = True
            break
        time.sleep(0.05)
    if not opened:
        raise OSError("the clipboard is in use by another program; try again")
    try:
        if not u32.EmptyClipboard():
            raise OSError("the clipboard could not be emptied (error %d)" % ctypes.get_last_error())
        for fmt, data in ((fmt_png, png), (CF_DIB, dib)):
            h = glob(data)
            if not u32.SetClipboardData(fmt, h):
                err = ctypes.get_last_error()
                k32.GlobalFree(h)
                raise OSError("the clipboard refused the image (error %d)" % err)
    finally:
        u32.CloseClipboard()
    if not u32.IsClipboardFormatAvailable(fmt_png):
        raise OSError("the image is not on the clipboard after copying")


def _wx_clipboard_png(png):
    import io
    img = wx.Image(io.BytesIO(png), wx.BITMAP_TYPE_PNG)
    comp = wx.DataObjectComposite()
    obj = wx.CustomDataObject(wx.DataFormat("PNG"))
    obj.SetData(png)
    comp.Add(obj, True)
    comp.Add(wx.BitmapDataObject(wx.Bitmap(img)))
    if not wx.TheClipboard.Open():
        raise OSError("the clipboard is in use by another program; try again")
    try:
        if not wx.TheClipboard.SetData(comp):
            raise OSError("the clipboard refused the image")
        wx.TheClipboard.Flush()  # stays on the clipboard after the program is closed
    finally:
        wx.TheClipboard.Close()


def copy_png_to_clipboard(win, png):
    """Returns how it was copied; raises with the reason if it was not."""
    if os.name == "nt":
        try:
            top = wx.GetTopLevelParent(win) if win else None
            hwnd = top.GetHandle() if top else None
            _win_clipboard_png(hwnd, png)
            return "Win32 clipboard (PNG + DIB)"
        except Exception as ex:
            _log("Win32 clipboard failed (%s); trying wx" % ex)
    _wx_clipboard_png(png)
    return "wx clipboard"


class _Toast(wx.PopupWindow):
    def __init__(self, parent, text):
        wx.PopupWindow.__init__(self, parent, wx.BORDER_SIMPLE)
        self.SetBackgroundColour(wx.Colour("#1F2937"))
        st = wx.StaticText(self, label=text)
        st.SetForegroundColour(wx.Colour("#FFFFFF"))
        st.SetFont(ui_font(9, 500))
        s = wx.BoxSizer(wx.VERTICAL)
        s.Add(st, 0, wx.ALL, self.FromDIP(10))
        self.SetSizerAndFit(s)


def toast(win, text, ms=2600):
    """Short confirmation at the bottom of the window, gone by itself."""
    try:
        top = wx.GetTopLevelParent(win)
        t = _Toast(top, text)
        r = top.GetScreenRect()
        w, h = t.GetSize()
        t.SetPosition(wx.Point(r.x + (r.width - w) // 2, r.y + r.height - h - top.FromDIP(56)))
        t.Show()

        def gone(t=t):
            try:
                if t:
                    t.Destroy()
            except Exception:
                pass
        wx.CallLater(ms, gone)
    except Exception:
        pass


def fit_axes(fig, pad_pt=3.0):
    """Margins that hold every label, title and tick of the plot, and no more:
    the axes (all of them, keeping their arrangement) are moved and scaled
    so that everything drawn around them just fits inside the figure.
    tight_layout cannot do this for axes placed by hand, which is why saved
    images lost the x axis title or a peak label at the edge."""
    from matplotlib.transforms import Bbox
    axes = [a for a in fig.axes if a.get_visible()]
    if not axes:
        return
    for _ in range(3):
        try:
            # measure only (no redraw of the window, which would flash the
            # export layout on screen)
            r = fig.canvas.get_renderer()
            tb = Bbox.union([a.get_tightbbox(r) for a in axes])
            ab = Bbox.union([a.get_window_extent(r) for a in axes])
        except Exception:
            return
        W, H = fig.bbox.width, fig.bbox.height
        pad = pad_pt * fig.dpi / 72.0
        l, rgt = ab.x0 - tb.x0 + pad, tb.x1 - ab.x1 + pad
        b, top = ab.y0 - tb.y0 + pad, tb.y1 - ab.y1 + pad
        tx0, tx1, ty0, ty1 = l, W - rgt, b, H - top
        if tx1 - tx0 < 0.25 * W or ty1 - ty0 < 0.25 * H or ab.width <= 0 or ab.height <= 0:
            return
        if abs(tx0 - ab.x0) < 0.5 and abs(tx1 - ab.x1) < 0.5 and abs(ty0 - ab.y0) < 0.5 and abs(ty1 - ab.y1) < 0.5:
            return
        sx, sy = (tx1 - tx0) / ab.width, (ty1 - ty0) / ab.height
        for a in axes:
            e = a.get_window_extent(r)
            x0 = tx0 + (e.x0 - ab.x0) * sx
            y0 = ty0 + (e.y0 - ab.y0) * sy
            a.set_position([x0 / W, y0 / H, e.width * sx / W, e.height * sy / H])


def _hook_save():
    """Forget the cached settings whenever anything is saved (the file time
    alone can miss two writes within the same clock tick)."""
    orig = getattr(T, "_save", None)
    if orig is None or getattr(orig, "_cache_hook", False):
        return

    def _save(update, _orig=orig):
        try:
            return _orig(update)
        finally:
            _SET["mtime"] = None
            _SET["data"] = {}
    _save._cache_hook = True
    T._save = _save


_hook_save()


def export_settings():
    d = dict(EXPORT_DEFAULTS)
    try:
        saved = settings().get("figure_export")
        if isinstance(saved, dict):
            d.update({k: v for k, v in saved.items() if k in d})
    except Exception:
        pass
    return d


def ask_export_settings(win):
    """Size of copied and saved images: journal size (default 8.5 x 6 cm,
    single column, text 9 pt: readable when pasted into Word) or as on
    screen."""
    d = export_settings()
    sizes = [8.5, 12.0, 17.5]
    v = ask_form(win, "Image size", [
        dict(key="mode", label="Size", kind="choice", choices=["Journal size (below)", "As on screen"],
             value=0 if d["mode"] == "journal" else 1, width=190),
        dict(key="preset", label="Width preset", kind="choice",
             choices=["Single column, 8.5 cm", "1.5 columns, 12 cm", "Two columns, 17.5 cm", "Custom (below)"],
             value=sizes.index(d["width_cm"]) if d["width_cm"] in sizes else 3, width=190),
        dict(key="w", label="Width", value="%g" % d["width_cm"], unit="cm"),
        dict(key="h", label="Height", value="%g" % d["height_cm"], unit="cm"),
        dict(key="text", label="Text size", kind="choice", choices=["7 pt", "8 pt", "9 pt", "10 pt", "11 pt", "12 pt"],
             value=int(round(d["text"])) - 7, width=90),
        dict(key="dpi", label="Resolution", kind="choice", choices=["300 dpi", "600 dpi"],
             value=1 if d["dpi"] >= 600 else 0, width=90)],
        note="Used by Copy image and Save image as in every plot.",
        ok="Apply")
    if v is None:
        return
    try:
        w = float(v["w"].replace(",", "."))
        h = float(v["h"].replace(",", "."))
    except ValueError:
        w, h = d["width_cm"], d["height_cm"]
    if v["preset"] < 3:
        w = sizes[v["preset"]]
    new = {"mode": "journal" if v["mode"] == 0 else "screen", "width_cm": max(2.0, min(w, 40.0)),
           "height_cm": max(1.5, min(h, 40.0)), "text": 7.0 + v["text"], "dpi": 600 if v["dpi"] == 1 else 300}
    T._save({"figure_export": new})
    _status(win, "Images: %s" % ("as on screen" if new["mode"] == "screen" else "%.1f x %.1f cm, text %g pt, %d dpi" % (
        new["width_cm"], new["height_cm"], new["text"], new["dpi"])))


def forward_wheel(win, e):
    """Scroll the column of tiles (RowScroller) that holds win."""
    p = win
    while p is not None and not isinstance(p, RowScroller):
        p = p.GetParent()
    if p is None:
        return
    rot, delta = e.GetWheelRotation(), e.GetWheelDelta() or 120
    if e.GetWheelAxis() != wx.MOUSE_WHEEL_VERTICAL or not rot:
        return
    lines = e.GetLinesPerAction() or 3
    x, y = p.GetViewStart()
    p.Scroll(-1, max(0, y - int(round(rot / float(delta) * lines))))


def _tile_title(w):
    try:
        if isinstance(w, SplitBox):
            return ", ".join(_tile_title(p) for p in w.all_panes if p) or "tiles"
        return getattr(w, "title", "") or "tile"
    except RuntimeError:
        return "tile"


def close_tile(win):
    """Hide a tile (plot or table) of a view until it is shown again
    (right click > Closed tiles)."""
    child, parent = win, win.GetParent()
    while parent is not None and not isinstance(parent, (SplitBox, StackBox)):
        child, parent = parent, parent.GetParent()
    if parent is None:
        return
    hook = getattr(parent, "on_close_tile", None)  # its owner may handle it (deconvolution results)
    if hook is not None and hook(child):
        return
    parent.close(child)
    if isinstance(parent, SplitBox) and not parent.panes:
        close_tile(parent)  # nothing left in this row: close the row
    _status(win, "Tile closed (right click > Closed tiles shows it again)")


def _stack_of(win):
    p = win
    while p is not None and not isinstance(p, StackBox):
        p = p.GetParent()
    return p


def _extra_rows(row):
    """The rows of an extra result (MSTab.set_extra): one window, a list of
    windows or None, as a list."""
    if row is None:
        return []
    return [r for r in (row if isinstance(row, (list, tuple)) else [row]) if r is not None]


def closed_tiles(win):
    """(title, reopen) for the closed tiles of the view that holds win."""
    sb = _stack_of(win)
    if sb is None:
        return []
    out = []
    for p in list(sb.closed):
        if p:
            out.append((_tile_title(p), lambda p=p: sb.reopen(p)))
    for row in sb.all_panes:
        if isinstance(row, SplitBox) and row:
            for p in list(row.closed):
                if p:
                    def reopen(row=row, p=p):
                        row.reopen(p)
                        if row in sb.closed:
                            sb.reopen(row)
                    out.append((_tile_title(p), reopen))
    return out


def tile_menu(win):
    """Close this tile, and the closed tiles to show again."""
    items = [(None, None), ("Close this tile", lambda: close_tile(win))]
    closed = closed_tiles(win)
    if closed:
        def all_():
            for t, fn in closed:
                fn()
        items.append(("Closed tiles", [(t, fn) for t, fn in closed] + [(None, None), ("Show them all", all_)]))
    sb = _stack_of(win)
    if sb is not None:
        journal = sb.mode() == "journal" and not sb.resized()
        sc = sb.scale()
        sizes = []
        if sb.has_fit():
            sizes.append(((sb.fit_label or "First %d tiles fill the window" % sb.fit_first)
                          + (" (default)" if sb.default_mode == "two" else ""), sb.two,
                          sb.mode() == "two" and not sb.resized()))
        sizes += [("Half size" + (" (default)" if sb.default_mode == "journal" else ""), lambda: sb.journal(0.5),
                   journal and abs(sc - 0.5) < 0.01),
                  ("Journal proportions", lambda: sb.journal(1.0), journal and abs(sc - 1.0) < 0.01),
                  ("Fit every tile in the window", sb.fit_all, sb.mode() == "fit" and not sb.resized())]
        items.append(("Tile heights", sizes))
    return items


def tile_shape(pane, width):
    """Height for journal proportions: the widest plot of the row gets the
    width : height of copied images (image settings); tables none."""
    cfg = export_settings()
    aspect = float(cfg["height_cm"]) / float(cfg["width_cm"])
    if isinstance(pane, SplitBox):
        if not pane.panes or not any(isinstance(p, PlotCard) for p in pane.panes):
            return None
        tw = sum(pane.weights) or 1.0
        wmax = max(w for p, w in zip(pane.panes, pane.weights)) / tw
        return width * wmax * aspect
    if isinstance(pane, (PlotCard, PairCard)):
        return width * aspect
    return None


def build_menu(win, items):
    """wx.Menu from a list of entries: (label, function), (label, function,
    ticked) for a tick item, (label, [entries]) for a submenu, (label, None)
    for a greyed line, (None, None) for a separator."""
    menu = wx.Menu()
    last_sep = True
    for it in items:
        label, what = it[0], it[1]
        if label is None:
            if not last_sep:
                menu.AppendSeparator()
                last_sep = True
            continue
        last_sep = False
        if isinstance(what, (list, tuple)):
            menu.AppendSubMenu(build_menu(win, what), label)
            continue
        import shortcuts
        label = shortcuts.menu_label(label)  # e.g. "Full view<tab>Home"
        if len(it) > 2:
            mi = menu.AppendCheckItem(wx.ID_ANY, label)
            mi.Check(bool(it[2]))
        else:
            mi = menu.Append(wx.ID_ANY, label)
        if what is None:
            mi.Enable(False)
        else:
            win.Bind(wx.EVT_MENU, lambda ev, f=what: f(), mi)
    n = menu.GetMenuItemCount()
    if n and menu.FindItemByPosition(n - 1).IsSeparator():
        menu.Delete(menu.FindItemByPosition(n - 1))
    return menu


def popup_menu(win, menu, pos=None):
    """Open a native menu without a stale tooltip covering its entries."""
    try:
        wx.ToolTip.Enable(False)
        return win.PopupMenu(menu) if pos is None else win.PopupMenu(menu, pos)
    finally:
        wx.ToolTip.Enable(True)


def _spaced_locator(spacing_pt=32.0):
    """Round-number y ticks about as far apart as the x ticks look: the
    number of ticks follows the height of the plot (about one per 30 pt, 2
    to 6), so short tiles and copied images do not get crowded ticks."""
    from matplotlib.ticker import MaxNLocator

    class SpacedLocator(MaxNLocator):
        def __call__(self):
            try:
                ax = self.axis.axes
                h_pt = ax.get_position().height * ax.figure.get_size_inches()[1] * 72.0
                self.set_params(nbins=int(max(2, min(6, round(h_pt / spacing_pt)))))
            except Exception:
                pass
            return MaxNLocator.__call__(self)

    # at least three labelled ticks, so the scale can always be read
    return SpacedLocator(nbins=4, steps=[1, 2, 2.5, 5, 10], min_n_ticks=3)


def _x_locator(steps=(1, 2, 2.5, 5, 10)):
    """At most 8 labelled ticks on the x axis, fewer where the plot is too narrow for their numbers
    (a copied journal sized image of a mass spectrum had its 5 digit numbers run into each other).
    steps: the tick steps allowed (times: without 2.5, which gave 2.5, 5.0, 7.5 min). dense (set on the
    locator): about twice as many numbers (the plots of a tile of several plots), still apart."""
    from matplotlib.ticker import MaxNLocator

    class WidthLocator(MaxNLocator):
        dense = False

        def __call__(self):
            try:
                ax = self.axis.axes
                # from the width of the figure, not of the axes: fit_axes moves the axes after the
                # ticks are measured, and a tick added then had its number cut off at the edge
                w_pt = 0.8 * ax.figure.get_size_inches()[0] * 72.0
                lo, hi = self.axis.get_view_interval()
                digits = len("%d" % max(abs(lo), abs(hi))) + (2 if abs(hi - lo) < 10 else 0)
                fs = float(self.axis._major_tick_kw.get("labelsize", 9) or 9)
                each = digits * 0.62 * fs + (1.2 if self.dense else 1.6) * fs  # a number and the gap to the next
                self.set_params(nbins=int(max(2, min(16 if self.dense else 8, w_pt // each))))
            except Exception:
                pass
            return MaxNLocator.__call__(self)
    return WidthLocator(nbins=8, steps=list(steps))


_TEXT_W = {}


def _text_w_pt(s, size, italic=False):
    """Width of a text in points at size (the plot font), measured without drawing."""
    key = (s, round(float(size), 2), bool(italic))
    w = _TEXT_W.get(key)
    if w is None:
        from matplotlib.textpath import text_to_path
        from matplotlib.font_manager import FontProperties
        w = text_to_path.get_text_width_height_descent(
            s, FontProperties(size=size, style="italic" if italic else "normal"), ismath=False)[0]
        if len(_TEXT_W) > 2000:
            _TEXT_W.clear()
        _TEXT_W[key] = w
    return w


_SUP = str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻")


def _side_formatter(card):
    """y numbers of a plot in a tile of several plots: divided by the power of ten of the axis, which is
    written in the y title of the plot (no ×10ⁿ above the axis); the title follows every zoom."""
    import weakref
    from matplotlib.ticker import ScalarFormatter

    class SideFormatter(ScalarFormatter):
        def set_locs(self, locs):
            ScalarFormatter.set_locs(self, locs)
            c = ref()
            if c is not None:
                c._side_power(self.orderOfMagnitude)

    ref = weakref.ref(card)
    f = SideFormatter(useOffset=False, useMathText=True)
    f.set_powerlimits((-3, 4))
    f._side = True
    return f


def _unit_formatter():
    """x numbers of a plot with its unit at the right end of the axis (on the line of the numbers): a number
    that the unit would cover is left out (worked out again at every drawing, so also after a zoom)."""
    from matplotlib.ticker import ScalarFormatter

    class UnitFormatter(ScalarFormatter):
        def __call__(self, x, pos=None):
            s = ScalarFormatter.__call__(self, x, pos)
            try:
                ax = self.axis.axes
                u = getattr(ax, "_unit_art", None)
                if not s or u is None or not u.get_visible() or not u.get_text() or u not in ax.texts:
                    return s
                k = ax.figure.dpi / 72.0
                ticks = self.axis.majorTicks
                fs = ticks[0].label1.get_fontsize() if ticks else 9.0
                px = ax.transData.transform((x, 0.0))[0]
                half = 0.5 * _text_w_pt(s, fs) * k
                uw = _text_w_pt(u.get_text(), u.get_fontsize(), u.get_style() == "italic") * k
                if px + half > ax.bbox.x1 - uw - 0.6 * fs * k:
                    return ""
            except Exception:
                pass
            return s

    f = UnitFormatter(useOffset=False)
    f.set_scientific(False)
    f._unit = True
    return f


def _titles_state(fig):
    """(title, text, size) of every plot title, to put them back after fit_titles."""
    out = []
    for ax in fig.axes:
        for t in (ax.title, ax._left_title, ax._right_title):
            out.append((t, t.get_text(), t.get_fontsize()))
    return out


def _restore_titles(state):
    for t, txt, fs in state:
        t.set_text(txt)
        t.set_fontsize(fs)


def fit_titles(fig, min_size=5.5):
    """Titles above a plot wider than the figure (the description of a deconvolution in a narrow
    copied image) get smaller type, down to min_size points, then a second line; they pushed the
    plot aside."""
    try:
        r = fig.canvas.get_renderer()
        W = fig.bbox.width
    except Exception:
        return
    for ax in fig.axes:
        try:
            e = ax.get_window_extent(r)
        except Exception:
            continue
        pad = 4.0 * fig.dpi / 72.0
        rooms = {id(ax._right_title): e.x1 - pad, id(ax._left_title): W - e.x0 - pad,
                 id(ax.title): 2.0 * min(0.5 * (e.x0 + e.x1), W - 0.5 * (e.x0 + e.x1)) - 2 * pad}
        for t in (ax.title, ax._left_title, ax._right_title):
            txt = t.get_text()
            if not txt or not t.get_visible():
                continue
            room = rooms[id(t)]  # a title runs from its anchor at the edge of the plot to the edge of the figure
            try:
                while t.get_window_extent(r).width > room and t.get_fontsize() > min_size:
                    t.set_fontsize(max(min_size, t.get_fontsize() - 0.5))
                if t.get_window_extent(r).width > room and "\n" not in txt:
                    # break after a comma or before a parenthesis, near the middle
                    cut = [i + 1 for i, ch in enumerate(txt) if ch == ","] + [i for i, ch in enumerate(txt) if ch == "("]
                    if cut:
                        i = min(cut, key=lambda j: abs(j - len(txt) / 2.0))
                        t.set_text(txt[:i].rstrip() + "\n" + txt[i:].lstrip())
            except Exception:
                pass


def style_axes(ax, xlabel="", ylabel=""):
    from matplotlib.ticker import AutoMinorLocator, ScalarFormatter, MaxNLocator
    for s in ax.spines.values():
        s.set_visible(True)
        s.set_linewidth(0.8)
        s.set_color(INK)
    ax.xaxis.set_ticks_position("bottom")
    ax.yaxis.set_ticks_position("left")
    # ticks outside the frame, only where the numbers are (bottom and left)
    ax.tick_params(which="both", direction="out", top=False, right=False, color=INK, width=0.8)
    ax.tick_params(which="major", length=4.0, labelsize=9, labelcolor=INK, pad=2.5)
    ax.tick_params(which="minor", length=2.2)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.xaxis.set_major_locator(_x_locator())
    # whole numbers on the x axis: a zoomed mass or m/z axis read "1 2 3 ... +1.135e4" (an offset)
    xfmt = ScalarFormatter(useOffset=False)
    xfmt.set_scientific(False)
    ax.xaxis.set_major_formatter(xfmt)
    ax.yaxis.set_major_locator(_spaced_locator())
    fmt = ScalarFormatter(useMathText=True)
    fmt.set_powerlimits((-3, 4))
    ax.yaxis.set_major_formatter(fmt)
    ax.yaxis.get_offset_text().set_fontsize(8.5)
    ax.set_xlabel(xlabel, fontsize=10, color=INK, labelpad=2)
    ax.set_ylabel(ylabel, fontsize=10, color=INK, labelpad=3)
    ax.set_facecolor("white")


def draw_peaks(ax, t, y, peaks, color, selected=None):
    """Shade integrated peaks, draw their baselines and label retention times."""
    arts = []
    order = sorted(range(len(peaks)), key=lambda k: -peaks[k]["height"])
    labelled = set(order[:8])
    if selected is not None:
        labelled.add(selected)
    for n, p in enumerate(peaks):
        i0, i1 = p["i0"], p["i1"]
        ts = t[i0:i1 + 1]
        ys = y[i0:i1 + 1]
        if len(ts) < 2:
            continue
        base = p["b0"] + (p["b1"] - p["b0"]) * (ts - ts[0]) / max(ts[-1] - ts[0], 1e-12)
        sel = selected is not None and n == selected
        arts.append(ax.fill_between(ts, base, np.maximum(ys, base), color=color, alpha=0.32 if sel else 0.16,
                                    lw=0, zorder=1))
        ln, = ax.plot([ts[0], ts[-1]], [p["b0"], p["b1"]], color=INK, lw=0.7, alpha=0.8, zorder=2)
        ln._no_scale = True
        arts.append(ln)
        for x, yb in ((ts[0], p["b0"]), (ts[-1], p["b1"])):
            tk, = ax.plot([x, x], [yb, y[i0] if x == ts[0] else y[i1]], color=INK, lw=0.6, alpha=0.6, zorder=2)
            tk._no_scale = True
        if n in labelled:
            arts.append(ax.annotate("%.2f" % p["rt"], (p["rt"], p["apex_y"]), xytext=(0, 3),
                                    textcoords="offset points", ha="center", va="bottom", fontsize=7,
                                    color=INK if sel else "#444444", fontweight="bold" if sel else "normal",
                                    zorder=4))
    return arts


def apex_of_profile(x, y, i):
    """m/z of the apex of a profile peak (Gaussian through the points above
    half height), not the highest data point: on a TOF spectrum the points
    are 3 to 5 ppm apart, so the highest point can be 2 ppm off the apex."""
    try:
        import ms_deconv
        lo, hi = max(0, i - 25), min(len(x), i + 26)
        return ms_deconv.refine_peak_mass(np.column_stack([x[lo:hi], y[lo:hi]]), float(x[i]))
    except Exception:
        return float(x[i])


def label_maxima(ax, x, y, n=6, fmt="%.1f", min_sep=None, xlim=None, refine=None):
    """Label the n highest local maxima inside the view (refine(x, y, i):
    position of the label, e.g. the apex of a profile peak)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if len(x) < 3:
        return []
    if xlim is not None:
        m = (x >= xlim[0]) & (x <= xlim[1])
        x, y = x[m], y[m]
        if len(x) < 3:
            return []
    import ms_deconv  # the selection runs in the library (Python fallback in ms_deconv)
    chosen = ms_deconv.label_maxima_index(x, y, n, min_sep)
    if not chosen:
        return []
    lo, hi = (xlim if xlim is not None else (x[0], x[-1]))
    out = []
    for i in chosen:
        xi = refine(x, y, i) if refine is not None else x[i]
        f = (xi - lo) / max(hi - lo, 1e-12)
        ha = "left" if f < 0.05 else "right" if f > 0.95 else "center"
        out.append(ax.annotate(fmt % xi, (xi, y[i]), xytext=(2 if ha == "left" else -2 if ha == "right" else 0, 3),
                               textcoords="offset points", ha=ha, va="bottom", fontsize=7, color=INK,
                               zorder=4))
    return out


# ==========================================================================
# widgets
# ==========================================================================
T.ICONS.setdefault("export", '<path d="M12 4v11M7 10l5 5 5-5M5 20h14"/>')
T.ICONS.setdefault("integrate", '<path d="M3 20h18M4 19c3 0 3-13 6-13s3 13 6 13"/><path d="M8 19v-6M12 19v-10"/>')
T.ICONS.setdefault("help", '<circle cx="12" cy="12" r="9"/><path d="M9.5 9.5a2.5 2.5 0 1 1 3.5 2.3c-.6.3-1 .9-1 1.6V14M12 17.2v.1"/>')
TOOL_ICONS = {
    "select": '<path d="M6 3l12 9-5.5 1.2L15.5 20l-2.6 1.2-3-6.8L6 18z"/>',
    "zoom": '<circle cx="10.5" cy="10.5" r="6.2"/><path d="M15.2 15.2L21 21M8 10.5h5M10.5 8v5"/>',
    "pan": '<path d="M12 2.5v19M2.5 12h19M12 2.5L9.3 5.2M12 2.5l2.7 2.7M12 21.5l-2.7-2.7M12 21.5l2.7-2.7'
           'M2.5 12l2.7-2.7M2.5 12l2.7 2.7M21.5 12l-2.7-2.7M21.5 12l-2.7 2.7"/>',
    "full": '<path d="M4 9V4h5M20 9V4h-5M4 15v5h5M20 15v5h-5"/>',
    "background": '<rect x="3.5" y="5.5" width="17" height="13" rx="1.5"/><path d="M3.5 13l7.5-7.5M6 18.5l13-13'
                  'M12 18.5l8.5-8.5"/>',
    "auto": '<path d="M2 20h20M3 19.5c2.5 0 3-11 5.5-11s3 11 5.5 11 2-5 4-5 2 5 3 5"/>'
            '<path d="M18 2.5l.9 2 2 .9-2 .9-.9 2-.9-2-2-.9 2-.9z"/>',
    "drag_peak": '<path d="M2 19.5h20M5 19.5c2.6 0 3.6-12 7-12s4.4 12 7 12"/><path d="M5 23v-7M19 23v-7"/>',
    "click_peak": '<path d="M2 19.5h20M4 19.5c3 0 4.2-12 8-12s5 12 8 12"/><circle cx="12" cy="3.6" r="1.8" fill="{c}"/>',
    "split": '<path d="M2 19.5h20M3 19.5c2 0 3-10 5.2-10s2.3 6 3.8 6 1.8-9 4.3-9 2.7 13 4.7 13"/>'
             '<path d="M12 3v2M12 7.5v2M12 12v2M12 16.5v3"/>',
    "delete": '<path d="M2 19.5h20M3 19.5c3 0 4-9 7-9s4 9 7 9"/><path d="M16 3.5l5 5M21 3.5l-5 5"/>',
    "clear": '<path d="M4 7h16M9 7V4.5h6V7M6.5 7l1 13h9l1-13M10 11v6M14 11v6"/>',
    "xic": '<path d="M2 20h11M4 20V10M7.5 20V5M11 20v-7"/><path d="M14 13h7.5M18.5 9.5l3.5 3.5-3.5 3.5"/>',
    "measure": '<path d="M4 18V7M20 18V7M4 12.5h16M7.5 9l-3.5 3.5L7.5 16M16.5 9l3.5 3.5-3.5 3.5"/>',
    "label": '<path d="M3 11.5V4h7.5L21 14.5 13.5 22z"/><circle cx="7.3" cy="8.3" r="1.5"/>',
    "panel": '<rect x="3" y="4.5" width="18" height="15" rx="2"/><path d="M15 4.5v15M17.3 9h1.4M17.3 12h1.4"/>',
    "calibrate": '<path d="M3 20h18M5 20V9M9 20V4M13 20v-8M17 20V7"/><path d="M3 16l6-5 4 3 8-7" '
                 'stroke-dasharray="2 2"/>',
}
for _k in ("panel", "calibrate", "xic", "measure", "label", "full"):
    T.ICONS.setdefault(_k, TOOL_ICONS[_k])


def flat(parent, label, kind="secondary", icon=None, tooltip=None, height=30, min_width=0, handler=None, padx=12):
    b = T._cls("FlatButton")(parent, wx.ID_ANY, label, icon=icon, kind=kind, tooltip=tooltip, height=height,
                             padx=padx, min_width=min_width)
    if handler is not None:
        b.Bind(wx.EVT_BUTTON, handler)
    return b


class Segmented(wx.Control):
    """Two or more segments in a pill; the active one is a white chip."""

    def __init__(self, parent, labels, on_change):
        wx.Control.__init__(self, parent, style=wx.BORDER_NONE)
        self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
        self.labels, self.on_change = labels, on_change
        self.active, self.hover = 0, -1
        self.enabled = [True] * len(labels)
        self.SetFont(ui_font(9.5, 600))
        dc = wx.ClientDC(self)
        dc.SetFont(self.GetFont())
        self.widths = [dc.GetTextExtent(l)[0] + self.FromDIP(28) for l in labels]
        T.keep_dip(self, "widths", then=self.Refresh)
        sz = wx.Size(sum(self.widths) + self.FromDIP(8), self.FromDIP(34))
        self.SetMinSize(sz)
        self.SetInitialSize(sz)
        self.Bind(wx.EVT_PAINT, self._paint)
        self.Bind(wx.EVT_LEFT_DOWN, self._click)
        self.Bind(wx.EVT_MOTION, self._motion)
        self.Bind(wx.EVT_LEAVE_WINDOW, lambda e: self._set_hover(-1))

    def _seg_at(self, x):
        xx = self.FromDIP(4)
        for i, w in enumerate(self.widths):
            if xx <= x < xx + w:
                return i
            xx += w
        return -1

    def _set_hover(self, i):
        if i != self.hover:
            self.hover = i
            self.Refresh()

    def _motion(self, e):
        self._set_hover(self._seg_at(e.GetX()))

    def _click(self, e):
        i = self._seg_at(e.GetX())
        if i >= 0 and self.enabled[i] and i != self.active:
            self.active = i
            self.Refresh()
            self.on_change(i)

    def set_active(self, i):
        self.active = i
        self.Refresh()

    def enable_segment(self, i, flag):
        self.enabled[i] = flag
        self.Refresh()

    def _paint(self, e):
        dc = wx.AutoBufferedPaintDC(self)
        dc.SetBackground(wx.Brush(self.GetParent().GetBackgroundColour()))
        dc.Clear()
        gc = T.crisp(dc)
        w, h = self.GetClientSize()
        gc.SetBrush(wx.Brush(wx.Colour("#E9EEF6")))
        gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(wx.Colour(C["line"])).Width(1.0)))
        # Windows 11 corners (the outer frame 5 px, the chosen view 4 px), not pills
        gc.DrawRoundedRectangle(0.5, 0.5, w - 1, h - 1, self.FromDIP(5))
        x = self.FromDIP(4)
        p = self.FromDIP(3)
        for i, (lab, sw) in enumerate(zip(self.labels, self.widths)):
            if i == self.active:
                gc.SetPen(wx.TRANSPARENT_PEN)
                gc.SetBrush(wx.Brush(wx.Colour(14, 28, 48, 22)))
                gc.DrawRoundedRectangle(x, p + 1.5, sw, h - 2 * p, self.FromDIP(4))
                gc.SetBrush(wx.Brush(wx.Colour("#FFFFFF")))
                gc.DrawRoundedRectangle(x, p, sw, h - 2 * p, self.FromDIP(4))
                col = C["accent_text"]
            else:
                col = C["muted"] if self.enabled[i] else "#A9B1BC"
                if i == self.hover and self.enabled[i]:
                    col = C["text"]
            gc.SetFont(ui_font(9.5, 700 if i == self.active else 500), wx.Colour(col))
            tw, th = gc.GetTextExtent(lab)
            gc.DrawText(lab, x + (sw - tw) / 2.0, (h - th) / 2.0)
            x += sw


class ToolButton(wx.Control):
    """Icon and label; a mode button stays pressed while its tool is active."""

    def __init__(self, parent, key, label, icon, tooltip, mode=True):
        wx.Control.__init__(self, parent, style=wx.BORDER_NONE)
        self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
        self.key, self.label, self.icon, self.mode = key, label, icon, mode
        self.active = self.hover = self.down = False
        self._bmps = {}
        self.SetFont(ui_font(9, 600))
        dc = wx.ClientDC(self)
        dc.SetFont(self.GetFont())
        tw = dc.GetTextExtent(label)[0] if label else 0
        self.full_w = self.FromDIP(16) + (self.FromDIP(6) + tw if label else 0) + self.FromDIP(16)
        T.keep_dip(self, "full_w")
        self.text = label
        sz = wx.Size(self.full_w, self.FromDIP(32))
        self.SetMinSize(sz)
        self.SetInitialSize(sz)
        self.SetToolTip(tooltip)
        for ev, fn in ((wx.EVT_PAINT, self._paint), (wx.EVT_ENTER_WINDOW, lambda e: self._h(True)),
                       (wx.EVT_LEAVE_WINDOW, lambda e: self._h(False)), (wx.EVT_LEFT_DOWN, self._down),
                       (wx.EVT_LEFT_UP, self._up)):
            self.Bind(ev, fn)
        self.on_click = None

    def _h(self, v):
        self.hover = v
        if not v:
            self.down = False
        self.Refresh()

    def _down(self, e):
        self.down = True
        self.Refresh()

    def _up(self, e):
        was = self.down
        self.down = False
        self.Refresh()
        if was and self.on_click:
            wx.CallAfter(self.on_click, self.key)

    def set_active(self, v):
        if v != self.active:
            self.active = v
            self.Refresh()

    def show_label(self, show):
        """Icon and label, or the icon only (a narrow window)."""
        self.label = self.text if show else ""
        w = self.full_w if show else self.FromDIP(32)
        self.SetMinSize(wx.Size(w, self.FromDIP(32)))
        self.SetSize(wx.Size(w, self.FromDIP(32)))
        self.Refresh()

    def _bmp(self, col):
        if col not in self._bmps:
            body = TOOL_ICONS.get(self.icon) or T.ICONS.get(self.icon, "")
            self._bmps[col] = T._bitmap_from_svg(body, self.FromDIP(16), col)
        return self._bmps[col]

    def _paint(self, e):
        dc = wx.AutoBufferedPaintDC(self)
        dc.SetBackground(wx.Brush(self.GetParent().GetBackgroundColour()))
        dc.Clear()
        gc = T.crisp(dc)
        w, h = self.GetClientSize()
        r = self.FromDIP(4)  # Windows 11 corners
        if self.active:
            gc.SetBrush(wx.Brush(wx.Colour(C["accent_bg"])))
            gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(wx.Colour(59, 130, 217, 150)).Width(1.0)))
            gc.DrawRoundedRectangle(0.5, 0.5, w - 1, h - 1, r)
        elif self.down or self.hover:
            gc.SetBrush(wx.Brush(wx.Colour(C["down"] if self.down else C["hover"])))
            gc.SetPen(wx.TRANSPARENT_PEN)
            gc.DrawRoundedRectangle(0.5, 0.5, w - 1, h - 1, r)
        col = C["accent_text"] if (self.active or self.hover) else "#3B4450"
        bmp = self._bmp(col)
        x = self.FromDIP(8)
        if bmp:
            gc.DrawBitmap(bmp, x, (h - self.FromDIP(16)) / 2.0, self.FromDIP(16), self.FromDIP(16))
        if self.label:
            gc.SetFont(ui_font(9, 700 if self.active else 600), wx.Colour(col))
            tw, th = gc.GetTextExtent(self.label)
            gc.DrawText(self.label, x + self.FromDIP(21), (h - th) / 2.0)


class ToolStrip(wx.Panel):
    """Row of tools above the plots. items: ("mode", key, label, icon, tip),
    ("action", key, label, icon, tip) or ("sep", group title)."""

    def __init__(self, parent, items, on_mode, on_action):
        wx.Panel.__init__(self, parent, style=wx.BORDER_NONE)
        self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
        self.SetBackgroundColour(wx.Colour(C["panel"]))
        self.on_mode, self.on_action = on_mode, on_action
        self.buttons = {}
        s = wx.WrapSizer(wx.HORIZONTAL)
        s.AddSpacer(self.FromDIP(6))
        for it in items:
            if it[0] == "sep":
                line = wx.Panel(self, size=wx.Size(1, self.FromDIP(22)))
                line.SetBackgroundColour(wx.Colour(C["line"]))
                s.Add(line, 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT | wx.RIGHT | wx.TOP | wx.BOTTOM, self.FromDIP(6))
                if it[1]:
                    st = wx.StaticText(self, label=it[1].upper())
                    st.SetFont(ui_font(7.5, 700))
                    st.SetForegroundColour(wx.Colour(C["faint"]))
                    s.Add(st, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT, self.FromDIP(4))
                continue
            kind, key, label, icon, tip = it
            b = ToolButton(self, key, label, icon, tip, mode=(kind == "mode"))
            b.on_click = self._clicked
            self.buttons[key] = b
            s.Add(b, 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT | wx.TOP | wx.BOTTOM, self.FromDIP(3))
        self.SetSizer(s)
        self.Bind(wx.EVT_PAINT, self._paint)
        self.Bind(wx.EVT_SIZE, lambda e: (self.Layout(), self.Refresh(), e.Skip()))

    def Layout(self):
        self.fit_labels()
        return wx.Panel.Layout(self)

    def fit_labels(self):
        """The tools in one row: in a narrow window the last tools lose
        their labels first (icon and tooltip stay)."""
        w = self.GetClientSize()[0]
        if w <= 0 or not self.buttons:
            return
        btns = [b for b in self.buttons.values() if b.IsShown() and b.text]
        need = sum(it.CalcMin().width for it in self.GetSizer().GetChildren() if it.IsShown())
        need += sum(b.full_w - b.GetMinSize()[0] for b in btns)
        short = set()
        for b in reversed(btns):
            if need <= w:
                break
            need -= b.full_w - self.FromDIP(32)
            short.add(b)
        for b in btns:
            if bool(b.label) == (b in short):
                b.show_label(b not in short)

    def _clicked(self, key):
        b = self.buttons[key]
        if b.mode:
            self.set_mode(key)
            self.on_mode(key)
        else:
            self.on_action(key)

    def set_mode(self, key):
        for k, b in self.buttons.items():
            if b.mode:
                b.set_active(k == key)

    def _paint(self, e):
        dc = wx.AutoBufferedPaintDC(self)
        dc.SetBackground(wx.Brush(wx.Colour(C["bg"])))
        dc.Clear()
        gc = T.crisp(dc)
        w, h = self.GetClientSize()
        gc.SetBrush(wx.Brush(wx.Colour(C["panel"])))
        gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(wx.Colour(C["line"])).Width(1.0)))
        gc.DrawRoundedRectangle(0.5, 0.5, w - 1.5, h - 1.5, self.FromDIP(12))


class SectionHeader(wx.Panel):
    def __init__(self, parent, title, colour, first=False, foldable=True):
        wx.Panel.__init__(self, parent, style=wx.BORDER_NONE)
        self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
        self.title, self.colour, self.first = title, colour, first
        self.collapsed = False
        self.foldable = foldable
        self.on_click = None
        self.SetMinSize(wx.Size(-1, self.FromDIP(34)))
        self.Bind(wx.EVT_PAINT, self._paint)
        if foldable:
            self.Bind(wx.EVT_LEFT_UP, lambda e: self.on_click and self.on_click())
            self.SetCursor(wx.Cursor(wx.CURSOR_HAND))

    def _paint(self, e):
        dc = wx.AutoBufferedPaintDC(self)
        dc.SetBackground(wx.Brush(wx.Colour(C["panel"])))
        dc.Clear()
        gc = T.crisp(dc)
        w, h = self.GetClientSize()
        if not self.first:
            gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(wx.Colour(C["line"])).Width(1.0)))
            gc.StrokeLine(0, 0.5, w, 0.5)
        gc.SetPen(wx.TRANSPARENT_PEN)
        gc.SetBrush(wx.Brush(wx.Colour(self.colour)))
        bh = self.FromDIP(16)
        gc.DrawRoundedRectangle(0, (h - bh) / 2.0 + 2, self.FromDIP(3), bh, 1.5)
        gc.SetFont(ui_font(9.5, 700), wx.Colour(C["text"]))
        tw, th = gc.GetTextExtent(self.title)
        gc.DrawText(self.title, self.FromDIP(14), (h - th) / 2.0 + 2)
        if not self.foldable:
            return
        # chevron: down when open, right when collapsed
        cx, cy, k = w - self.FromDIP(18), h / 2.0 + 2, self.FromDIP(4)
        gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(wx.Colour(C["muted"])).Width(1.4 * self.FromDIP(100) / 100.0)
                               .Cap(wx.CAP_ROUND)))
        if self.collapsed:
            pts = [wx.Point2D(cx - k / 2.0, cy - k), wx.Point2D(cx + k / 2.0, cy), wx.Point2D(cx - k / 2.0, cy + k)]
        else:
            pts = [wx.Point2D(cx - k, cy - k / 2.0), wx.Point2D(cx, cy + k / 2.0), wx.Point2D(cx + k, cy - k / 2.0)]
        gc.StrokeLines(pts)


class SidePanel(wx.ScrolledWindow):
    """Control column: sections of labelled rows, buttons and notes."""

    def __init__(self, parent, width=300, foldable=True):
        wx.ScrolledWindow.__init__(self, parent, style=wx.VSCROLL | wx.BORDER_NONE)
        self.foldable = foldable
        self.SetBackgroundColour(C["panel"])
        self.SetScrollRate(0, self.FromDIP(12))
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        outer = wx.BoxSizer(wx.HORIZONTAL)
        outer.Add(self.sizer, 1, wx.EXPAND)
        sb = wx.SystemSettings.GetMetric(wx.SYS_VSCROLL_X, self)
        outer.AddSpacer(sb if sb > 0 else self.FromDIP(17))  # room for the scroll bar
        self.SetSizer(outer)
        self.width = width
        self.fit_screen()
        self._first = True
        self.pad = self.FromDIP(14)
        T.keep_dip(self, "pad")

    def fit_screen(self):
        """Width for the window size: narrower on a small screen (the right
        margin of the widest rows is used then). True when it changed."""
        sb = max(wx.SystemSettings.GetMetric(wx.SYS_VSCROLL_X, self), 0)
        w = screen_dip(self, self.width, self.width - 28) + sb
        if self.GetMinSize().width == w:
            return False
        self.SetMinSize(wx.Size(w, -1))
        return True

    def section(self, title, colour, collapsed=None):
        """Starts a section; its header shows or hides the rows below it."""
        hdr = SectionHeader(self, title, colour, self._first, foldable=self.foldable)
        self.sizer.Add(hdr, 0, wx.EXPAND)
        self._first = False
        if not hasattr(self, "_headers"):
            self._headers = []
        self._headers.append(hdr)
        if not self.foldable:
            return hdr
        hdr.on_click = lambda h=hdr: self.toggle(h)
        saved = (T._load().get("side_collapsed") or {}).get(title)
        want = saved if saved is not None else bool(collapsed)
        if want:
            wx.CallAfter(self.toggle, hdr, True, False)
        return hdr

    def _section_items(self, hdr):
        items = list(self.sizer.GetChildren())
        idx = [i for i, it in enumerate(items) if it.GetWindow() is hdr]
        if not idx:
            return []
        out = []
        for it in items[idx[0] + 1:]:
            if isinstance(it.GetWindow(), SectionHeader):
                break
            out.append(it)
        return out

    def toggle(self, hdr, collapse=None, remember=True):
        try:
            collapse = (not hdr.collapsed) if collapse is None else collapse
        except RuntimeError:
            return
        if collapse == hdr.collapsed:
            return
        items = self._section_items(hdr)
        if collapse:
            hdr._was_shown = [it.IsShown() for it in items]
            for it in items:
                it.Show(False)
        else:
            for it, was in zip(items, getattr(hdr, "_was_shown", [True] * len(items))):
                it.Show(was)
        hdr.collapsed = collapse
        hdr.Refresh()
        if remember:
            d = dict(T._load().get("side_collapsed") or {})
            d[hdr.title] = collapse
            T._save({"side_collapsed": d})
        self.Layout()
        self.FitInside()

    def _label(self, text, muted=False):
        st = wx.StaticText(self, label=text)
        st.SetFont(ui_font(9, 400))
        st.SetForegroundColour(wx.Colour(C["muted"] if muted else C["text"]))
        return st

    def text(self, value="", width=64, tooltip=None):
        t = wx.TextCtrl(self, value=value, size=wx.Size(self.FromDIP(width), -1), style=wx.TE_PROCESS_ENTER)
        t.SetFont(ui_font(9, 400))
        if tooltip:
            t.SetToolTip(tooltip)
        return t

    def row(self, label, *items, stretch_first=False):
        """label on the left, controls (or strings for units) on the right."""
        s = wx.BoxSizer(wx.HORIZONTAL)
        if label:
            s.Add(self._label(label), 1 if not stretch_first else 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT,
                  self.FromDIP(6))
        for it in items:
            if isinstance(it, str):
                s.Add(self._label(it, muted=True), 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT | wx.RIGHT,
                      self.FromDIP(4))
            else:
                s.Add(it, 1 if (stretch_first and it is items[0]) else 0, wx.ALIGN_CENTER_VERTICAL)
        self.sizer.Add(s, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, self.pad)
        return s

    def full(self, ctrl):
        self.sizer.Add(ctrl, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, self.pad)
        return ctrl

    def buttons(self, *btns):
        s = wx.BoxSizer(wx.HORIZONTAL)
        for n, b in enumerate(btns):
            s.Add(b, 1, wx.EXPAND | (wx.LEFT if n else 0), self.FromDIP(6))
        self.sizer.Add(s, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, self.pad)
        return s

    def note(self, text):
        st = wx.StaticText(self, label=text)
        st.SetFont(ui_font(8.5, 400))
        st.SetForegroundColour(wx.Colour(C["muted"]))
        st.Wrap(self.FromDIP(270))
        self.sizer.Add(st, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, self.pad)
        return st

    def check(self, label, value=False):
        cb = wx.CheckBox(self, label=label)
        cb.SetFont(ui_font(9, 400))
        cb.SetValue(value)
        self.sizer.Add(cb, 0, wx.LEFT | wx.RIGHT | wx.TOP, self.pad)
        return cb

    def choice(self, items):
        c = wx.Choice(self, choices=items)
        c.SetFont(ui_font(9, 400))
        if items:
            c.SetSelection(0)
        return c

    def spin(self, value, lo, hi, width=64):
        s = wx.SpinCtrl(self, value=str(value), min=lo, max=hi, initial=value,
                        size=wx.Size(self.FromDIP(width), -1))
        s.SetFont(ui_font(9, 400))
        return s

    def gap(self, h=14):
        self.sizer.AddSpacer(self.FromDIP(h))


class PeakTable(wx.Panel):
    """Card with the integration results."""
    own_refresh = True   # (paints itself again when resized: place_panes)
    COLS = [("#", 34), ("Trace", 120), ("RT (min)", 98), ("Start", 62), ("End", 62), ("Height", 86), ("Area", 96),
            ("Area %", 62)]

    def __init__(self, parent, on_select, on_delete):
        wx.Panel.__init__(self, parent, style=wx.BORDER_NONE)
        self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
        self.SetBackgroundColour(C["bg"])
        self.on_select, self.on_delete = on_select, on_delete
        self.title, self.info = "Peaks", ""
        self.list = wx.ListCtrl(self, style=wx.LC_REPORT | wx.LC_SINGLE_SEL | wx.BORDER_NONE)
        self.list.SetFont(ui_font(9, 400))
        for n, (name, w) in enumerate(self.COLS):
            self.list.InsertColumn(n, name, wx.LIST_FORMAT_RIGHT if n > 1 else wx.LIST_FORMAT_LEFT,
                                   width=self.FromDIP(w))
        s = wx.BoxSizer(wx.VERTICAL)
        s.AddSpacer(self.FromDIP(30))
        s.Add(self.list, 1, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, self.FromDIP(10))
        self.SetSizer(s)
        self.Bind(wx.EVT_PAINT, self._paint)
        self.Bind(wx.EVT_SIZE, lambda e: (self.Layout(), refresh_card(self), e.Skip()))
        self.list.Bind(wx.EVT_LIST_ITEM_SELECTED, lambda e: self.on_select(e.GetIndex()))
        self.Bind(wx.EVT_MOUSEWHEEL, lambda e: forward_wheel(self, e))
        self.list.Bind(wx.EVT_LIST_KEY_DOWN, self._key)
        self.list.Bind(wx.EVT_CONTEXT_MENU, self._menu)

    def _key(self, e):
        if e.GetKeyCode() in (wx.WXK_DELETE, wx.WXK_BACK):
            i = self.list.GetFirstSelected()
            if i >= 0:
                self.on_delete(i)
        else:
            e.Skip()

    def _menu(self, e):
        i = self.list.GetFirstSelected()
        items = []
        if i >= 0:
            items.append(("Delete peak %d" % (i + 1), lambda: self.on_delete(i)))
        tab = self.GetParent()
        while tab is not None and not isinstance(tab, TabBase):
            tab = tab.GetParent()
        if tab is not None:
            items += [("Export the peak table\u2026", tab.on_export_table),
                      ("Remove every peak", tab.on_clear_peaks)]
        items += tile_menu(self)
        m = build_menu(self, items)
        popup_menu(self, m)
        m.Destroy()

    def set_peaks(self, peaks, units):
        self.list.DeleteAllItems()
        self.list.SetColumn(5, self._col("Height (%s)" % units if units else "Height", 5))
        self.list.SetColumn(6, self._col("Area (%s\u00b7s)" % units if units else "Area", 6))
        for n, p in enumerate(peaks):
            i = self.list.InsertItem(n, str(n + 1))
            vals = [p.get("trace", ""), "%.3f" % p["rt"], "%.3f" % p["t0"], "%.3f" % p["t1"], _g(p["height"]),
                    _g(p["area"]), "%.2f" % p["area_pct"]]
            for c, v in enumerate(vals, 1):
                self.list.SetItem(i, c, v)
        traces = sorted(set(p.get("trace", "") for p in peaks))
        self.info = ("%d peak%s in %s" % (len(peaks), "s" if len(peaks) != 1 else "", ", ".join(traces))
                     if peaks else "")
        self.Refresh()

    def _col(self, text, n):
        c = wx.ListItem()
        c.SetText(text)
        c.SetAlign(wx.LIST_FORMAT_RIGHT)
        c.SetMask(wx.LIST_MASK_TEXT | wx.LIST_MASK_FORMAT | wx.LIST_MASK_WIDTH)
        c.SetWidth(self.list.GetColumnWidth(n))
        return c

    def select(self, i):
        if 0 <= i < self.list.GetItemCount():
            self.list.Select(i)
            self.list.EnsureVisible(i)

    def _paint(self, e):
        paint_card(self, self.FromDIP(30), table_frame, self._draw_header,
                   (self.title, self.info, self.GetDPIScaleFactor()))

    def _draw_header(self, gc, w, font=ui_font):
        gc.SetFont(font(10, 700), wx.Colour(C["text"]))
        tw, th = gc.GetTextExtent(self.title)
        gc.DrawText(self.title, self.FromDIP(14), self.FromDIP(9))
        gc.SetFont(font(9, 400), wx.Colour(C["muted"]))
        txt = self.info or "No peaks"
        sw, sh = gc.GetTextExtent(txt)
        gc.DrawText(txt, self.FromDIP(14) + tw + self.FromDIP(10), self.FromDIP(9) + (th - sh) / 2.0)


def _g(v):
    if v == 0:
        return "0"
    a = abs(v)
    if a >= 1e6 or a < 1e-2:
        return "%.4g" % v
    if a >= 1000:
        return "%.0f" % v
    return "%.4g" % v


# ==========================================================================
# sizes that follow the screen (every view uses these rules)
# ==========================================================================
COMPACT_DIP = (1600, 900)  # a window narrower or lower than this (DIP, e.g. a laptop) uses the compact sizes
PLOT_MIN_DIP = 280  # least height of a plot tile in a view ...
PLOT_MIN_SHARE = 0.42  # ... but at most this share of the view height (two plots fill most of a small view)


def window_dip(win):
    """Client size (DIP) of the window that holds win; before that window
    has its size, the work area of its screen."""
    top = wx.GetTopLevelParent(win) if win else None
    try:
        s = float(win.GetDPIScaleFactor()) or 1.0
    except Exception:
        s = 1.0
    w = h = 0
    if top is not None:
        w, h = top.GetClientSize()
    if w < 600 * s or h < 400 * s:  # not sized yet
        try:
            i = wx.Display.GetFromWindow(top) if top is not None else 0
            a = wx.Display(max(0, i)).GetClientArea()
            w, h = a.width, a.height
        except Exception:
            return COMPACT_DIP
    return w / s, h / s


def is_compact(win):
    """The window is small (a laptop screen): compact side panels and bars."""
    w, h = window_dip(win)
    return w < COMPACT_DIP[0] or h < COMPACT_DIP[1]


def screen_dip(win, normal, compact):
    """A fixed size in pixels: normal DIP, or compact DIP in a small window."""
    return win.FromDIP(compact if is_compact(win) else normal)


def plot_floor(win, view_h):
    """Least height (pixels) of a plot tile in a view view_h pixels high:
    readable on a laptop without making large screens any different."""
    return int(min(win.FromDIP(PLOT_MIN_DIP), PLOT_MIN_SHARE * view_h))


def _is_plot_tile(p):
    if isinstance(p, (PlotCard, PairCard)):
        return True
    return isinstance(p, SplitBox) and any(isinstance(q, PlotCard) for q in p.panes)


def _at_once(method):
    """A method that changes the sizes of tiles, shown at once (change_at_once)."""
    def run(self, *a, **k):
        out = []
        change_at_once(self.GetTopLevelParent(), lambda: out.append(method(self, *a, **k)))
        return out[0] if out else None
    run.__name__, run.__doc__ = method.__name__, method.__doc__
    return run


def paint_now(win):
    """What is waiting to be painted in the window of win (and everything in it) painted now, in one pass."""
    try:
        top = win.GetTopLevelParent()
        if wx.Platform == "__WXMSW__":
            import ctypes
            ctypes.windll.user32.RedrawWindow(ctypes.c_void_p(top.GetHandle()), None, None,
                                              0x0080 | 0x0100)   # RDW_ALLCHILDREN | RDW_UPDATENOW
        else:
            top.Update()
    except (RuntimeError, OSError, AttributeError):
        pass


def refresh_own(win):
    """Only win's own pixels painted again (gaps between its tiles), not the windows inside it: they repaint
    themselves when their size changes, and moved ones keep their pixels."""
    if wx.Platform == "__WXMSW__":
        try:
            import ctypes
            ctypes.windll.user32.InvalidateRect(ctypes.c_void_p(win.GetHandle()), None, True)
            return
        except (RuntimeError, OSError, AttributeError):
            pass
    win.Refresh()


def live_drag_on():
    from ms_plot_render import _LIVE
    return _LIVE["on"]


def refresh_card(win):
    """A card resized: painted again; while an edge is dragged its list (a native control) is left to repaint
    what it needs itself."""
    if live_drag_on():
        refresh_own(win)
    else:
        win.Refresh()


def place_panes(box, rects):
    """The panes of box at their rectangles [(win, (x, y, w, h))], then box painted again. While an edge is
    dragged only what changes is painted: the panes given a new size and box's own gaps (moved panes keep their
    pixels)."""
    if not live_drag_on():
        for win, r in rects:
            win.SetSize(*r)
        box.Refresh()
        return
    deferred = box.BeginRepositioningChildren()   # (moved together, at the end)
    try:
        for win, r in rects:
            old = win.GetSize()
            win.SetSize(*r)
            if (r[2], r[3]) != (old.width, old.height) and not getattr(win, "own_refresh", False):
                win.Refresh()   # (cards and rows paint themselves again when resized)
    finally:
        if deferred:
            box.EndRepositioningChildren()
    refresh_own(box)


class LiveDrag(object):
    """An edge being dragged resizes at once: apply(arg) runs for the newest mouse position (at most every
    STEP_S; a later position waits for a timer) and is painted at once (input never waits for the painting).
    Plots show their last picture stretched while the mouse moves and are drawn sharp when it rests or is
    released."""
    STEP_S = 0.006

    def __init__(self, win, apply):
        from ms_plot_render import set_live
        self.win, self.apply, self.arg = win, apply, None
        self._step_t = self._rest_t = None
        self._last = 0.0
        self._busy = False
        set_live(True)

    def move(self, arg):
        self.arg = arg
        if self._busy:
            return
        wait = self.STEP_S - (time.perf_counter() - self._last)
        if wait <= 0:
            if self._step_t is not None:
                self._step_t.Stop()
                self._step_t = None
            self._step()
        elif self._step_t is None:
            self._step_t = wx.CallLater(max(1, int(wait * 1000 + 0.5)), self._step)

    def _step(self):
        self._step_t = None
        self._last = time.perf_counter()
        self._busy = True
        try:
            self.apply(self.arg)
            paint_now(self.win)
        except RuntimeError:
            return  # closed
        finally:
            self._busy = False
        if self._rest_t is None:
            self._rest_t = wx.CallLater(250, self._sharp)
        else:
            self._rest_t.Start(250)

    def _sharp(self):
        self._rest_t = None
        from ms_plot_render import flush_live
        try:
            change_at_once(self.win.GetTopLevelParent(), flush_live)
        except RuntimeError:
            pass

    def end(self):
        """Stops; the plots get their sizes (drawn by draw_pending_plots when called inside change_at_once)."""
        from ms_plot_render import set_live
        for t in (self._step_t, self._rest_t):
            if t is not None:
                t.Stop()
        self._step_t = self._rest_t = None
        set_live(False)


class SplitBox(wx.Panel):
    """Panes in a column (or a row) separated by gaps. Drag a gap to resize
    the two panes next to it (they follow the mouse)."""

    own_refresh = True   # (paints itself again when resized: place_panes)

    def __init__(self, parent, orient=wx.VERTICAL, gap=10, min_size=60):
        wx.Panel.__init__(self, parent, style=wx.BORDER_NONE)
        self.SetBackgroundColour(C["bg"])
        self.orient = orient
        self.gap = self.FromDIP(gap)
        self.min = self.FromDIP(min_size)
        T.keep_dip(self, "gap", "min", "view_h", then=lambda: self.relayout())
        self.panes, self.weights = [], []  # the panes shown (closed ones left out)
        self.all_panes, self.all_weights = [], []
        self.closed = set()
        self.on_change = None  # callback(weights) while and after a drag
        self._drag = None
        self._live = None
        self._delta = 0
        self._hover = None
        for ev, fn in ((wx.EVT_SIZE, lambda e: self.relayout()), (wx.EVT_LEFT_DOWN, self._down),
                       (wx.EVT_MOTION, self._motion), (wx.EVT_LEFT_UP, self._up),
                       (wx.EVT_MOUSE_CAPTURE_LOST, self._lost), (wx.EVT_LEAVE_WINDOW, self._leave),
                       (wx.EVT_PAINT, self._paint), (wx.EVT_MOUSEWHEEL, lambda e: forward_wheel(self, e))):
            self.Bind(ev, fn)

    def set_panes(self, panes, weights):
        self.all_panes = list(panes)
        self.all_weights = [max(float(w), 0.05) for w in weights]
        self.closed = {p for p in self.closed if p in self.all_panes}
        self._visible()
        self.relayout()

    def _visible(self):
        self.panes = [p for p in self.all_panes if p not in self.closed]
        self.weights = [w for p, w in zip(self.all_panes, self.all_weights) if p not in self.closed]
        for p in self.all_panes:
            try:
                p.Show(p not in self.closed)
            except RuntimeError:
                pass

    @_at_once
    def close(self, pane):
        self._sync()
        self.closed.add(pane)
        self._visible()
        self.relayout()

    @_at_once
    def reopen(self, pane):
        self._sync()  # widths dragged since the tile was closed are kept
        self.closed.discard(pane)
        self._visible()
        self.relayout()

    def _sync(self):
        for p, w in zip(self.panes, self.weights):
            if p in self.all_panes:
                self.all_weights[self.all_panes.index(p)] = w

    def weight_of(self, win, default):
        for p, w in zip(self.panes, self.weights):
            if p is win:
                return w
        return default

    def _vertical(self):
        return self.orient == wx.VERTICAL

    def _sizes(self):
        n = len(self.panes)
        if not n:
            return []
        w, h = self.GetClientSize()
        total = max(0, (h if self._vertical() else w) - self.gap * (n - 1))
        tw = sum(self.weights) or 1.0
        sizes = [int(total * x / tw) for x in self.weights]
        sizes[-1] += total - sum(sizes)
        return sizes

    def relayout(self):
        w, h = self.GetClientSize()
        pos = 0
        rects = []
        for win, size in zip(self.panes, self._sizes()):
            if self._vertical():
                rects.append((win, (0, pos, w, max(size, 1))))
            else:
                rects.append((win, (pos, 0, max(size, 1), h)))
            pos += size + self.gap
        place_panes(self, rects)

    def _coord(self, e):
        return e.GetY() if self._vertical() else e.GetX()

    def _gap_at(self, c):
        pos = 0
        sizes = self._sizes()
        for i, size in enumerate(sizes[:-1]):
            pos += size
            if pos - 2 <= c <= pos + self.gap + 2:
                return i
            pos += self.gap
        return None

    def _down(self, e):
        i = self._gap_at(self._coord(e))
        if i is None:
            e.Skip()
            return
        self._drag = (i, self._coord(e), self._sizes(), list(self.weights))
        self._delta = 0
        if not self.HasCapture():
            self.CaptureMouse()
        self._live = LiveDrag(self, self._apply_drag)

    def _new_sizes(self, delta):
        i, c0, sizes = self._drag[:3]
        pair = sizes[i] + sizes[i + 1]
        if pair < 2 * self.min:
            return sizes[i], sizes[i + 1]  # too small for two panes of the least size: unchanged
        a = min(max(sizes[i] + delta, self.min), pair - self.min)
        return a, pair - a

    def _apply_drag(self, delta):
        """The two panes at the sizes for a drag by delta (from the sizes when the drag began)."""
        if self._drag is None:
            return
        i, c0, sizes, w0 = self._drag
        a, b = self._new_sizes(delta)
        self.weights = list(w0)
        pair = w0[i] + w0[i + 1]
        if a + b > 0 and (a, b) != (sizes[i], sizes[i + 1]):
            self.weights[i] = pair * a / float(a + b)
            self.weights[i + 1] = pair - self.weights[i]
        self.relayout()
        if self.on_change:
            self.on_change(self.weights)

    def _motion(self, e):
        c = self._coord(e)
        if self._drag is not None:
            self._delta = c - self._drag[1]
            if self._live is not None:
                self._live.move(self._delta)
            return
        i = self._gap_at(c)
        cur = (wx.CURSOR_SIZENS if self._vertical() else wx.CURSOR_SIZEWE) if i is not None else wx.CURSOR_ARROW
        self.SetCursor(wx.Cursor(cur))
        if i != self._hover:
            self._hover = i
            refresh_own(self)   # (the grip is drawn in the gap: the tiles stay)

    @_at_once
    def _finish(self, delta):
        self._apply_drag(delta)
        self._drag = None
        if self._live is not None:
            self._live.end()
            self._live = None
        self.Refresh()   # (while dragging only the parts that changed were painted)

    def _up(self, e):
        if self.HasCapture():
            self.ReleaseMouse()
        if self._drag is not None:
            self._finish(self._coord(e) - self._drag[1])
        else:
            e.Skip()

    def _lost(self, e):
        if self._drag is not None:
            self._finish(self._delta)  # kept as last shown

    def _leave(self, e):
        if self._drag is None and self._hover is not None:
            self._hover = None
            self.SetCursor(wx.Cursor(wx.CURSOR_ARROW))
            refresh_own(self)

    def _paint(self, e):
        dc = wx.PaintDC(self)
        if self._hover is None or self._drag is not None:
            return
        sizes = self._sizes()
        pos = sum(sizes[:self._hover + 1]) + self.gap * self._hover
        w, h = self.GetClientSize()
        L, t = self.FromDIP(40), max(2, self.FromDIP(4))
        dc.SetPen(wx.TRANSPARENT_PEN)
        dc.SetBrush(wx.Brush(wx.Colour(C["line2"])))
        mid = pos + (self.gap - t) // 2
        if self._vertical():
            dc.DrawRoundedRectangle((w - L) // 2, mid, L, t, t // 2)
        else:
            dc.DrawRoundedRectangle(mid, (h - L) // 2, t, L, t // 2)


# ==========================================================================
# shared tab logic: chromatogram views, ranges, integration
# ==========================================================================
class StackBox(wx.Panel):
    """Tiles in a column, each with its own height. Drag the gap below a tile
    to make it taller or shorter: the column may then be longer than the
    window and scrolls (it sits in a RowScroller). Double click a gap, or
    Fit tiles in the right click menu, and every tile fits the window again.
    Heights are kept as fractions of the window height, so they follow when
    the window is resized."""

    def __init__(self, parent, gap=10, min_size=100):
        wx.Panel.__init__(self, parent, style=wx.BORDER_NONE)
        self.SetBackgroundColour(C["bg"])
        self.gap = self.FromDIP(gap)
        self.min = self.FromDIP(min_size)
        T.keep_dip(self, "gap", "min", "view_h", then=lambda: self.relayout())
        self.panes, self.weights = [], []  # the panes shown (closed ones left out)
        self.all_panes, self.all_weights = [], []
        self.closed = set()
        self.frac = {}  # pane -> height / window height, after the user resized tiles
        self.shape = None  # callback(pane, width) -> height wanted for journal proportions, or None
        self.fit_first = 0  # mode "two": this many first tiles fill the window, the others below at that size
        # mode "two" instead by share: callback(pane) -> its share of the window, or None for a tile below (the
        # tiles with a share fill the window, the others follow below); fit_label names the mode in the menu
        self.fit_share = None
        self.fit_label = None
        self.fit_floor = 0  # DIP: with fit_share, least height of each share (a small window then scrolls)
        self.mode_key, self.default_mode = "tile_mode", "journal"  # setting of the tile mode
        self.view_h = self.FromDIP(600)
        self.on_change = None  # callback() while and after a drag (the column length changed)
        self._drag = None
        self._live = None
        self._delta = 0
        self._hover = None
        for ev, fn in ((wx.EVT_SIZE, lambda e: self.relayout()), (wx.EVT_LEFT_DOWN, self._down),
                       (wx.EVT_LEFT_DCLICK, self._dclick), (wx.EVT_MOTION, self._motion),
                       (wx.EVT_LEFT_UP, self._up), (wx.EVT_MOUSE_CAPTURE_LOST, self._lost),
                       (wx.EVT_LEAVE_WINDOW, self._leave), (wx.EVT_PAINT, self._paint),
                       (wx.EVT_MOUSEWHEEL, lambda e: forward_wheel(self, e))):
            self.Bind(ev, fn)

    def set_panes(self, panes, weights):
        self.all_panes = list(panes)
        self.all_weights = [max(float(w), 0.05) for w in weights]
        self.closed = {p for p in self.closed if p in self.all_panes}
        self.frac = {p: f for p, f in self.frac.items() if p in self.all_panes}
        self._update_visible()
        self.relayout()

    def _update_visible(self):
        self.panes = [p for p in self.all_panes if p not in self.closed]
        self.weights = [w for p, w in zip(self.all_panes, self.all_weights) if p not in self.closed]
        for p in self.all_panes:
            try:
                p.Show(p not in self.closed)
            except RuntimeError:
                pass

    @_at_once
    def close(self, pane):
        self.closed.add(pane)
        self._update_visible()
        if self.on_change:
            self.on_change()
        self.relayout()

    @_at_once
    def reopen(self, pane):
        self.closed.discard(pane)
        self._update_visible()
        if self.on_change:
            self.on_change()
        self.relayout()

    def weight_of(self, win, default):
        for p, w in zip(self.all_panes, self.all_weights):
            if p is win:
                return w
        return default

    def mode(self):
        m = settings().get(self.mode_key, self.default_mode)
        return m if (m != "two" or self.has_fit()) else "journal"

    def has_fit(self):
        """Mode "two" can be chosen: some tiles fill the window."""
        if self.fit_share is not None:
            return any(self._shares())
        return bool(self.fit_first)

    def _shares(self, bases=None):
        """Share of the window of each tile shown in mode "two" (0 for the tiles below). fit_share may give
        (share, pixels): the tile gets these pixels (its title and axes) and its share of what is left;
        bases, a list, receives these pixels."""
        if self.fit_share is None:
            out = [1.0 if i < self.fit_first else 0.0 for i in range(len(self.panes))]
            if bases is not None:
                bases[:] = [0] * len(out)
            return out
        out, base = [], []
        for p in self.panes:
            try:
                s = self.fit_share(p)
            except Exception:
                s = None
            s, b = s if isinstance(s, tuple) else (s, 0)
            ok = bool(s) and float(s) > 0
            out.append(float(s) if ok else 0.0)
            base.append(int(b) if ok else 0)
        if bases is not None:
            bases[:] = base
        return out

    @staticmethod
    def scale():
        """Tile height as a share of the journal proportions (0.5 = half)."""
        try:
            return min(1.5, max(0.2, float(settings().get("tile_scale", 0.5))))
        except (TypeError, ValueError):
            return 0.5

    def heights(self, H=None, W=None):
        """Tile heights for a window H high (view_h) and a column W wide (its width)."""
        H = self.view_h if H is None else H
        Wc = self.GetClientSize()[0] if W is None else W
        n = len(self.panes)
        if not n:
            return []
        avail = max(0, H - self.gap * (n - 1))
        out = [None] * n
        for i, p in enumerate(self.panes):  # tiles the user sized
            f = self.frac.get(p)
            if f is not None:
                out[i] = int(max(self.min, f * H))
        if self.mode() == "two":
            # the first tiles fill the window (e.g. chromatogram and spectrum; with fit_share the tiles given a
            # share, each its share of the window); the others follow below at the same height, tables smaller
            bases = []
            shares = self._shares(bases)
            top = [i for i in range(n) if shares[i] > 0]
            k = len(top)
            unit = int(max(self.min, (H - self.gap * (k - 1)) / float(k)))
            fixed = sum(out[i] for i in top if out[i] is not None)
            free = [i for i in top if out[i] is None]
            if free:
                rest = H - self.gap * (k - 1) - fixed
                ts = sum(shares[i] for i in free)
                left = max(0, rest - sum(bases[i] for i in free))
                each = max(left / ts, self.FromDIP(self.fit_floor) if self.fit_share is not None else 0)
                for i in free:
                    out[i] = int(max(self.min, bases[i] + each * shares[i]))
                if self.fit_share is not None and sum(out[i] for i in free) < rest:
                    out[free[-1]] += rest - sum(out[i] for i in free)  # (the last one ends at the window's end)
            W = max(Wc, self.FromDIP(300))
            for i in range(n):
                if i in top:
                    continue
                if out[i] is None:
                    try:
                        want = self.shape(self.panes[i], W) if self.shape is not None else 0
                    except Exception:
                        want = 0
                    out[i] = int(max(self.min, 0.30 * H)) if want is None else unit
            return out
        if self.shape is not None and self.mode() == "journal":
            # journal proportions: each tile as high as a copied image of it
            # would be for its width (the view scrolls); tables get a fixed share
            W = max(Wc, self.FromDIP(300))
            sc = self.scale()
            low = plot_floor(self, H)  # on a small screen a plot is never lower than this
            for i, p in enumerate(self.panes):
                if out[i] is None:
                    try:
                        want = self.shape(p, W)
                    except Exception:
                        want = None
                    if want is None:  # tables: their own share, not scaled (rows stay readable)
                        out[i] = int(max(self.min, 0.30 * H))
                        continue
                    # at most about 60 % of the window (so the next tile shows
                    # below it), both scaled by the tile size (half by default)
                    # a tile can ask for its own size (HRMS: the deconvolution
                    # result rows at full size while the others are at half)
                    sc_i = getattr(p, "tile_scale", None) or sc
                    if isinstance(p, PairCard):
                        # several plots in one tile: half the window at half size (the chromatograms and the
                        # spectra of both polarities in one window), the whole window at journal proportions
                        out[i] = int(max(self.min, low, sc_i * (H - self.gap)))
                        continue
                    out[i] = int(max(self.min, low, sc_i * min(want, 0.62 * H)))
            return out
        free = [i for i in range(n) if out[i] is None]
        rest = max(0, avail - sum(h for h in out if h is not None))
        # the others share the rest by weight; a tile below the minimum gets
        # the minimum and the rest is shared again (so they fit if possible)
        while free:
            tw = sum(self.weights[i] for i in free) or 1.0
            small = [i for i in free if rest * self.weights[i] / tw < self.min]
            if not small or len(small) == len(free):
                for i in free:
                    out[i] = int(max(self.min, rest * self.weights[i] / tw))
                break
            for i in small:
                out[i] = self.min
                rest = max(0, rest - self.min)
            free = [i for i in free if i not in small]
        if self.shape is None:  # a view without tile modes (Compare): its plots readable, the rest scrolls
            low = plot_floor(self, H)
            for i, p in enumerate(self.panes):
                if self.frac.get(p) is None and _is_plot_tile(p):
                    out[i] = max(out[i], low)
        return out

    def total(self, H=None, W=None):
        hs = self.heights(H, W)
        return sum(hs) + self.gap * max(0, len(hs) - 1)

    def relayout(self):
        w = self.GetClientSize()[0]
        y = 0
        rects = []
        for win, hgt in zip(self.panes, self.heights()):
            rects.append((win, (0, y, w, max(hgt, 1))))
            y += hgt + self.gap
        place_panes(self, rects)

    @_at_once
    def fit_all(self):
        self.frac = {}
        T._save({self.mode_key: "fit"})
        if self.on_change:
            self.on_change()
        self.relayout()

    @_at_once
    def two(self):
        self.frac = {}
        T._save({self.mode_key: "two"})
        if self.on_change:
            self.on_change()
        self.relayout()

    @_at_once
    def journal(self, scale=None):
        self.frac = {}
        upd = {self.mode_key: "journal"}
        if scale is not None:
            upd["tile_scale"] = float(scale)
        T._save(upd)
        if self.on_change:
            self.on_change()
        self.relayout()

    def resized(self):
        return bool(self.frac)

    def _gap_at(self, c):
        pos = 0
        for i, hgt in enumerate(self.heights()):
            pos += hgt
            if pos - 2 <= c <= pos + self.gap + 2:
                return i
            pos += self.gap
        return None

    def _visible(self):
        """Part of the column inside the window (column coordinates)."""
        y0 = -self.GetPosition()[1]
        return y0, y0 + self.view_h

    def _down(self, e):
        i = self._gap_at(e.GetY())
        if i is None:
            e.Skip()
            return
        self._drag = (i, e.GetY(), self.heights())
        self._delta = 0
        if not self.HasCapture():
            self.CaptureMouse()
        self._live = LiveDrag(self, self._apply_drag)

    def _dclick(self, e):
        if self._gap_at(e.GetY()) is not None:
            self.fit_all()
        else:
            e.Skip()

    def _new_height(self, delta):
        i, c0, hs = self._drag
        top = sum(hs[:i]) + self.gap * i
        v0, v1 = self._visible()
        h = max(self.min, hs[i] + delta)
        return int(min(h, max(self.min, v1 - top - self.gap)))  # the gap stays in the window

    def _apply_drag(self, delta):
        """The tile at the height for a drag by delta (the others as when the drag began)."""
        if self._drag is None:
            return
        i, c0, hs = self._drag
        H = float(max(self.view_h, 1))
        self.frac = {p: hgt / H for p, hgt in zip(self.panes, hs)}
        self.frac[self.panes[i]] = self._new_height(delta) / H
        if self.on_change:
            self.on_change()
        self.relayout()

    def _motion(self, e):
        if self._drag is not None:
            self._delta = e.GetY() - self._drag[1]
            if self._live is not None:
                self._live.move(self._delta)
            return
        i = self._gap_at(e.GetY())
        self.SetCursor(wx.Cursor(wx.CURSOR_SIZENS if i is not None else wx.CURSOR_ARROW))
        if i != self._hover:
            self._hover = i
            refresh_own(self)   # (the grip is drawn in the gap: the tiles stay)

    @_at_once
    def _finish(self, delta):
        self._apply_drag(delta)
        self._drag = None
        if self._live is not None:
            self._live.end()
            self._live = None
        self.Refresh()   # (while dragging only the parts that changed were painted)

    def _up(self, e):
        if self.HasCapture():
            self.ReleaseMouse()
        if self._drag is not None:
            self._finish(e.GetY() - self._drag[1])
        else:
            e.Skip()

    def _lost(self, e):
        if self._drag is not None:
            self._finish(self._delta)  # kept as last shown

    def _leave(self, e):
        if self._drag is None and self._hover is not None:
            self._hover = None
            self.SetCursor(wx.Cursor(wx.CURSOR_ARROW))
            refresh_own(self)

    def _paint(self, e):
        dc = wx.PaintDC(self)
        if self._hover is None or self._drag is not None:
            return
        hs = self.heights()
        pos = sum(hs[:self._hover + 1]) + self.gap * self._hover
        w = self.GetClientSize()[0]
        L, t = self.FromDIP(40), max(2, self.FromDIP(4))
        dc.SetPen(wx.TRANSPARENT_PEN)
        dc.SetBrush(wx.Brush(wx.Colour(C["line2"])))
        dc.DrawRoundedRectangle((w - L) // 2, pos + (self.gap - t) // 2, L, t, t // 2)


class RowScroller(wx.ScrolledWindow):
    """Holds the column of tiles (a StackBox) with a vertical scroll bar:
    tiles fit the window until one is made taller (or there are too many
    to fit at a readable height); then the column scrolls. The mouse wheel
    scrolls it anywhere over the tiles."""

    def __init__(self, parent, min_row=110):
        wx.ScrolledWindow.__init__(self, parent, style=wx.VSCROLL | wx.BORDER_NONE)
        self.SetBackgroundColour(C["bg"])
        self.SetScrollRate(0, self.FromDIP(24))
        self.min_row = self.FromDIP(min_row)
        T.keep_dip(self, "min_row", then=lambda: (self.SetScrollRate(0, self.FromDIP(24)), self.refit()))
        self.box = None
        self._busy = False
        self.Bind(wx.EVT_SIZE, self._size)

    def set_box(self, box):
        self.box = box
        box.min = max(box.min, self.min_row)
        box.on_change = self.refit
        self.refit()

    def _size(self, e):
        self.refit()
        e.Skip()

    def refit(self):
        if self.box is None or self._busy:
            return
        if live_drag_on():
            self._refit_once()
            return
        self._busy = True
        try:
            w, h = self.GetClientSize()
            self.box.view_h = max(h, 1)
            x, y = self.CalcScrolledPosition(0, 0)
            self.box.SetSize(x, y, w, max(self.box.GetSize()[1], 1))  # width first: journal heights follow it
            H = max(h, self.box.total(h))
            self.SetVirtualSize(w, H)
            w2, h2 = self.GetClientSize()  # the scroll bar may have taken some width
            if (w2, h2) != (w, h):
                self.box.view_h = max(h2, 1)
                self.box.SetSize(x, y, w2, max(self.box.GetSize()[1], 1))
                H = max(h2, self.box.total(h2))
                self.SetVirtualSize(w2, H)
            x, y = self.CalcScrolledPosition(0, 0)
            self.box.SetSize(x, y, w2, H)
            self.box.relayout()
        finally:
            self._busy = False


    def _refit_once(self):
        """As refit, with the column sized once (an edge being dragged): its length worked out for the new width
        first. The same sizes and places as refit."""
        self._busy = True
        try:
            box = self.box
            w, h = self.GetClientSize()
            box.view_h = max(h, 1)
            H = max(h, box.total(h, w))
            self.SetVirtualSize(w, H)
            w2, h2 = self.GetClientSize()  # the scroll bar may have taken some width
            if (w2, h2) != (w, h):
                box.view_h = max(h2, 1)
                H = max(h2, box.total(h2, w2))
                self.SetVirtualSize(w2, H)
            x, y = self.CalcScrolledPosition(0, 0)
            if tuple(box.GetRect()) != (x, y, w2, H):
                box.SetSize(x, y, w2, H)   # (laid out by its size event)
            else:
                box.relayout()
        finally:
            self._busy = False


class ChromView(object):
    """A chromatogram card and the traces drawn in it."""

    def __init__(self, tab, card, key):
        self.tab, self.card, self.key = tab, card, key
        self.traces = []
        card.on_range = lambda a, b, s: tab.chrom_range(a, b, s)
        card.on_pick = lambda x, y: tab.chrom_pick(x, y, self)
        card.on_menu = lambda x, y: tab.chrom_menu(x, y, self)
        card.on_view = lambda: tab.chrom_view_changed(self)
        card.readout_fmt = lambda x, y: tab.chrom_readout(x, y, self)
        card.tool_getter = lambda: tab.tool
        card.tools_supported = {"background", "drag_peak", "click_peak", "split", "delete", "calibrate"}
        card.on_tool = lambda tool, kind, a, b: tab.chrom_tool(self, tool, kind, a, b)


class TabBase(wx.Panel):
    ylabel = "Intensity"
    units = ""

    def __init__(self, parent, frame):
        wx.Panel.__init__(self, parent, style=wx.BORDER_NONE)
        self.SetBackgroundColour(C["bg"])
        self.frame = frame
        self.views = []
        self.peaks = []
        self.sel_peak = None
        self.avg = None  # (t0, t1) selected range
        self.bg = None  # (b0, b1) background range
        self.link = None
        self.tool = "select"
        self.side = SidePanel(self, 318)

    def tree_changed(self):
        """The traces or spectra of this view changed: the tree of the files panel is drawn again (once, after
        the change)."""
        fp = getattr(self.frame, "files", None)
        if fp is not None and not getattr(fp, "_tree_due", False):
            fp._tree_due = True
            wx.CallAfter(fp.tree_refresh)

    def status(self, text):
        try:
            self.frame.SetStatusText(text, 0)
        except Exception:
            pass

    def all_traces(self):
        return [tr for v in self.views for tr in v.traces]

    # ---------------------------------------------------------------- tools
    TOOL_HINTS = {
        "select": "Select: drag a time range or click a time",
        "zoom": "Zoom: drag a box on any plot",
        "pan": "Pan: drag a plot to move the view",
        "background": "Background: drag the background range",
        "drag_peak": "Drag peak: drag from the start to the end of a peak",
        "click_peak": "Click peak: click near the top of a peak",
        "split": "Split: click inside a peak to divide it",
        "delete": "Delete: click a peak to remove it",
        "xic": "Mass chromatogram: click a peak in a spectrum",
        "measure": "Measure: click two peaks in a spectrum",
        "label": "Label: click a peak to pin or unpin its label",
    }

    def tool_items(self, spectrum_tools):
        items = [("mode", "select", "", "select", "Select: drag a time range or click a time (Esc)"),
                 ("mode", "zoom", "", "zoom", "Zoom: drag a box on any plot"),
                 ("mode", "pan", "", "pan", "Pan: drag to move the view"),
                 ("action", "full", "", "full", "Full view of every plot"),
                 ("sep", ""),
                 ("mode", "background", "Background", "background", self.TOOL_HINTS["background"]),
                 ("sep", ""),
                 ("action", "auto", "Auto", "auto", "Integrate every peak in the visible window"),
                 ("mode", "drag_peak", "Drag", "drag_peak", self.TOOL_HINTS["drag_peak"]),
                 ("mode", "click_peak", "Click", "click_peak", self.TOOL_HINTS["click_peak"]),
                 ("mode", "split", "Split", "split", self.TOOL_HINTS["split"]),
                 ("mode", "delete", "Delete", "delete", self.TOOL_HINTS["delete"]),
                 ("action", "clear", "", "clear", "Remove every integrated peak")]
        if spectrum_tools:
            items += [("sep", "")] + spectrum_tools
        return items

    def build_tools(self, spectrum_tools):
        self.tools = ToolStrip(self, self.tool_items(spectrum_tools), self.set_tool, self.tool_action)
        self.tools.set_mode("select")
        import ms_ribbon
        if ms_ribbon.enabled(self.frame):
            self.tools.Hide()  # (its tools are in the ribbon of the window)
        import shortcuts
        shortcuts.decorate_tools(self)  # e.g. "Zoom: drag a box on any plot (Z)"
        self.Bind(wx.EVT_CHAR_HOOK, self._keys)
        return self.tools

    def _keys(self, e):
        if e.GetKeyCode() == wx.WXK_ESCAPE and self.tool != "select":
            self.tools.set_mode("select")
            self.set_tool("select")
        elif e.GetKeyCode() == wx.WXK_F9:
            self.toggle_side()
        else:
            e.Skip()

    def toggle_side(self, show=None):
        """Show or hide the side panel (top bar button, F9 or right click).
        The panel is secondary: everything is also in the right click menus."""
        show = (not self.side.IsShown()) if show is None else bool(show)
        if hasattr(self.frame, "set_side"):
            self.frame.set_side(show)
        else:
            self.side.Show(show)
            self.Layout()

    def view_menu(self):
        """Last entries of every plot menu: the side panel (tiles: see tile_menu)."""
        return [(None, None), ("Hide the side panel (F9)" if self.side.IsShown() else "Show the side panel (F9)",
                               lambda: self.toggle_side())]

    def set_tool(self, key):
        self.tool = key
        self.status(self.TOOL_HINTS.get(key, ""))
        self.tool_changed()

    def tool_changed(self):
        pass

    def tool_action(self, key):
        if key == "panel":
            self.toggle_side()
            return
        if key == "deconvolute" and hasattr(self, "on_deconvolve"):
            self.on_deconvolve()
            return
        if key == "full":
            for card in self.all_cards():
                card.reset_view()
        elif key == "auto":
            self.on_integrate()
        elif key == "clear":
            self.on_clear_peaks()

    def all_cards(self):
        return [v.card for v in self.views]

    def _trace_near(self, v, x, y):
        if not v.traces:
            return None
        if len(v.traces) == 1 or y is None:
            return v.traces[0]
        best, dist = v.traces[0], None
        for tr in v.traces:
            yd = float(np.interp(x, tr["t"], tr["y"] * tr.get("_scale", 1.0)))
            d = abs(yd - y)
            if dist is None or d < dist:
                best, dist = tr, d
        return best

    def _peak_at(self, v, x):
        """Index of an integrated peak of this view that contains time x."""
        names = [tr["name"] for tr in v.traces]
        hits = [i for i, p in enumerate(self.peaks) if p.get("trace") in names and p["t0"] <= x <= p["t1"]]
        if not hits:
            return None
        return min(hits, key=lambda i: abs(self.peaks[i]["rt"] - x))

    def _add_peaks(self, new, trace, tr=None):
        if tr is None:
            tr = self._find_trace(trace)[1]
        for p in new:
            p["trace"] = trace
            if tr is not None:
                p["key"] = tr.get("key", trace)
                p["sig"] = _trace_sig(tr)
        self.peaks = LI.finish(self.peaks + list(new))
        self.sel_peak = self.peaks.index(new[0]) if new else None
        self.refresh_peaks()
        if self.sel_peak is not None:
            self.table.select(self.sel_peak)

    def _show_integrated(self, trace, p, quiet=False):
        """A peak of this trace with the apex of p is integrated already (the
        same peak twice would halve its area %): select it and return True."""
        same = [i for i, q in enumerate(self.peaks) if q.get("trace") == trace and abs(q["rt"] - p["rt"]) < 1e-9]
        if not same:
            return False
        self.select_peak(same[0])
        self.table.select(same[0])
        if not quiet:
            self.status("The peak at %.3f min is integrated already; delete it first" % p["rt"])
        return True

    def chrom_tool(self, v, tool, kind, a, b):
        if tool == "calibrate":
            if kind == "drag" and hasattr(self, "calibrate_range"):
                self.calibrate_range(a, b)
            return
        if not v.traces:
            return
        if tool == "drag_peak" and kind == "drag":
            selected = self.int_trace.GetStringSelection()
            tr = next((t for t in v.traces if t["name"] == selected), None)
            if tr is None:
                tr = self._trace_near(v, (a + b) / 2.0, None)
            p = LI.manual_peak(tr["t"], tr["y"], a, b)
            if p is None:
                self.status("The dragged range is too short")
                return
            if self._show_integrated(tr["name"], p):
                return
            self._add_peaks([p], tr["name"], tr)
            self.status("%s: peak at %.3f min, area %s" % (tr["name"], p["rt"], _g(p["area"])))
        elif tool == "click_peak":
            tr = self._trace_near(v, a, b)
            p = LI.peak_at(tr["t"], tr["y"], a)
            if p is None:
                self.status("No peak found here; use Drag instead")
                return
            if self._show_integrated(tr["name"], p, quiet=True):
                return
            self._add_peaks([p], tr["name"], tr)
            self.status("Peak at %.3f min (%.3f to %.3f min): area %s" % (p["rt"], p["t0"], p["t1"], _g(p["area"])))
        elif tool in ("split", "delete"):
            i = self._peak_at(v, a)
            if i is None:
                self.status("Click inside an integrated peak")
                return
            if tool == "delete":
                self.delete_peak(i)
                return
            p = self.peaks[i]
            _, tr = self._find_trace(p["trace"])
            parts = LI.split_peak(tr["t"], tr["y"], p, a) if tr is not None else None
            if not parts:
                self.status("Click further inside the peak to split it")
                return
            del self.peaks[i]
            self._add_peaks(list(parts), p["trace"])

    # ------------------------------------------------------- chromatograms
    def normalise(self):
        return bool(getattr(self, "norm", None) and self.norm.GetValue())

    def chrom_extras(self, ax, scale_of, view):
        return []  # extra legend entries [(colour, label)]

    def empty_text(self, view):
        return "Open a data file"

    def plot_view(self, v, keep_view=False):
        card = v.card
        ax = card.ax
        view = ax.get_xlim() if (keep_view and card.full) else None
        card.reset()
        norm = self.normalise() and len(v.traces) > 1
        style_axes(ax, "Time (min)", "Relative intensity (%)" if norm else self.ylabel)
        if not v.traces:
            ax.text(0.5, 0.5, self.empty_text(v), transform=ax.transAxes, ha="center", va="center",
                    fontsize=9, color=MUTED)
            card.full = None
            card.set_legend([])
            card.draw()
            return
        scale_of = {}
        lo, hi = None, None
        for tr in v.traces:
            y = tr["y"]
            top = float(np.nanmax(y)) if len(y) else 0.0
            s = 100.0 / top if (norm and top > 0) else 1.0
            scale_of[tr["name"]] = s
            tr["_line"], = ax.plot(tr["t"], y * s, color=tr["color"], lw=LW_TRACE, label=tr["label"], zorder=3)
            tr["_scale"] = s
            if len(tr["t"]):
                lo = tr["t"][0] if lo is None else min(lo, tr["t"][0])
                hi = tr["t"][-1] if hi is None else max(hi, tr["t"][-1])
            mine = [i for i, p in enumerate(self.peaks) if p.get("trace") == tr["name"]]
            if mine:
                sp = [dict(self.peaks[i], b0=self.peaks[i]["b0"] * s, b1=self.peaks[i]["b1"] * s,
                           apex_y=self.peaks[i]["apex_y"] * s) for i in mine]
                sel = mine.index(self.sel_peak) if self.sel_peak in mine else None
                draw_peaks(ax, tr["t"], y * s, sp, tr["color"], sel)
        if self.avg and self.avg[1] != self.avg[0]:
            ax.axvspan(self.avg[0], self.avg[1], color=AVG_COL, alpha=0.13, lw=0, zorder=0)._range_mark = True
        if self.bg:
            ax.axvspan(self.bg[0], self.bg[1], facecolor=BG_COL, alpha=0.16, lw=0, zorder=0, hatch="///",
                       edgecolor="#FFFFFF")._range_mark = True
        extra = self.chrom_extras(ax, scale_of, v) or []
        card.set_legend([(tr["color"], tr["label"]) for tr in v.traces] + list(extra))
        if card.pair is not None and lo is not None:  # the plots of a tile share the time axis: the full view
            for o in self.views:
                if o is not v and o.card.pair is card.pair:
                    for tr in o.traces:
                        if len(tr["t"]):
                            lo, hi = min(lo, tr["t"][0]), max(hi, tr["t"][-1])
        card.full = (lo, hi)
        ax.set_xlim(*(view if view else card.full))
        card.autoscale_y()
        card.draw()

    def plot_chroms(self, keep_view=False):
        for v in self.views:
            self.plot_view(v, keep_view)

    def busy_averaging(self):
        """True (and says so) while a spectrum is still being averaged: a new
        range then changes nothing, so the shading, the fields, the linked view
        and undo stay with the spectra shown."""
        if getattr(self, "_averaging", False):
            self.status("Still averaging the previous range")
            return True
        return False

    def chrom_range(self, a, b, shift):
        if self.busy_averaging():
            self.plot_chroms(keep_view=True)  # the dragged span goes, the range shown stays
            return
        if shift:
            self.bg = (a, b)
        else:
            self.avg = (a, b)
        self.fill_range_fields()
        self.plot_chroms(keep_view=True)
        self.range_changed(shift)

    def chrom_pick(self, x, y, view):
        pass

    def chrom_readout(self, x, y, view):
        return "%.3f min   %s" % (x, _g(y))

    def chrom_view_changed(self, v):
        if self.link is not None and self.link.GetValue():
            xl = v.card.ax.get_xlim()
            for o in self.views:
                if o is not v and o.card.full:
                    o.card.ax.set_xlim(*xl)
                    o.card.autoscale_y()
                    o.card.draw()

    def chrom_menu(self, x, y, v):
        items = [("Average a time range\u2026", self.ask_average),
                 ("Integrate the peaks in the visible window\u2026", lambda: self.ask_integrate(v))]
        if self.avg and self.avg[1] != self.avg[0]:
            items.append(("Add the selected range as a peak", self.on_add_peak))
        if self.peaks:
            items += [("Export the peak table\u2026", self.on_export_table),
                      ("Remove every integrated peak", self.on_clear_peaks)]
        items += self.extra_chrom_menu(x, y, v)
        items += [(None, None), ("Full view", v.card.reset_view),
                  ("Copy chromatogram data", lambda: self.copy_chrom(v)),
                  ("Export chromatogram\u2026", lambda: self.export_chrom(v))]
        return items + self.view_menu()

    def extra_chrom_menu(self, x, y, v):
        return []

    def ask_average(self):
        """Time range to average and background (small window)."""
        if self.busy_averaging():
            return
        fields = [dict(key="rng", label="Average", kind="pair", value=(self.t0.GetValue(), self.t1.GetValue()),
                       unit="min"),
                  dict(key="bg", label="Background", kind="pair", value=(self.b0.GetValue(), self.b1.GetValue()),
                       unit="min", tip="Empty = no background"),
                  dict(key="use_bg", label="Subtract the background", kind="check", value=self.use_bg.GetValue())]
        if getattr(self, "binw", None) is not None and self.binw.IsShown():
            fields.append(dict(key="binw", label="Bin width", kind="text", value=self.binw.GetValue(), unit="m/z"))
        v = ask_form(self, "Average a time range", fields, ok="Average")
        if v is None:
            return
        self.t0.SetValue(v["rng"][0])
        self.t1.SetValue(v["rng"][1])
        self.b0.SetValue(v["bg"][0])
        self.b1.SetValue(v["bg"][1])
        self.use_bg.SetValue(v["use_bg"])
        if "binw" in v:
            self.binw.SetValue(v["binw"])
        self.average_fields()

    def average_fields(self):
        """Average the range typed in the fields (overridden per view)."""
        pass

    def ask_integrate(self, v=None):
        names = [tr["name"] for tr in self.all_traces()]
        if not names:
            self.status("Nothing to integrate yet")
            return
        cur = self.int_trace.GetStringSelection()
        if v is not None and v.traces and cur not in [tr["name"] for tr in v.traces]:
            cur = v.traces[0]["name"]
        vals = ask_form(self, "Integrate peaks", [
            dict(key="trace", label="Trace", kind="choice", choices=names,
                 value=names.index(cur) if cur in names else 0, width=190),
            dict(key="thr", label="Minimum height", kind="text", value=self.int_thr.GetValue(), unit="% of top"),
            dict(key="w", label="Minimum width", kind="text", value=self.int_w.GetValue(), unit="s")],
            note="Peaks of this trace in the visible window are replaced.", ok="Integrate", colour=T.GROUP["red"])
        if vals is None:
            return
        self.int_trace.SetSelection(vals["trace"])
        self.int_thr.SetValue(vals["thr"])
        self.int_w.SetValue(vals["w"])
        self.on_integrate()

    def _trace_table(self, v):
        if not v.traces:
            return None, None
        t = v.traces[0]["t"]
        cols = [t] + [tr["y"] if len(tr["y"]) == len(t) else np.interp(t, tr["t"], tr["y"]) for tr in v.traces]
        head = ["Time (min)"] + [tr["name"] for tr in v.traces]
        return head, np.column_stack(cols)

    def copy_chrom(self, v):
        head, data = self._trace_table(v)
        if data is None:
            return
        txt = "\t".join(head) + "\n" + "\n".join("\t".join("%.6g" % x for x in row) for row in data)
        if _clip(txt):
            self.status("Chromatogram data copied (%d rows)" % len(data))
        else:
            self.status("The clipboard is in use by another program; try again")

    def export_chrom(self, v):
        head, data = self._trace_table(v)
        if data is None:
            return
        tag = getattr(v, "file_tag", "")
        path = self.frame.ask_save("Export chromatogram", self.frame.base_name() + "_chromatogram%s.csv" % tag,
                                   "CSV (*.csv)|*.csv|Text (*.txt)|*.txt")
        if not path:
            return
        sep = "," if path.lower().endswith(".csv") else "\t"
        try:
            with open(path, "w", encoding="utf-8") as fh:
                fh.write(sep.join(head) + "\n")
                for row in data:
                    fh.write(sep.join("%.6g" % x for x in row) + "\n")
        except Exception as ex:
            _report_error(self, _write_failed(path), ex)
            return
        self.status("Saved " + path)

    # --------------------------------------------------------- integration
    def build_integration(self):
        sp = self.side
        sp.section("Peak integration", T.GROUP["red"])
        self.int_trace = sp.choice([])
        self.int_trace.SetMinSize(wx.Size(self.FromDIP(160), -1))
        sp.row("Trace", self.int_trace)
        self.int_thr = sp.text("2")
        sp.row("Minimum height", self.int_thr, "% of top")
        self.int_w = sp.text("2")
        sp.row("Minimum width", self.int_w, "s")
        sp.buttons(flat(sp, "Integrate", "primary", icon="integrate", handler=self.on_integrate),
                   flat(sp, "Add range", handler=self.on_add_peak,
                        tooltip="Integrate the selected range as one peak"))
        sp.buttons(flat(sp, "Clear", handler=self.on_clear_peaks),
                   flat(sp, "Export table", icon="export", handler=self.on_export_table))

    def integration_get(self):
        """Integration settings of a method preset (method_presets.py)."""
        return {"integration_min_height": _num(self.int_thr, None), "integration_min_width": _num(self.int_w, None)}

    def integration_set(self, v):
        for key, ctrl in (("integration_min_height", self.int_thr), ("integration_min_width", self.int_w)):
            if key in v:
                ctrl.ChangeValue("%g" % v[key])

    def update_trace_choice(self):
        names = [tr["name"] for tr in self.all_traces()]
        cur = self.int_trace.GetStringSelection()
        self.int_trace.Set(names)
        if names:
            self.int_trace.SetSelection(names.index(cur) if cur in names else 0)
        # integrated peaks follow their trace: when it changes (other
        # wavelength of the slider, smoothing, TIC or base peak, calibration)
        # each peak is integrated again between the same limits
        by_key, by_name = {}, {}
        for tr in self.all_traces():
            by_key[tr.get("key", tr["name"])] = tr
            by_name[tr["name"]] = tr
        kept, changed = [], False
        for p in self.peaks:
            tr = by_key.get(p.get("key")) or by_name.get(p.get("trace"))
            if tr is None:
                changed = True
                continue
            sig = _trace_sig(tr)
            if p.get("sig") == sig and p.get("trace") == tr["name"]:
                kept.append(p)
                continue
            q = LI.manual_peak(tr["t"], tr["y"], p["t0"], p["t1"])
            changed = True
            if q is None:
                continue
            q.update(trace=tr["name"], key=tr.get("key", tr["name"]), sig=sig)
            kept.append(q)
        if changed:
            self.peaks = LI.finish(kept)
            if self.sel_peak is not None and self.sel_peak >= len(self.peaks):
                self.sel_peak = None
            self.table.set_peaks(self.peaks, self.units)

    def _find_trace(self, name):
        for v in self.views:
            for tr in v.traces:
                if tr["name"] == name:
                    return v, tr
        return None, None

    def _int_source(self):
        v, tr = self._find_trace(self.int_trace.GetStringSelection())
        if tr is None and self.views and self.views[0].traces:
            v, tr = self.views[0], self.views[0].traces[0]
        return v, tr

    def integrate_view(self, v):
        if v.traces:
            names = [tr["name"] for tr in self.all_traces()]
            cur = self.int_trace.GetStringSelection()
            if cur not in [tr["name"] for tr in v.traces]:
                self.int_trace.SetSelection(names.index(v.traces[0]["name"]))
        self.on_integrate()

    def on_integrate(self, e=None):
        v, tr = self._int_source()
        if tr is None:
            self.status("Nothing to integrate yet")
            return
        x0, x1 = v.card.ax.get_xlim()
        new = LI.auto_integrate(tr["t"], tr["y"], _num(self.int_thr, 2.0), _num(self.int_w, 2.0),
                                smooth_points=0, t_range=(x0, x1))
        for p in new:
            p["trace"] = tr["name"]
            p["key"] = tr.get("key", tr["name"])
            p["sig"] = _trace_sig(tr)
        # replaces this trace's peaks inside the window; other peaks stay
        keep = [p for p in self.peaks if p.get("trace") != tr["name"] or p["t1"] < x0 or p["t0"] > x1]
        self.peaks = LI.finish(keep + new)
        self.sel_peak = None
        self.refresh_peaks()
        self.status("%d peaks integrated in %s (%.2f to %.2f min)" % (len(new), tr["name"], x0, x1))

    def on_add_peak(self, e=None):
        v, tr = self._int_source()
        if tr is None:
            return
        if not self.avg or self.avg[1] == self.avg[0]:
            self.status("Select a time range on the chromatogram first")
            return
        p = LI.manual_peak(tr["t"], tr["y"], *self.avg)
        if p is None:
            self.status("The selected range is too short")
            return
        if self._show_integrated(tr["name"], p):
            return
        self._add_peaks([p], tr["name"], tr)

    def on_clear_peaks(self, e=None):
        self.peaks, self.sel_peak = [], None
        self.refresh_peaks()

    def delete_peak(self, i):
        if 0 <= i < len(self.peaks):
            del self.peaks[i]
            self.peaks = LI.finish(self.peaks)
            self.sel_peak = None
            self.refresh_peaks()

    def refresh_peaks(self):
        self.table.set_peaks(self.peaks, self.units)
        self.plot_chroms(keep_view=True)

    def select_peak(self, i):
        if not (0 <= i < len(self.peaks)):
            return
        self.sel_peak = i
        self.plot_chroms(keep_view=True)
        self.peak_selected(self.peaks[i])

    def peak_selected(self, p):
        pass

    def on_export_table(self, e=None):
        if not self.peaks:
            self.status("No peaks to export")
            return
        path = self.frame.ask_save("Export peak table", self.frame.base_name() + "_peaks.csv", "CSV (*.csv)|*.csv")
        if path:
            try:
                LI.write_table(path, self.peaks, self.units or "counts")
            except Exception as ex:
                _report_error(self, _write_failed(path), ex)
                return
            self.status("Saved " + path)

    # ------------------------------------------------------ range fields
    def fill_range_fields(self):
        if self.avg and self.avg[1] != self.avg[0]:
            self.t0.SetValue(_fmt_t(self.avg[0]))
            self.t1.SetValue(_fmt_t(self.avg[1]))
        if self.bg:
            self.b0.SetValue(_fmt_t(self.bg[0]))
            self.b1.SetValue(_fmt_t(self.bg[1]))

    def read_range_fields(self):
        """Range and background typed in the fields. One time typed (from or
        to) is that one time, none is no range. Returns a note for the status
        bar when the background has only one time (then none is subtracted)."""
        a, b = _num(self.t0), _num(self.t1)
        if a is not None and b is not None:
            self.avg = (min(a, b), max(a, b))
        elif a is not None or b is not None:
            t = a if a is not None else b
            self.avg = (t, t)
        else:
            self.avg = None
        c, d = _num(self.b0), _num(self.b1)
        self.bg = (min(c, d), max(c, d)) if (c is not None and d is not None) else None
        if (c is None) != (d is None):
            return "The background needs a start and an end time: no background subtracted"
        return None

    def range_changed(self, shift):
        pass


def _trace_sig(tr):
    """Identity of a trace's data (a changed signature means the peak must be
    integrated again)."""
    y = tr.get("y")
    try:
        return (tr.get("name"), len(y), round(float(np.nansum(y)), 6))
    except Exception:
        return (tr.get("name"), 0, 0.0)


def _clip(text):
    """Text to the clipboard; False when another program holds it."""
    for attempt in range(10):  # another program may hold the clipboard for a moment
        if wx.TheClipboard.Open():
            try:
                ok = wx.TheClipboard.SetData(wx.TextDataObject(text))
                if ok:
                    wx.TheClipboard.Flush()  # stays after the program is closed
            finally:
                wx.TheClipboard.Close()
            return bool(ok)
        time.sleep(0.05)
    return False


def _write_failed(path):
    """First line of the message when a file could not be written."""
    return ("%s could not be saved in %s: it may be open in another program"
            % (os.path.basename(path), os.path.dirname(path) or "this folder"))


def _mouse_notes(sp, lines):
    return None


def _pol_shared(events, e):
    """True when another scan event has the polarity of event e (scan + SIM of one polarity)."""
    try:
        pol = events[e].get("polarity")
        return pol in ("+", "-") and sum(1 for ev in events if ev.get("polarity") == pol) > 1
    except Exception:
        return False


def _pol_short(pol, e, events=None):
    """'(+)', '(\u2212)'; '(+, event 2)' when another event has the same polarity
    (trace names, file names and peaks must tell the events apart)."""
    if events is not None and _pol_shared(events, e):
        return {"+": "(+, event %d)", "-": "(\u2212, event %d)"}[pol] % (e + 1)
    return {"+": "(+)", "-": "(\u2212)"}.get(pol, "(event %d)" % (e + 1))


def _pol_name(pol, e):
    return {"+": "Positive ions", "-": "Negative ions"}.get(pol, "Scan event %d" % (e + 1))


# ==========================================================================
# spectrum tiles: drawing, pinned labels, measuring, export. Shared by the
# Mass spectrometry view (MSTab) and the m/z tiles of the Compare view
# (lcms_compare.CompareTab). A spectrum is described by a dict "col":
# {"spec" (m/z, intensity) or None, "spec_card" (PlotCard), "pins" (m/z
# pinned), "measure" (m1, m2 or None, text), "m1" (first peak measured)}.
# ==========================================================================
SPEC_WILDCARD = "Text, m/z and intensity (*.txt)|*.txt|JCAMP-DX (*.jdx)|*.jdx"


def shared_mz_range(data, cols):
    """The m/z range of the full view of the spectra of cols (scan events of data): the acquired range,
    the same for the plots of one tile (they share the m/z axis)."""
    lo = hi = None
    for c in cols:
        e = c["e"]
        ev = data.events[e] if data is not None and e is not None and e < len(data.events) else {}
        sp = c.get("spec")
        a = ev.get("mz_low") or (float(sp[:, 0].min()) if sp is not None and len(sp) else None)
        b = ev.get("mz_high") or (float(sp[:, 0].max()) if sp is not None and len(sp) else None)
        if a is None or b is None:
            continue
        lo, hi = (a, b) if lo is None else (min(lo, a), max(hi, b))
    if lo is None:
        return 0.0, 1.0
    return (lo, hi) if hi > lo else (lo - 1.0, hi + 1.0)


def ask_mz_range(win, card):
    """A small window for the m/z range of a spectrum (the other spectra of its tile follow)."""
    if not card.full:
        return
    x0, x1 = card.ax.get_xlim()
    v = ask_form(win, "m/z range", [dict(key="lo", label="From", value="%.6g" % x0, unit="m/z"),
                                    dict(key="hi", label="To", value="%.6g" % x1, unit="m/z")], ok="Show")
    if v is None:
        return
    a, b = _num_s(v["lo"], None), _num_s(v["hi"], None)
    if a is None or b is None or a == b:
        win.status("Enter the lower and the upper m/z")
        return
    set_mz_range(card, min(a, b), max(a, b), lock=True)


def _tile_cards(card):
    pair = getattr(card, "pair", None)
    return list(pair.cards) if pair is not None and getattr(pair, "cards", None) else [card]


def set_mz_range(card, a, b, lock=False):
    """The m/z range of a spectrum; lock: kept (double click and new spectra come back to it, "Full m/z range"
    releases it), for every plot of the tile."""
    if lock:
        for c in _tile_cards(card):
            c.xlock = (float(a), float(b))
    card.ax.set_xlim(a, b)
    card.ylock = None
    card.autoscale_y()
    card.draw()
    card._view_changed()


def mz_range_items(win, card):
    """Menu entries of a spectrum: m/z range (a small window, kept until the full m/z range is chosen) and the full
    m/z range."""
    def full():
        for c in _tile_cards(card):
            c.xlock = None
        card.reset_view()
    return [("m/z range\u2026", lambda: ask_mz_range(win, card)), ("Full m/z range", full)]


def spec_nearest_peak(spec, xlim, x, pick_min):
    """m/z of the strongest point near x (within 1/150 of the view, at least pick_min), or None."""
    if spec is None or not len(spec):
        return None
    mz, it = spec[:, 0], spec[:, 1]
    x0, x1 = xlim
    w = max((x1 - x0) / 150.0, pick_min)
    m = np.where(np.abs(mz - x) <= w)[0]
    if not len(m):
        return None
    return float(mz[m[np.argmax(it[m])]])


def spec_peak_height(spec, m, pick_min):
    w = np.abs(spec[:, 0] - m) <= pick_min
    return float(spec[w, 1].max()) if np.any(w) else 0.0


def measure_text(m1, m2, sign):
    """Distance of two peaks; isotope spacing (charge) or adjacent charge states (charges and mass) for
    ions of the polarity sign (+1 or -1)."""
    H = 1.007276
    lo, hi = sorted((m1, m2))
    d = hi - lo
    if d <= 0:
        return "Same peak"
    if d < 1.3:
        z = int(round(1.0 / d))
        return "Δ %.3f m/z (isotope spacing: z ≈ %d)" % (d, z) if z <= 60 else "Δ %.3f m/z" % d
    zf = (lo - sign * H) / d
    z = int(round(zf))
    if z >= 1 and abs(zf - z) < 0.2:
        mass = z * (hi - sign * H)
        ch = "+" if sign > 0 else "−"
        return "Δ %.2f m/z: charges %d%s and %d%s, M ≈ %.1f Da" % (d, z + 1, ch, z, ch, mass)
    return "Δ %.2f m/z" % d


def draw_spectrum(card, spec, full, sticks, view=None, n_labels=8, fmt="%.1f", refine=None, color=TRACE):
    """The spectrum as a profile line or as sticks, over the full view (x0, x1) or the view kept, the y
    range made automatically and the n strongest maxima labelled; returns the labels."""
    ax = card.ax
    mz, it = spec[:, 0], spec[:, 1]
    if sticks:
        from ms_plot_render import SpectrumSticks
        artist = SpectrumSticks(mz, it, color=color, linewidth=LW_STICK)
        ax.add_collection(artist, autolim=False)
        card.extra_scale = [artist.full_data()]
    else:
        plot_trace(ax, mz, it, color=color, lw=LW_SPEC)
    card.full = full
    ax.set_xlim(*(view if view else (card.xlock or card.full)))
    card.autoscale_y(pad=0.12)
    x0, x1 = ax.get_xlim()
    if n_labels <= 0:
        return []
    return label_maxima(ax, mz, it, n=n_labels, fmt=fmt, min_sep=(x1 - x0) / 24.0, xlim=(x0, x1), refine=refine)


def _clear_of(ax, ann, placed):
    """A label at a point placed where it does not cover the labels placed before: right of the point, else left
    of it, else lower down on the right. Returns placed with it."""
    try:
        r = ax.figure.canvas.get_renderer()
        x0, y0 = ann.xyann
        tries = [(x0, y0, "left"), (-x0, y0, "right")]
        bb = ann.get_window_extent(r)
        h = bb.height * 72.0 / ax.figure.dpi + 2.0
        tries += [(x0, y0 - k * h, "left") for k in range(1, 5)]
        for x, y, ha in tries:
            ann.xyann = (x, y)
            ann.set_ha(ha)
            bb = ann.get_window_extent(r)
            if not any(bb.overlaps(p) for p in placed):
                break
        else:
            ann.set_ha("left")
            ann.xyann = (x0, y0)
        return placed + [ann.get_window_extent(r)]
    except Exception:
        return placed


def draw_charge_marks(ax, spec, sm):
    """The charge states of a found mass (Find mass) over the spectrum: a triangle in the colour of the mass and
    "13+" over the states present, a grey triangle over the others. sm: M, c (mass per charge), w, zs, on,
    colour, sign."""
    from matplotlib.transforms import offset_copy
    mz, it = spec[:, 0], spec[:, 1]
    x0, x1 = ax.get_xlim()
    tm = offset_copy(ax.transData, fig=ax.figure, y=3.5, units="points")
    tl = offset_copy(ax.transData, fig=ax.figure, y=7.5, units="points")
    on = set(sm.get("on") or [])
    sign = "+" if sm.get("sign", 1) > 0 else "−"
    for z in sm.get("zs") or []:
        m = (sm["M"] + z * sm["c"]) / z
        if not (min(x0, x1) <= m <= max(x0, x1)):
            continue
        sel = np.abs(mz - m) <= max(0.3, sm.get("w", 5.0) / z)
        yv = float(it[sel].max()) if np.any(sel) else 0.0
        col = sm["colour"] if z in on else ABSENT_COL
        ln, = ax.plot([m], [yv], ls="none", marker="v", ms=5, mew=0, color=col, transform=tm, zorder=9,
                      clip_on=False)
        ln._no_scale = True
        if z in on:
            ax.text(m, yv, "%d%s" % (z, sign), transform=tl, fontsize=7.5, color=col, ha="center", va="bottom",
                    zorder=9, clip_on=False)


def draw_spec_marks(ax, col, auto, pin_fmt, height):
    """Pinned labels (an automatic label next to one is removed) and the measurement of a spectrum;
    height(m): top of the peak at m."""
    x0, x1 = ax.get_xlim()
    pins = [m for m in col.get("pins", []) if x0 <= m <= x1]
    for a in auto:
        if any(abs(a.xy[0] - m) < (x1 - x0) / 60.0 for m in pins):
            a.remove()
    for m in pins:
        ax.annotate(pin_fmt % m, (m, height(m)), xytext=(0, 4), textcoords="offset points",
                    ha="center", va="bottom", fontsize=7.5, color=C["accent_text"], fontweight="bold", zorder=5,
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=C["accent"], lw=0.6, alpha=0.95))
    meas = col.get("measure")
    if meas:
        m1, m2, text = meas
        if m2 is None:
            ln = ax.axvline(m1, color=C["accent"], lw=0.9, ls=(0, (3, 2)))
            ln._no_scale = True
        else:
            h = max(height(m1), height(m2))
            ytop = h * 1.04
            for m in (m1, m2):
                tk, = ax.plot([m, m], [height(m), ytop], color=C["accent"], lw=0.8,
                              ls=(0, (2, 2)))
                tk._no_scale = True
            br, = ax.plot([m1, m2], [ytop, ytop], color=C["accent"], lw=1.2)
            br._no_scale = True
            ax.annotate(text, ((m1 + m2) / 2.0, ytop), xytext=(0, 3), textcoords="offset points", ha="center",
                        va="bottom", fontsize=7.5, color=C["accent_text"], zorder=6,
                        bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.9))


def spec_click(col, tool, pk, pin_fmt, measure, status):
    """Label (pin or unpin the label of the peak at m/z pk) or Measure (first peak, then the second:
    measure(m1, m2) gives the text) on a spectrum; True when it is to be drawn again."""
    if tool == "label":
        pins = col.setdefault("pins", [])
        near = [m for m in pins if abs(m - pk) < 1e-6]
        if near:
            pins.remove(near[0])
        else:
            pins.append(pk)
        return True
    if tool == "measure":
        if col.get("m1") is None:
            col["m1"] = pk
            col["measure"] = (pk, None, "")
            status(("First peak m/z " + pin_fmt + "; now click the second peak") % pk)
        else:
            m1 = col.pop("m1")
            text = measure(m1, pk)
            col["measure"] = (m1, pk, text)
            status(text)
        return True
    return False


def spec_pin_item(col, pk, pin_fmt, label):
    """Menu entry Label m/z (or Remove the label of m/z) of the peak at pk; label(m) pins or unpins it."""
    pinned = any(abs(m - pk) < 1e-6 for m in col.get("pins", []))
    return ((("Remove the label of m/z " if pinned else "Label m/z ") + pin_fmt) % pk, lambda m=pk: label(m))


def spec_clear_items(col, replot):
    """Menu entries that remove every pinned label and the measurement of a spectrum."""
    items = []
    if col.get("pins"):
        items.append(("Remove every label", lambda: (col.update(pins=[]), replot())))
    if col.get("measure"):
        items.append(("Remove the measurement", lambda: (col.update(measure=None, m1=None), replot())))
    return items


def spec_display_item(sticks, set_style):
    """Menu entry Display: profile (line) or sticks; set_style(0 or 1)."""
    return ("Display", [("Profile (line)", lambda: set_style(0), not sticks),
                        ("Sticks", lambda: set_style(1), sticks)])


def copy_spectrum(win, spec):
    """The spectrum to the clipboard as text (m/z and intensity, tab separated)."""
    if spec is None:
        return
    if _clip("\n".join("%.4f\t%.2f" % (a, b) for a, b in spec)):
        win.status("Spectrum copied (%d points)" % len(spec))
    else:
        win.status("The clipboard is in use by another program; try again")


def write_spectrum_file(path, spec, name, pol, rng):
    """A spectrum as text (m/z and intensity) or, for .jdx and .dx, as JCAMP-DX with its polarity and time
    range."""
    if path.lower().endswith((".jdx", ".dx")):
        t0, t1 = rng if rng else (None, None)
        lcms_data.write_spectrum_jdx(path, spec, name, pol or "+", t0, t1)
    else:
        lcms_data.write_spectrum_txt(path, spec)


def open_spectrum_in_deconvolute(win, owner, spec, name, sign, rng):
    """A spectrum written next to the data (owner: the file, with out_dir and write_spectrum) and opened in
    the full Deconvolute window (separate); returns the path written or None."""
    import ms_deconv
    # UniDec writes <name>_unidecfiles\\<name>_...: a name short enough for
    # the Windows path limit (long file names in deep folders)
    folder, base = ms_deconv._unidec_place(owner.out_dir(), name)
    path = os.path.join(folder, base + ".txt")
    try:
        os.makedirs(folder, exist_ok=True)
        if sign < 0 and hasattr(lcms_data, "write_spectrum_jdx"):
            # negative ions: JCAMP-DX carries the polarity, so the Deconvolute
            # window starts in negative mode (loss of H+ as the charge carrier)
            path = os.path.join(folder, base + ".jdx")
            rng = rng or (None, None)
            lcms_data.write_spectrum_jdx(path, spec, name, "-", rng[0], rng[1])
        else:
            owner.write_spectrum(path, spec)
    except Exception as ex:
        _report_error(win, _write_failed(path), ex)
        return None
    open_in_unidec(path)
    win.status("Opening in the Deconvolute window: " + path)
    return path


# ==========================================================================
# MS tab: the tiles of the scan events one below the other (positive, negative)
# ==========================================================================
class MSTab(TabBase):
    ylabel = "Intensity"
    units = ""
    MAXCOLS = 4
    MZ_FMT = "%.1f"  # automatic peak labels
    PIN_FMT = "%.2f"  # pinned labels and messages
    APEX_LABELS = False  # peak labels at the apex of profile peaks (HRMS) instead of the highest point
    apex_fn = staticmethod(apex_of_profile)  # apex of a labelled profile peak (HRMS: hrms_calib.profile_apex)
    READ_FMT = "m/z %.2f   %s"  # hover readout
    PICK_MIN = 0.3  # smallest m/z window for clicking a peak
    TOL_DEFAULT = "0.5"
    FULL_WINDOW = True  # offer the separate Deconvolute window (UniDec GUI) for this kind of data
    OPEN_HINT = "Open a data file or drop it here"
    FIND_MASS = True  # Find mass (right click a chromatogram): the mass search tiles

    def __init__(self, parent, frame):
        TabBase.__init__(self, parent, frame)
        self.data = None
        self.cols = []
        self.xics = {}
        self.xic_hidden = set()  # (event, xic_key) of mass chromatograms hidden in the files panel tree
        self.xic_views = []
        self.xic_rows = []
        self._xic_sig = None
        # mass searches (Find mass): [{"id", "masses", "ws", "names" (a row each), "carrier", "mzmin", "events",
        # "res": {e: [result of each row]}}] (search_rows)
        self.searches = []
        self.search_views = []
        self.search_rows = []
        self._search_sig = None
        self._search_job = None
        self.pick_t = None
        self.scroller = RowScroller(self)
        self.rows = StackBox(self.scroller)
        self.rows.on_close_tile = self._close_search_tile
        self.rows.shape = tile_shape
        if self.FIND_MASS:  # LCMS Analysis (HRMS Analysis keeps its own tile heights)
            # by default every plot tile (chromatograms, spectra, deconvolution, mass searches) fills the window,
            # the tables below: the spectrum with the fit and charge states of a deconvolved mass stays in view
            # next to the deconvolved spectrum (with + and - one below the other, since 4.51, it was scrolled out)
            self.rows.fit_share = self._fit_share
            self.rows.fit_label = "Plots fill the window"
            self.rows.fit_floor = 60  # (plots lower than this: a very small window scrolls instead)
            self.rows.mode_key, self.rows.default_mode = "lcms_tile_mode", "two"
        self.scroller.set_box(self.rows)
        self.bottom = SplitBox(self.rows, wx.HORIZONTAL)
        self.table = PeakTable(self.bottom, self.select_peak, self.delete_peak)
        self.bottom.set_panes([self.table], [1])
        self.extras = {}  # key -> (row or None, table or None, weight); deconvolution, calibration
        self.panels = []  # objects with overlay(col, ax) drawn onto the spectra

        sp = self.side
        sp.section("Chromatograms", T.GROUP["blue"])
        self.kind = sp.choice(["TIC (total ion current)", "BPC (base peak)", "Mass chromatograms only"])
        self.kind.SetMinSize(wx.Size(self.FromDIP(160), -1))
        sp.row("Trace", self.kind)
        self.xic_fields, self.xic_labels = [], []
        for i in range(2):
            t = sp.text("", 150, tooltip="m/z values, separated by commas (803.2+-0.3, 889.09+-10ppm)")
            srow = sp.row("m/z list %s" % ("(+)" if i == 0 else "(\u2212)"), t)
            self.xic_fields.append(t)
            self.xic_labels.append(srow.GetChildren()[0].GetWindow())
        self.tol = sp.text(self.TOL_DEFAULT)
        self.tol_unit = sp.choice(["m/z", "ppm"])
        self.tol_unit.SetMinSize(wx.Size(self.FromDIP(64), -1))
        sp.row("Default window \u00b1", self.tol, self.tol_unit)
        self.xic_mode = sp.choice(["Stacked tiles", "Overlaid on the TIC"])
        self.xic_mode.SetMinSize(wx.Size(self.FromDIP(160), -1))
        sp.row("Mass chrom. shown", self.xic_mode)
        self.smooth = sp.spin(0, 0, 51)
        sp.row("Smoothing", self.smooth, "points")
        self.norm = sp.check("Scale each trace to 100 %", False)
        self.link = sp.check("Link the time axes", True)
        sp.buttons(flat(sp, "Update chromatograms", "primary", icon="replot", handler=self.on_update))

        sp.section("Mass spectra", T.GROUP["yellow"])
        self.t0, self.t1 = sp.text("", 58), sp.text("", 58)
        sp.row("Average", self.t0, "to", self.t1, "min")
        self.b0, self.b1 = sp.text("", 58), sp.text("", 58)
        sp.row("Background", self.b0, "to", self.b1, "min")
        self.use_bg = sp.check("Subtract background", True)
        self.binw = sp.text("0.05")
        sp.row("Bin width", self.binw, "m/z")
        self.spec_style = sp.choice(["Profile", "Sticks"])
        self.spec_style.SetMinSize(wx.Size(self.FromDIP(120), -1))
        sp.row("Spectrum display", self.spec_style)
        sp.buttons(flat(sp, "Average", handler=self.on_average))
        self.use_ev = sp.choice([])
        self.use_ev.SetMinSize(wx.Size(self.FromDIP(160), -1))
        sp.row("Spectrum to use", self.use_ev)
        sp.buttons(flat(sp, "Export\u2026", icon="export", handler=lambda e: self.on_export_spectrum()))
        self.extra_sections(sp)
        self.build_integration()

        self.build_tools(self.spectrum_tools())
        left = wx.BoxSizer(wx.VERTICAL)
        left.Add(self.tools, 0, wx.EXPAND | wx.BOTTOM, self.FromDIP(10))
        left.Add(self.scroller, 1, wx.EXPAND)
        hs = wx.BoxSizer(wx.HORIZONTAL)
        hs.Add(left, 1, wx.EXPAND | wx.ALL, self.FromDIP(12))
        hs.Add(self.side, 0, wx.EXPAND)
        self.SetSizer(hs)

        self.kind.Bind(wx.EVT_CHOICE, self.on_update)
        self.xic_mode.Bind(wx.EVT_CHOICE, self.on_update)
        self.spec_style.Bind(wx.EVT_CHOICE, lambda e: self._style_changed())
        self.norm.Bind(wx.EVT_CHECKBOX, lambda e: self.plot_chroms(keep_view=True))
        for c in self.xic_fields + [self.tol]:
            c.Bind(wx.EVT_TEXT_ENTER, self.on_update)
        # a list typed by the user is read again; the others keep the windows of
        # their mass chromatograms when the default window changes
        self._xic_dirty = set()
        for i, c in enumerate(self.xic_fields):
            c.Bind(wx.EVT_TEXT, lambda ev, i=i: (self._xic_dirty.add(i), ev.Skip()))
        for c in (self.t0, self.t1, self.b0, self.b1, self.binw):
            c.Bind(wx.EVT_TEXT_ENTER, self.on_average)
        self.smooth.Bind(wx.EVT_SPINCTRL, self.on_update)
        self._build_cols(1)
        self.plot_chroms()
        for col in self.cols:
            self.plot_spec(col)

    # ----------------------------------------------------- extension hooks
    def spectrum_tools(self):
        return [("mode", "xic", "Mass chrom.", "xic", self.TOOL_HINTS["xic"]),
                ("mode", "measure", "Measure", "measure", self.TOOL_HINTS["measure"]),
                ("mode", "label", "Label", "label", self.TOOL_HINTS["label"]),
                ("sep", ""),
                ("action", "deconvolute", "Deconvolute\u2026", "deconvolve", "Deconvolute the spectrum shown")]

    def extra_sections(self, sp):
        """Side panel sections between Mass spectra and Peak integration."""
        import deconv_tab
        self.dec = deconv_tab.DeconvPanel(self, hr=self.PICK_MIN < 0.1)
        self.panels.append(self.dec)

    def reveal(self, win):
        """Scroll the column of tiles so that win is at the top."""
        def go():
            try:
                y = win.GetPosition()[1]
                ux, uy = self.scroller.GetScrollPixelsPerUnit()
                if uy:
                    self.scroller.Scroll(0, max(0, int(y / uy)))
            except RuntimeError:
                pass
        wx.CallAfter(go)

    def keep_in_view(self, wins):
        """Scroll the column of tiles, if one of wins is not wholly in view, so that they are (as much as fits,
        the most of the tiles above them shown too): a mass selected in a deconvolution shows the spectrum with
        its charge states."""
        def go():
            try:
                ux, uy = self.scroller.GetScrollPixelsPerUnit()
                if not uy:
                    return
                spans = [(w.GetPosition()[1], w.GetPosition()[1] + w.GetSize()[1]) for w in wins if w and w.IsShown()]
                if not spans:
                    return
                y0, y1 = min(a for a, b in spans), max(b for a, b in spans)
                vy = self.scroller.GetViewStart()[1] * uy
                vh = self.scroller.GetClientSize()[1]
                if y0 >= vy and y1 <= vy + vh:
                    return
                top = y0 if y1 - y0 > vh else max(0, y1 - vh)
                self.scroller.Scroll(0, max(0, int(top / uy)))
            except RuntimeError:
                pass
        wx.CallAfter(go)

    def set_extra(self, key, row=None, table=None, weight=3.4):
        """Show (or with row and table None, hide) the result rows of the
        deconvolution or calibration, plus the table in the bottom row. row
        is one row or a list of rows, shown in that order below the spectra
        (deconvolution results: the newest first)."""
        old = self.extras.get(key)
        if row is None and table is None:
            self.extras.pop(key, None)
            for w in (_extra_rows(old[0]) + [old[1]] if old else []):
                try:
                    if w is not None:
                        w.Hide()
                except RuntimeError:
                    pass  # already destroyed
        else:
            self.extras[key] = (row, table, weight)
        self._layout_rows()

    # ------------------------------------------------------------ layout
    def _build_cols(self, n):
        """The tiles of n scan events: with one event a chromatogram tile and a spectrum tile, with more one
        tile of chromatograms and one of spectra, each with a plot per event one below the other (PairCard);
        mass chromatograms in between (_build_xic_rows)."""
        keep = {self.bottom} | {r for row, t, w in self.extras.values() for r in _extra_rows(row)}
        for w in list(self.rows.all_panes):
            if w not in keep:
                w.Destroy()
        self.cols, self.views, self.xic_views, self.xic_rows = [], [], [], []
        self.search_views, self.search_rows = [], []
        self._xic_sig = self._search_sig = None
        self.chrom_tile = PairCard(self.rows, "Chromatograms") if n > 1 else None
        self.spec_tile = PairCard(self.rows, "Mass spectra") if n > 1 else None
        for e in range(n):
            card = PlotCard(self.chrom_tile or self.rows, "Chromatogram", "", mode="range")
            view = ChromView(self, card, e)
            view.kind, view.xic, view.file_tag = "main", None, ""
            card.image_name = lambda v=view: self._chrom_image_name(v)
            spc = PlotCard(self.spec_tile or self.rows, "Mass spectrum", "", mode="zoom")
            col = {"e": e, "view": view, "spec_card": spc, "spec": None, "desc": "", "range": None}
            spc.on_menu = lambda x, y, c=col: self.spec_menu(x, y, c)
            spc.on_view = lambda c=col: self.plot_spec(c, keep_view=True)
            spc.readout_fmt = lambda x, y: self.READ_FMT % (x, _g(y))
            spc.image_name = lambda c=col: (self._spec_name(c["e"]) + "_spectrum") if self.data else None
            spc.tool_getter = lambda: self.tool
            spc.tools_supported = {"xic", "measure", "label"}
            spc.on_tool = lambda tool, kind, x, y, c=col: self.spec_tool(c, tool, x, y)
            self.cols.append(col)
        if n > 1:
            self.chrom_tile.set_cards([c["view"].card for c in self.cols])
            self.spec_tile.set_cards([c["spec_card"] for c in self.cols])
        # the tile of the chromatograms and the one of the spectra (reveal, tests)
        self.row_main, self.row_spec = self._chrom_tiles()[0] if n else None, self._spec_tiles()[0] if n else None
        self._layout_rows()

    def _fit_share(self, pane):
        """"Plots fill the window": (share, pixels for the title and axes) of a plot tile (a tile of both
        polarities: a share per plot), None for the tables below."""
        if not self.data or pane is self.bottom:
            return None
        if "dec" in self.extras:
            # of several deconvolution results only the one shown in the table (the others below)
            rows = _extra_rows(self.extras["dec"][0])
            if pane in rows:
                act = next((getattr(getattr(p, "active", None), "row", None) for p in self.panels
                            if getattr(p, "active", None) is not None), None)
                return (1.0, self.FromDIP(100)) if pane is (act if act in rows else rows[0]) else None
        if isinstance(pane, PairCard):
            return float(len(pane.cards) or 1), self.FromDIP(55)
        return 1.0, self.FromDIP(100)

    def _chrom_tiles(self):
        return [self.chrom_tile] if getattr(self, "chrom_tile", None) else [c["view"].card for c in self.cols]

    def _spec_tiles(self):
        return [self.spec_tile] if getattr(self, "spec_tile", None) else [c["spec_card"] for c in self.cols]

    def ev_colour(self, e):
        """Colour of the traces and spectra of scan event e: negative ions red when both polarities are shown."""
        return NEG_COL if len(self.cols) > 1 and self.pol(e) == "-" else TRACE

    def _layout_rows(self):
        defaults = {"main": 4.0, "xic": 3.0, "search": 3.0, "spec": 4.0, "table": 2.2}
        chroms, specs = self._chrom_tiles(), self._spec_tiles()
        srows = getattr(self, "search_rows", [])
        panes = chroms + self.xic_rows + srows + specs
        kinds = ["main"] * len(chroms) + ["xic"] * len(self.xic_rows) + ["search"] * len(srows) + \
            ["spec"] * len(specs)
        tables = [getattr(self, "hrms_details", self.table)]
        for key in ("dec", "cal"):
            if key in self.extras:
                row, table, w = self.extras[key]
                for r in _extra_rows(row):
                    panes.append(r)
                    kinds.append(key)
                defaults[key] = w
                if table is not None:
                    tables.append(table)
        panes.append(self.bottom)
        kinds.append("table")
        tw = [self.bottom.weight_of(t, 1.0) for t in tables]
        self.bottom.set_panes(tables, tw)
        weights = [self.rows.weight_of(p, defaults[k]) for p, k in zip(panes, kinds)]
        self.rows.set_panes(panes, weights)
        self.scroller.refit()
        self.views = [c["view"] for c in self.cols] + self.xic_views + list(getattr(self, "search_views", []))

    def _xic_groups(self, lists):
        """The mass chromatograms shown, in tiles: the same m/z in several events in one tile. [[(e, x)]], in
        the order of the lists (the first of each event, then the second, ...)."""
        groups, where = [], {}
        for k in range(max([len(x) for x in lists] + [0])):
            for e, mlist in enumerate(lists):
                if k >= len(mlist):
                    continue
                x = mlist[k]
                key = round(float(x["mz"]), 4)
                g = where.get(key)
                if g is None or any(e2 == e for e2, x2 in groups[g]):
                    where[key] = g = len(groups)
                    groups.append([])
                groups[g].append((e, x))
        return groups

    def _tile_members(self):
        """What each tile of the column now holds: a set of ("main", e), ("xic", e, m/z), ("spec", e),
        ("extra", j) or ("bottom",) (old_layouts)."""
        evs = [c["e"] for c in self.cols]
        out = [{("main", e) for e in evs}] if getattr(self, "chrom_tile", None) else [{("main", e)} for e in evs]
        lists = [self.xic_shown(e) for e in evs] if self.xic_mode.GetSelection() == 0 else [[] for e in evs]
        for g in self._xic_groups(lists):
            out.append({("xic", evs[i], round(float(x["mz"]), 4)) for i, x in g})
        for s in self._searches_shown():
            out.append({("search", s["id"], e) for e in s["events"]})
        out += [{("spec", e) for e in evs}] if getattr(self, "spec_tile", None) else [{("spec", e)} for e in evs]
        n_extra = len(self.rows.all_panes) - len(out) - 1  # deconvolution and calibration rows
        out += [{("extra", j)} for j in range(max(0, n_extra))] + [{("bottom",)}]
        return out

    def old_layouts(self, old, version=1):
        """Tile layouts of an analysis saved before the plots of the scan events shared tiles
        (ms_project.set_views), as the tiles are now. 4.5 and older (version 1): the chromatograms of the
        events were one row of tiles side by side, the spectra another and each row of mass chromatograms held
        one tile per event; after the column and the bottom row came these rows (in the order they were made)
        and the rows of the deconvolution results. Version 2: a tile per plot, one below the other.
        [] when it cannot be told."""
        n = len(self.cols)
        evs = [c["e"] for c in self.cols]
        lists = [self.xic_shown(e) for e in evs] if self.xic_mode.GetSelection() == 0 else [[] for e in evs]
        new = self._tile_members()
        n_extra = sum(1 for m in new if next(iter(m))[0] == "extra")
        if not old or not n:
            return []
        stack = old[0]
        units = []  # per tile of the old column: [(member, closed)]
        if version >= 2:
            units += [[(("main", e), False)] for e in evs]
            for k in range(max([len(x) for x in lists] + [0])):
                units += [[(("xic", e, round(float(lists[i][k]["mz"]), 4)), False)] for i, e in enumerate(evs)
                          if k < len(lists[i])]
            units += [[(("spec", e), False)] for e in evs]
            rest = list(old[2:])
        else:
            if len(old) < 4 or len(old[2]["weights"]) != n or len(old[3]["weights"]) != n:
                return []
            main, spec, rest = old[2], old[3], list(old[4:])
            k = max([len(x) for x in lists] + [0])
            # the rows of mass chromatograms: k rows of one tile per event, one after the other
            p = next((i for i in range(len(rest) - k + 1) if all(len(r["weights"]) == n for r in rest[i:i + k])),
                     None)
            if p is None:
                return []
            xrows, rest = rest[p:p + k], rest[:p] + rest[p + k:]
            units.append([(("main", e), i in main["closed"]) for i, e in enumerate(evs)])
            for j in range(k):
                units.append([(("xic", e, round(float(lists[i][j]["mz"]), 4)), i in xrows[j]["closed"])
                              for i, e in enumerate(evs) if j < len(lists[i])])
            units.append([(("spec", e), i in spec["closed"]) for i, e in enumerate(evs)])
        units += [[(("extra", j), False)] for j in range(n_extra)] + [[(("bottom",), False)]]
        if len(stack["weights"]) != len(units):
            return []
        frac = {int(i): v for i, v in stack["frac"].items()}
        closed = set()
        for i, u in enumerate(units):
            closed |= {m for m, c in u if c or i in stack["closed"]}
        out = dict(weights=[], closed=[], frac={})
        for j, mem in enumerate(new):
            hit = [i for i, u in enumerate(units) if any(m in mem for m, c in u)]
            out["weights"].append(stack["weights"][hit[0]] if hit else 3.0)
            if hit and mem <= closed:
                out["closed"].append(j)
            fr = [frac[i] for i in hit if i in frac]
            if fr:
                out["frac"][j] = min(0.95, sum(fr))
        return [out, old[1]] + rest

    def xic_shown(self, e):
        """The mass chromatograms of an event that are drawn (not hidden in the files panel tree)."""
        return [x for x in self.xics.get(e, []) if (e, xic_key(x)) not in self.xic_hidden]

    def set_xic_hidden(self, e, x, hidden):
        k = (e, xic_key(x))
        if hidden:
            self.xic_hidden.add(k)
        else:
            self.xic_hidden.discard(k)
        self.on_update()

    def _build_xic_rows(self):
        stacked = self.xic_mode.GetSelection() == 0
        lists = [self.xic_shown(c["e"]) for c in self.cols]
        sig = (stacked, tuple(tuple(xic_key(x) for x in l) for l in lists))
        if sig == self._xic_sig:
            return
        self._xic_sig = sig
        for row in self.xic_rows:
            row.Destroy()
        self.xic_rows, self.xic_views = [], []
        if stacked:
            # a tile per m/z: the same m/z of several events one below the other in it
            for g in self._xic_groups(lists):
                tile = PairCard(self.rows, ("m/z " + self.PIN_FMT) % g[0][1]["mz"]) if len(g) > 1 else None
                cards = []
                for i, x in g:
                    col = self.cols[i]
                    card = PlotCard(tile or self.rows, ("m/z " + self.PIN_FMT + " %s") % (
                        x["mz"], self.ev_short(col["e"])), "", mode="range")
                    v = ChromView(self, card, col["e"])
                    v.kind, v.xic = "xic", x
                    v.file_tag = "_mz%g_%s" % (x["mz"], {"+": "pos", "-": "neg"}.get(self.pol(col["e"]), "ev"))
                    card.image_name = lambda v=v: self._chrom_image_name(v)
                    self.xic_views.append(v)
                    cards.append(card)
                if tile is not None:
                    tile.set_cards(cards)
                    self.xic_rows.append(tile)
                else:
                    self.xic_rows += cards
        self._layout_rows()

    def _chrom_image_name(self, v):
        if not self.data:
            return None
        return self.data.name + "_chromatogram" + (v.file_tag or "_" + {"+": "pos", "-": "neg"}.get(self.pol(v.key), "ev"))

    # -------------------------------------------------------------- data
    def set_data(self, data):
        self.data = data
        self.peaks, self.sel_peak = [], None
        self.avg = self.bg = None
        self.pick_t = None
        for c in (self.t0, self.t1, self.b0, self.b1):
            c.SetValue("")
        self._remember_spectrum_selection()
        n = min(data.n_events, self.MAXCOLS) if data else 1
        if self._search_job is not None:
            self._search_job["cancel"] = True  # (its result belongs to the file read before)
        self.searches = []
        self._build_cols(n)
        self.xics = {e: [] for e in range(n)}
        for i, (t, lab) in enumerate(zip(self.xic_fields, self.xic_labels)):
            t.SetValue("")
            show = data is not None and i < n
            t.Show(show)
            lab.Show(show)
            if show:
                lab.SetLabel("m/z list %s" % self.ev_short(i))
        self.side.Layout()
        self.side.FitInside()
        self.use_ev.Set([self.data.event_label(e) for e in range(n)] if data else [])
        if data:
            self.use_ev.SetSelection(0)
        for p in self.panels:
            p.reset()
        self.on_update()
        for col in self.cols:
            self.plot_spec(col)

    def pol(self, e):
        return self.data.events[e]["polarity"] if self.data else None

    def ev_short(self, e):
        return _pol_short(self.pol(e), e, self.data.events if self.data else None)

    def on_update(self, e=None):
        self.tree_changed()
        if not self.data:
            self.plot_chroms()
            return
        if getattr(self, "_averaging", False):  # the reader is busy in the averaging thread
            self.status("Still averaging the previous range")
            return
        dw, dppm = self.default_window()
        dirty = getattr(self, "_xic_dirty", set())
        for i, fld in enumerate(self.xic_fields):
            if i in self.xics and i in dirty:  # typed in the side panel
                self.xics[i] = parse_xics(fld.GetValue(), dw, dppm)
        dirty.clear()
        if getattr(self, "_xic_default", None) not in (None, (dw, dppm)):
            # a new default window is for new mass chromatograms: the lists show
            # the windows of the others explicitly from now on
            for i in range(len(self.xic_fields)):
                self._set_xic_field(i)
        self._xic_default = (dw, dppm)
        n = self.smooth.GetValue()
        k = self.kind.GetSelection()
        stacked = self.xic_mode.GetSelection() == 0
        self._build_xic_rows()
        self._build_search_rows()
        self._search_traces()
        for col in self.cols:
            ev = col["e"]
            short = self.ev_short(ev)
            traces = []
            if k in (0, 1):
                t, y = self.data.chromatogram(ev, "tic" if k == 0 else "bpc")
                lab = "TIC" if k == 0 else "BPC"
                traces.append({"name": "%s %s" % (lab, short), "label": lab, "t": t, "y": LI.smooth(y, n),
                               "key": "main:%d" % ev})
            xl = self.xic_shown(ev)
            labs = [("m/z " + self.PIN_FMT) % x["mz"] for x in xl]
            for j, x in enumerate(xl):
                t, y = self.data.chromatogram(ev, "xic", x["mz"], xic_tol(x))
                lab = labs[j]
                if labs.count(lab) > 1:
                    # the same m/z with two windows: the names (peaks, integration) tell them apart
                    lab += " " + xic_window_text(x)
                tr = {"name": "%s %s" % (lab, short), "key": "xic:%d:%r" % (ev, xic_key(x)),
                      "label": lab, "t": t, "y": LI.smooth(y, n),
                      "color": PAL[(j + 1) % len(PAL)], "xic": x}
                if stacked:
                    for v in self.xic_views:
                        if v.key == ev and xic_key(v.xic) == xic_key(x):
                            if v.card.pair is not None:  # (with the same m/z of the other polarity: its colour)
                                tr["color"] = self.ev_colour(ev)
                            v.card.side_name, v.card.side_colour = short + " XIC", tr["color"]
                            v.traces = [tr]
                            v.card.set_title(("m/z " + self.PIN_FMT + " %s") % (x["mz"], short),
                                             "window " + xic_window_text(x))
                else:
                    traces.append(tr)
            for i, tr in enumerate(traces):
                tr["color"] = tr.get("color") or (TRACE if i == 0 else PAL[i % len(PAL)])
            if traces:
                traces[0]["color"] = self.ev_colour(ev) if k in (0, 1) else traces[0]["color"]
            col["view"].traces = traces
            col["view"].file_tag = "_" + {"+": "pos", "-": "neg"}.get(self.pol(ev), "ev%d" % (ev + 1))
            e_ = self.data.events[ev]
            sub = "m/z %g to %g" % (e_["mz_low"], e_["mz_high"]) if e_["mz_low"] else ""
            if n >= 3:
                sub += (", " if sub else "") + "smoothed (%d points)" % n
            # (a tile of both polarities: the y title of the plot)
            col["view"].card.side_name = short + " " + ("TIC", "BPC", "XIC")[min(k, 2)]
            col["view"].card.side_colour = self.ev_colour(ev)
            col["view"].card.set_title(_pol_name(self.pol(ev), ev) if self.pol(ev) or not getattr(
                self.data, "labels", None) else self.data.event_label(ev), sub)  # (labels: lcms_sources)
        self.update_trace_choice()
        self.plot_chroms(keep_view=True)

    def empty_text(self, v):
        if self.data and getattr(v, "kind", "") == "search":
            return "Not found"
        if self.data:
            return "No trace chosen"
        return self.OPEN_HINT

    def chrom_extras(self, ax, scale_of, v):
        if not self.data:
            return []
        if self.pick_t is not None and not (self.avg and self.avg[1] != self.avg[0]):
            ln = ax.axvline(self.pick_t, color=C["accent"], lw=0.9, ls=(0, (4, 3)), alpha=0.9, zorder=4)
            ln._no_scale = True
            ln._ui_only = True  # a marker of the view, not part of copied or saved images
        try:
            self._search_marks(ax, scale_of, v)
        except Exception as ex:
            _log("mass search marks: %r" % ex)
        idx = self.data.event_scans(v.key)
        sat = self.data.saturated[idx]
        if np.any(sat) and getattr(v, "kind", "main") == "main":
            from matplotlib.transforms import blended_transform_factory
            tr = blended_transform_factory(ax.transData, ax.transAxes)
            ln, = ax.plot(self.data.rt[idx][sat], np.full(int(sat.sum()), 0.988), ls="none", marker="|",
                          ms=4, mew=0.7, color="#D2456F", alpha=0.55, transform=tr, zorder=5)
            ln._no_scale = True
            return [("#D2456F", "saturated")]
        return []

    def extra_chrom_menu(self, x, y, v):
        if not self.data:
            return []
        items = [(None, None), ("Add a mass chromatogram %s\u2026" % self.ev_short(v.key),
                                lambda: self.ask_xic(v.key))]
        if self.FIND_MASS:
            items.append(("Find mass\u2026", lambda: self.ask_find_mass(v.key)))
        target = None
        if getattr(v, "kind", "") == "xic":
            target = v.xic
        else:  # overlaid: the mass chromatogram trace nearest to the click
            tr = self._trace_near(v, x, y)
            if tr is not None and tr.get("xic"):
                target = tr["xic"]
        if target is not None:
            lab = ("m/z " + self.PIN_FMT) % target["mz"]
            items.append(("Edit this mass chromatogram (%s)\u2026" % lab, lambda t=target: self.edit_xic(v.key, t)))
            items.append(("Remove this mass chromatogram (%s)" % lab, lambda t=target: self.remove_xic(v.key, t)))
        if self.xics.get(v.key):
            items.append(("Remove all mass chromatograms %s" % self.ev_short(v.key), lambda: self.clear_xic(v.key)))
        k, st, n = self.kind.GetSelection(), self.xic_mode.GetSelection(), self.smooth.GetValue()
        items += [(None, None),
                  ("Trace", [(lab, lambda i=i: self._set_choice(self.kind, i), k == i) for i, lab in
                             enumerate(["TIC (total ion current)", "BPC (base peak)", "Mass chromatograms only"])]),
                  ("Mass chromatograms shown", [(lab, lambda i=i: self._set_choice(self.xic_mode, i), st == i)
                                                for i, lab in enumerate(["Stacked tiles", "Overlaid on the TIC"])]),
                  ("Smoothing", [("Off", lambda: self._set_smooth(0), n == 0)] +
                   [("%d points" % m, lambda m=m: self._set_smooth(m), n == m) for m in (3, 5, 7, 9, 15, 25)]),
                  ("Scale each trace to 100 %", lambda: self._flip(self.norm, chrom=True), self.norm.GetValue()),
                  ("Link the time axes", lambda: self._flip(self.link), self.link.GetValue())]
        if self.bg:
            items.append(("Clear the background range", self.clear_bg))
        return items + self.frame.link_menu()

    def _set_choice(self, ctrl, i):
        ctrl.SetSelection(i)
        self.on_update()

    def _set_smooth(self, n):
        self.smooth.SetValue(n)
        self.on_update()

    def _flip(self, cb, chrom=False, spec=False):
        cb.SetValue(not cb.GetValue())
        if chrom:
            self.plot_chroms(keep_view=True)
        if spec:
            for c in self.cols:
                self.plot_spec(c, keep_view=True)

    # ------------------------------------------- method presets (method_presets.py)
    TRACES, XIC_SHOWN, STYLES = ("tic", "bpc", "xic"), ("stacked", "overlaid"), ("profile", "sticks")

    def method_get(self):
        """Settings of this view in a method preset (not the ranges, mass chromatogram lists or peaks)."""
        return dict(trace=self.TRACES[max(0, self.kind.GetSelection())], xic_window=_num(self.tol, None),
                    xic_window_unit=("m/z", "ppm")[max(0, self.tol_unit.GetSelection())],
                    xic_shown=self.XIC_SHOWN[max(0, self.xic_mode.GetSelection())], smoothing=self.smooth.GetValue(),
                    scale_100=self.norm.GetValue(), link_time_axes=self.link.GetValue(),
                    subtract_background=self.use_bg.GetValue(), bin_width=_num(self.binw, None),
                    spectrum_display=self.STYLES[max(0, self.spec_style.GetSelection())], **self.integration_get())

    def method_set(self, v):
        """The values given (the others stay) as if chosen in the panel: computed and drawn again."""
        for key, ctrl, keys in (("trace", self.kind, self.TRACES), ("xic_window_unit", self.tol_unit, ("m/z", "ppm")),
                                ("xic_shown", self.xic_mode, self.XIC_SHOWN),
                                ("spectrum_display", self.spec_style, self.STYLES)):
            if key in v:
                ctrl.SetSelection(keys.index(v[key]))
        for key, ctrl in (("xic_window", self.tol), ("bin_width", self.binw)):
            if key in v:
                ctrl.ChangeValue("%.10g" % v[key])
        for key, cb in (("scale_100", self.norm), ("link_time_axes", self.link), ("subtract_background", self.use_bg)):
            if key in v:
                cb.SetValue(bool(v[key]))
        if "smoothing" in v:
            self.smooth.SetValue(int(v["smoothing"]))
        self.integration_set(v)
        self.on_update()
        if self.data and self.avg and self.avg[1] != self.avg[0]:
            self.compute_average()  # bin width, background
        elif self.data and self.pick_t is not None:
            self.chrom_pick(self.pick_t, None, None)
        else:
            for c in self.cols:
                self.plot_spec(c, keep_view=True)

    def clear_bg(self):
        self.bg = None
        for c in (self.b0, self.b1):
            c.SetValue("")
        self.plot_chroms(keep_view=True)
        if self.avg and self.avg[1] != self.avg[0]:
            self.compute_average()

    def default_window(self):
        w = _num(self.tol, None)
        if w is None or w <= 0:
            w = float(self.TOL_DEFAULT)
        return w, self.tol_unit.GetSelection() == 1

    BINW_RANGE = (0.001, 10.0)  # bin widths accepted (m/z)

    def bin_width(self):
        """Bin width of the averaged and single spectra (m/z) from its field. A
        value that is not a number from 0.001 to 10 (0, negative, text, a
        width that would give millions of bins) is replaced in the field by
        the last valid one, with a note in the status bar."""
        txt = self.binw.GetValue().strip()
        v = _num_s(txt, None) if txt else None
        last = getattr(self, "_binw_ok", 0.05)
        lo, hi = self.BINW_RANGE
        self._binw_note = ""
        if v is None or not (lo <= v <= hi):
            if txt:  # (shown with the result: see _status_binw)
                self._binw_note = "bin width %r replaced by %g m/z: type a number from %g to %g" % (txt, last, lo, hi)
                self.status(self._binw_note[0].upper() + self._binw_note[1:])
            self.binw.ChangeValue("%g" % last)
            return last
        self._binw_ok = v
        return v

    def _status_binw(self, text):
        """Status text, with the note about a bin width that was replaced."""
        note = getattr(self, "_binw_note", "")
        self.status(text + ("; " + note if note else ""))

    def _mz_text(self, v):
        return "%g" % v if self.PICK_MIN >= 0.1 else "%.4f" % v

    def _xic_text(self, x):
        dw, dppm = self.default_window()
        txt = self._mz_text(x["mz"])
        if abs(x["w"] - dw) > 1e-12 or bool(x["ppm"]) != dppm:
            txt += "+-%g%s" % (x["w"], "ppm" if x["ppm"] else "")
        return txt

    def _set_xic_field(self, e):
        if e < len(self.xic_fields):
            # (ChangeValue: no text event, so the list is not read back as typed by the user)
            self.xic_fields[e].ChangeValue(", ".join(self._xic_text(x) for x in self.xics.get(e, [])))

    def ask_xic(self, e, mz=None, entry=None):
        """Dialog for a new mass chromatogram (or to edit one): m/z, window
        and unit, scan event."""
        dw, dppm = self.default_window()
        base = entry or {"mz": mz, "w": dw, "ppm": dppm}
        dlg = XicDialog(self, "Edit mass chromatogram" if entry else "Mass chromatogram",
                        "" if base["mz"] is None else self._mz_text(base["mz"]), base["w"], base["ppm"],
                        [self.data.event_label(c["e"]) for c in self.cols], e, allow_list=entry is None)
        try:
            if dlg.ShowModal() != wx.ID_OK:
                return
            vals = dlg.values()
        finally:
            dlg.Destroy()
        if not vals["mz"]:
            self.status("No m/z value entered")
            return
        events = [c["e"] for c in self.cols] if vals["event"] < 0 else [vals["event"]]
        # the window chosen here becomes the default for quick additions
        self.tol.SetValue("%g" % vals["w"])
        self.tol_unit.SetSelection(1 if vals["ppm"] else 0)
        if entry is not None:
            self.remove_xic(e, entry, update=False)
        for ev in events:
            for m in vals["mz"]:
                self.add_xic(ev, m, update=False, w=vals["w"], ppm=vals["ppm"])
        self.on_update()

    def edit_xic(self, e, entry):
        self.ask_xic(e, entry=entry)

    def add_xic(self, e, m, update=True, w=None, ppm=None):
        dw, dppm = self.default_window()
        x = {"mz": round(float(m), 2 if self.PICK_MIN >= 0.1 else 4), "w": dw if w is None else float(w),
             "ppm": dppm if ppm is None else bool(ppm)}
        vals = self.xics.setdefault(e, [])
        if not any(xic_key(v) == xic_key(x) for v in vals):
            vals.append(x)
        self._set_xic_field(e)
        if update:
            self.on_update()
        ev = self.data.events[e] if self.data and e < len(self.data.events) else {}
        lo, hi = ev.get("mz_low"), ev.get("mz_high")
        if lo and hi and not (lo - xic_tol(x) <= x["mz"] <= hi + xic_tol(x)):
            self.status(("m/z " + self.PIN_FMT + " is outside the scan range of %s (m/z %g to %g): its mass "
                         "chromatogram is empty") % (x["mz"], self.data.event_label(e), lo, hi))

    def remove_xic(self, e, x, update=True):
        self.xics[e] = [v for v in self.xics.get(e, []) if xic_key(v) != xic_key(x)]
        self.xic_hidden.discard((e, xic_key(x)))
        self._set_xic_field(e)
        if update:
            self.on_update()

    def clear_xic(self, e):
        self.xics[e] = []
        self.xic_hidden = set(k for k in self.xic_hidden if k[0] != e)
        self._set_xic_field(e)
        self.on_update()

    # ------------------------------------------------------ Find mass (ms_find_mass)
    def _search_by_id(self, sid):
        return next((s for s in self.searches if s.get("id") == sid), None)

    def _searches_shown(self):
        """The searches with their scan events that have plots here."""
        evs = {c["e"] for c in self.cols}
        out = []
        for s in self.searches:
            e2 = [e for e in s.get("events", []) if e in evs]
            if e2:
                out.append(s if e2 == list(s["events"]) else dict(s, events=e2))
        return out

    def _build_search_rows(self):
        """A tile per search: one plot, or (all scan events) one plot per event in one tile."""
        shown = self._searches_shown()
        sig = (tuple((s["id"], tuple(s["events"])) for s in shown), len(self.cols))
        if sig == self._search_sig:
            return
        self._search_sig = sig
        for row in self.search_rows:
            row.Destroy()
        self.search_rows, self.search_views = [], []
        for s in shown:
            tile = PairCard(self.rows, "Mass search") if len(s["events"]) > 1 else None
            cards = []
            for e in s["events"]:
                card = PlotCard(tile or self.rows, "Mass search", "", mode="range")
                v = ChromView(self, card, e)
                v.kind, v.xic, v.search_id = "search", None, s["id"]
                v.file_tag = "_mass_search_%s" % {"+": "pos", "-": "neg"}.get(self.pol(e), "ev%d" % (e + 1))
                card.image_name = lambda v=v: self._chrom_image_name(v)
                self.search_views.append(v)
                cards.append(card)
            if tile is not None:
                tile.set_cards(cards)
                self.search_rows.append(tile)
            else:
                self.search_rows += cards
        self._layout_rows()

    def _search_traces(self):
        """The traces, titles and notes of the search plots."""
        n = self.smooth.GetValue()
        names = set()
        for v in self.search_views:
            s = self._search_by_id(v.search_id)
            e = v.key
            short = self.ev_short(e)
            res = (s or {}).get("res", {}).get(e) or []
            rows = search_rows(s) if s else []
            one_w = len({w for M, w, nm in rows}) <= 1
            traces, missing, few = [], [], []
            for j, ((M, w, nm), r) in enumerate(zip(rows, res)):
                mt = _mass_text(M)
                lab = (nm or mt) + ("" if one_w else " \u00b1 %g" % w)
                if not r.get("regions"):
                    (few if r.get("few") else missing).append(lab)
                    continue
                name = "%s Da search %s" % (mt, short)
                if name in names:
                    name = "%s Da search %d %s" % (mt, s["id"], short)
                names.add(name)
                traces.append({"name": name, "label": lab, "t": r["t"], "y": LI.smooth(r["trace"], n),
                               "color": SEARCH_COLS[j % len(SEARCH_COLS)], "search": r, "mass_name": nm,
                               "key": "search:%d:%d:%s:%g" % (s["id"], e, mt, w)})
            v.traces = traces
            parts = ["\u00b1 %g Da" % rows[0][1]] if rows and one_w else []
            if s and self.pol(e) != "-" and s.get("carrier", "H+") != "H+":
                parts.append(s["carrier"])
            if missing:
                parts.append("%s not found" % ", ".join(missing))
            if few:
                parts.append("%s not found (too few charge states)" % ", ".join(few))
            card = v.card
            card.side_name, card.side_colour = short + " Search", self.ev_colour(e)
            card.set_title("Mass search" + (" " + short if len(self.cols) > 1 and card.pair is None else ""),
                           ", ".join(parts))

    def _search_note(self, s, e):
        """The harmonic note of each mass: another mass that fits the same m/z."""
        rel = {"half": "half this mass", "double": "double this mass",
               "one and a half times": "one and a half times this mass", "two thirds": "two thirds of this mass"}
        res = (s or {}).get("res", {}).get(e) or []
        rows = search_rows(s) if s else []
        lines = []
        for (M, w, nm), r in zip(rows, res):
            h = r.get("harmonic") if r.get("regions") else None
            if not h:
                continue
            txt = "These m/z also fit %s Da (%s), it may be the same species" % (
                _mass_text(round(h[0], 1)), rel.get(h[1], h[1]))
            lines.append(txt if len(rows) == 1 else "%s: %s" % (nm or _mass_text(M), txt[0].lower() + txt[1:]))
        return "\n".join(lines)

    def _search_marks(self, ax, scale_of, v):
        """The found regions shaded (search plots and the TIC of the event), the label of the best region and
        the harmonic note."""
        import ms_find_mass as F
        kind = getattr(v, "kind", "main")
        if kind == "search":
            v._labels = []  # (placed once the plot has its axis ranges: plot_view)
            for tr in v.traces:
                r, col = tr["search"], tr["color"]
                for a, b, g in F.shaded(r):
                    ax.axvspan(a, b, color=col, alpha=0.12, lw=0, zorder=0)
                g = F.best_region(r)
                if g is None:
                    continue
                y = tr["y"] * scale_of.get(tr["name"], 1.0)
                sel = (tr["t"] >= g["t0"]) & (tr["t"] <= g["t1"])
                if not np.any(sel):
                    continue
                i = int(np.argmax(np.where(sel, y, -np.inf)))
                st = g.get("states") or []
                ms = g.get("masses") or [g["mass"]]  # (a wide window can hold several species)
                txt = (tr["mass_name"] + ", " if tr.get("mass_name") else "") + (
                    " and ".join("%.1f" % m for m in ms) + " Da") + (
                    ", z %d to %d" % (min(st), max(st)) if st else "")
                ann = ax.annotate(txt, (tr["t"][i], y[i]), xytext=(8, -2), textcoords="offset points", fontsize=8,
                                  color=col, va="top", zorder=6)
                v._labels.append(ann)
            note = self._search_note(self._search_by_id(v.search_id), v.key)
            if note:
                ax.text(0.01, 0.96, note, transform=ax.transAxes, ha="left",
                        va="top", fontsize=8, color=NOTE_COL, zorder=7, linespacing=1.3,
                        bbox=dict(boxstyle="square,pad=0.2", fc="white", ec="none", alpha=0.8))
        elif kind == "main":
            for s in self._searches_shown():
                if v.key not in s["events"]:
                    continue
                for j, r in enumerate(s["res"].get(v.key) or []):
                    for a, b, g in F.shaded(r):
                        sp = ax.axvspan(a, b, color=SEARCH_COLS[j % len(SEARCH_COLS)], alpha=0.10, lw=0, zorder=0)
                        sp._range_mark = True

    def plot_view(self, v, keep_view=False):
        TabBase.plot_view(self, v, keep_view)
        if getattr(v, "kind", "") == "search" and getattr(v, "_labels", None):
            placed = []
            for ann in v._labels:  # the labels of the masses clear of each other
                placed = _clear_of(v.card.ax, ann, placed)
            v._labels = []
            v.card.canvas.draw_idle()

    def chrom_menu(self, x, y, v):
        if getattr(v, "kind", "") != "search":
            return TabBase.chrom_menu(self, x, y, v)
        s = self._search_by_id(v.search_id)
        items = []
        if s is not None:
            sid = s["id"]
            items.append(("Edit search…", lambda: self.ask_find_mass(v.key, search_id=sid)))
            rows = search_rows(s)
            if len(rows) > 1:
                items.append(("Remove a mass", [((nm + ", " if nm else "") + _mass_text(M) + " Da",
                                                 lambda k=k: self.remove_search_mass(sid, k))
                                                for k, (M, w, nm) in enumerate(rows)]))
        items += [("Find mass…", lambda: self.ask_find_mass(v.key)), (None, None),
                  ("Average a time range…", self.ask_average), (None, None),
                  ("Full view", v.card.reset_view), ("Copy chromatogram data", lambda: self.copy_chrom(v)),
                  ("Export chromatogram…", lambda: self.export_chrom(v))]
        return items + self.frame.link_menu() + self.view_menu()

    def _search_pick(self, view, x, y):
        """A click inside a found region: the averaged spectrum of its main part (trace >= 20 % of its top) with
        the charge states of the mass marked. False: not inside one."""
        import ms_find_mass as F
        best = None
        for tr in view.traces:
            for a, b, g in F.shaded(tr["search"]):
                if a <= x <= b:
                    yd = float(np.interp(x, tr["t"], tr["y"] * tr.get("_scale", 1.0)))
                    d = abs(yd - y) if y is not None else 0.0
                    if best is None or d < best[0]:
                        best = (d, tr, g)
        if best is None:
            return False
        if self.busy_averaging():
            return True
        tr, g = best[1], best[2]
        r = tr["search"]
        a, b = F.core(r, g, 0.2)
        self.pick_t = None
        self.avg = (a, b)
        self.fill_range_fields()
        self.plot_chroms(keep_view=True)
        self.compute_average()
        if a == b or any(c.get("range") == (a, b) for c in self.cols):
            for col in self.cols:
                if col["e"] == view.key and col["spec"] is not None:
                    # marked at the mass found (inside a wide window it is not the typed mass), the marks as
                    # narrow as the profile of the search
                    col["smark"] = {"M": float((g.get("masses") or [g.get("mass", r["M"])])[0]), "c": float(r["c"]),
                                    "w": min(float(r["w"]), 5.0),
                                    "zs": [int(z) for z in r["zs"]], "on": [int(z) for z in g.get("states") or []],
                                    "colour": tr["color"], "sign": -1 if r.get("polarity") == "-" else 1,
                                    "name": tr.get("mass_name") or ""}
                    self.plot_spec(col)
                    # the m/z range of the charge states (the other spectra of the tile follow)
                    card = col["spec_card"]
                    ev = self.data.events[col["e"]]
                    lo = max(ev.get("mz_low") or 0.0, float(r.get("mzmin", 0.0)) - 20.0)
                    hi = ev.get("mz_high") or (card.full[1] if card.full else None)
                    if card.full and hi and hi > lo:
                        set_mz_range(card, lo, hi)
            if a != b:
                self.frame.link_time("ms", rng=self.avg)
            self.status("Spectrum of %s, %.2f to %.2f min" % (tr.get("mass_name") or "%s Da" % _mass_text(r["M"]),
                                                              a, b))
        return True

    def ask_find_mass(self, e=None, search_id=None):
        """Find mass: the dialog (a table of masses), then the search (off the interface thread)."""
        if not self.data or not self.FIND_MASS:
            self.status("Open a data file first")
            return
        if self._search_job is not None:
            self.status("A mass search is running")
            return
        s = self._search_by_id(search_id)
        last = settings().get("find_mass")
        last = last if isinstance(last, dict) else {}
        n = len(self.cols)
        if s is not None:
            v = {"rows": [(_mass_text(M), "%g" % w, nm) for M, w, nm in search_rows(s)], "carrier": s["carrier"],
                 "mzmin": s["mzmin"]}
            ev = -1 if len(s["events"]) > 1 else s["events"][0]
        else:
            v = {k: last[k] for k in ("carrier", "mzmin") if k in last}
            rows = last.get("rows")
            v["rows"] = [tuple(str(x) for x in r[:3]) for r in rows if isinstance(r, (list, tuple)) and
                         len(r) >= 3] if isinstance(rows, list) else []
            ev = -1 if (last.get("event") == -1 and n > 1) else (e if e is not None and 0 <= e < n else 0)
        dlg = FindMassDialog(self, "Edit mass search" if s is not None else "Find mass", v,
                             [self.data.event_label(c["e"]) for c in self.cols], [self.pol(c["e"]) for c in self.cols],
                             ev)
        try:
            if dlg.ShowModal() != wx.ID_OK:
                return
            vals = dlg.values()
        finally:
            dlg.Destroy()
        T._save({"find_mass": {"rows": [list(r) for r in vals["table"]], "carrier": vals["carrier"],
                               "mzmin": vals["mzmin"], "event": -1 if vals["event"] < 0 else 0}})
        if not vals["rows"]:
            self.status("No mass entered")
            return
        self.find_mass(vals, replace=search_id if s is not None else None)

    def find_mass(self, vals, replace=None):
        """Search the masses vals["rows"] ([{"M", "w", "name"}]; carrier, mzmin shared) in the scan event
        vals["event"] (-1: every event shown) in a worker thread; the tile is made (or the search replace
        replaced) when it is done. Progress in the status bar and a small window with Cancel."""
        import ms_find_mass as F
        if self._search_job is not None or not self.data:
            return False
        events = [c["e"] for c in self.cols] if vals["event"] < 0 else [vals["event"]]
        rows = [{"M": float(r["M"]), "w": float(r.get("w", 5.0)), "name": str(r.get("name") or "")}
                for r in vals["rows"]]
        data, carrier, mzmin = self.data, vals["carrier"], float(vals["mzmin"])
        shorts = {e: self.ev_short(e) for e in events}
        n = max(1, len(events) * len(rows))
        job = {"vals": dict(vals, rows=rows), "events": events, "replace": replace, "data": data,
               "cancel": False, "frac": 0.0, "now": "", "done": False, "out": None, "err": None,
               "cancelled": False, "t0": time.time(), "win": None}

        def work():
            try:
                out, k = {}, 0
                for e in events:
                    lst = []
                    for r in rows:
                        what = r["name"] or "%s Da" % _mass_text(r["M"])
                        job["now"] = "%s %s" % (what, shorts[e]) if len(self.cols) > 1 else what
                        lst.append(F.search(data, e, r["M"], r["w"], carrier, mzmin, cancel=lambda: job["cancel"],
                                            progress=lambda x, k=k: job.__setitem__("frac", (k + x) / n)))
                        k += 1
                    out[e] = lst
                job["out"] = out
            except F.Cancelled:
                job["cancelled"] = True
            except Exception as ex:
                import traceback
                traceback.print_exc()
                job["err"] = ex
            finally:
                job["done"] = True
        self._search_job = job
        threading.Thread(target=work, name="find mass", daemon=True).start()
        if getattr(self, "_search_timer", None) is None:
            self._search_timer = wx.Timer(self)
            self.Bind(wx.EVT_TIMER, self._search_tick, self._search_timer)
        self._search_timer.Start(120)
        self.status("Find mass …")
        return True

    def _search_tick(self, e=None):
        job = self._search_job
        if job is None:
            self._search_timer.Stop()
            return
        if not job["done"]:
            txt = "Find mass: %s, %d %%" % (job["now"], int(100 * job["frac"]))
            self.status(txt)
            if job["win"] is None and time.time() - job["t0"] > 0.8 and not job["cancel"]:
                import deconv_tab
                job["win"] = deconv_tab.ProgressWindow(wx.GetTopLevelParent(self), "Find mass",
                                                       lambda: job.__setitem__("cancel", True),
                                                       cancel_tip="Stop the search")
            if job["win"] is not None:
                job["win"].update(int(1000 * job["frac"]), 1000, txt)
            return
        self._search_timer.Stop()
        self._search_job = None
        if job["win"] is not None:
            job["win"].finish()
        self._search_finished(job)

    def _search_finished(self, job):
        if job["cancelled"]:
            self.status("Mass search cancelled")
            return
        if job["err"] is not None:
            _report_error(self, "The mass search could not be finished", job["err"])
            return
        if self.data is not job["data"] or job["out"] is None:
            return  # (another file was read meanwhile)
        vals = job["vals"]
        rows = vals["rows"]
        old = self._search_by_id(job["replace"])
        sid = old["id"] if old is not None else max([s["id"] for s in self.searches] + [0]) + 1
        s = {"id": sid, "masses": [r["M"] for r in rows], "ws": [r["w"] for r in rows],
             "names": [r["name"] for r in rows], "carrier": vals["carrier"], "mzmin": float(vals["mzmin"]),
             "events": list(job["events"]), "res": job["out"]}
        if old is not None:
            self.searches = [s if o is old else o for o in self.searches]
        else:
            self.searches = self.searches + [s]
        self._searches_changed()
        found = [k for k in range(len(rows)) if any(lst[k]["regions"] for lst in s["res"].values())]
        import undo
        undo.commit(self, ("edit the mass search %s" if old is not None else "find mass %s") % search_text(s))
        if len(rows) == 1:
            self.status("%s %s" % (rows[0]["name"] or "%s Da" % _mass_text(rows[0]["M"]),
                                   "found" if found else "not found"))
        else:
            self.status("%d of %d masses found" % (len(found), len(rows)))

    def _searches_changed(self):
        self._build_search_rows()
        self._search_traces()
        self.update_trace_choice()
        self.plot_chroms(keep_view=True)
        self.tree_changed()

    def remove_search_mass(self, sid, k):
        """Row k of the search taken out (the others keep their results)."""
        s = self._search_by_id(sid)
        if s is None:
            return
        rows = search_rows(s)
        if not 0 <= k < len(rows):
            return
        keep = [j for j in range(len(rows)) if j != k]
        if not keep:
            self.remove_search(sid)
            return
        new = dict(s, masses=[rows[j][0] for j in keep], ws=[rows[j][1] for j in keep],
                   names=[rows[j][2] for j in keep],
                   res={e: [lst[j] for j in keep if j < len(lst)] for e, lst in s["res"].items()})
        new.pop("w", None)
        self.searches = [new if o is s else o for o in self.searches]
        self._searches_changed()
        self.status("%s removed from the mass search" % (rows[k][2] or "%s Da" % _mass_text(rows[k][0])))

    def remove_search(self, sid):
        self.searches = [s for s in self.searches if s.get("id") != sid]
        self._searches_changed()

    def _close_search_tile(self, pane):
        """Close this tile on a search tile: the search is removed (Undo brings it back)."""
        for v in self.search_views:
            if v.card is pane or (v.card.pair is not None and v.card.pair is pane):
                self.remove_search(v.search_id)
                self.status("Mass search closed")
                return True
        return False

    # ---------------------------------------------------------- spectra
    def chrom_pick(self, x, y, view):
        if not self.data:
            return
        if getattr(self, "_averaging", False):  # the reader is busy in the averaging thread
            self.status("Still averaging the previous range")
            return
        if getattr(view, "kind", "") == "search" and self._search_pick(view, x, y):
            return  # (inside a found region: the spectrum of its main part)
        self.avg = None
        self.pick_t = x
        self.t0.SetValue("%.3f" % x)
        self.t1.SetValue("")
        binw = self.bin_width()
        for col in self.cols:
            try:
                spec, k = self.data.scan_spectrum(col["e"], x, binw=binw)
            except Exception as ex:
                _report_error(self, "The spectrum at %.3f min could not be read" % x, ex)
                return
            # time and flags of that scan from the scans of its event: k numbers
            # the scans of the file, while rt and saturated hold them event by
            # event in the C++ reader (EngineFile), interleaved in lcms_data
            idx = self.data.event_scans(col["e"])
            j = idx[int(np.argmin(np.abs(self.data.rt[idx] - x)))] if len(idx) else k
            sat = getattr(self.data, "saturated", None)
            col["spec"] = spec
            self.spectrum_replaced(col)
            col["range"] = None
            col["measure"] = None
            col["m1"] = None
            col["desc"] = "scan %d at %.3f min%s" % (k + 1, self.data.rt[j],
                                                    ", saturated signal"
                                                    if sat is not None and len(sat) > j and sat[j] else "")
            self.plot_spec(col, keep_view=True)
        self._remember_spectrum_selection()
        self.plot_chroms(keep_view=True)
        self.spectra_changed()
        self.frame.link_time("ms", t=x)
        if getattr(self, "_binw_note", ""):
            self._status_binw("Spectra at %.3f min" % x)

    def spectrum_replaced(self, col):
        # Formula measurements belong to the previous spectrum, not its view.
        col.pop("overlay", None)
        col.pop("smark", None)  # (the charge states of a found mass: that spectrum only)

    def _remember_spectrum_selection(self):
        self._display_selection = (self.avg, self.bg, self.pick_t, self.use_bg.GetValue())

    def _restore_spectrum_selection(self):
        self.avg, self.bg, self.pick_t, use_bg = getattr(
            self, "_display_selection", (None, None, None, self.use_bg.GetValue()))
        self.use_bg.SetValue(use_bg)
        a = self.avg
        values = ((_fmt_t(a[0]), _fmt_t(a[1])) if a and a[0] != a[1] else
                  (_fmt_t(self.pick_t), "") if self.pick_t is not None else ("", ""))
        for ctrl, value in zip((self.t0, self.t1), values):
            ctrl.ChangeValue(value)
        for ctrl, value in zip((self.b0, self.b1), self.bg or (None, None)):
            ctrl.ChangeValue(_fmt_t(value) if value is not None else "")
        self.plot_chroms(keep_view=True)

    def spectra_changed(self):
        dec = getattr(self, "dec", None)
        if dec is not None:
            dec.refresh_inputs(quiet=True)

    def range_changed(self, shift):
        if not shift:
            self.pick_t = None
        if not shift or any(c["range"] for c in self.cols):
            self.compute_average()
        if not shift and self.avg and self.avg[1] != self.avg[0]:
            self.frame.link_time("ms", rng=self.avg)

    def average_fields(self):
        self.on_average()

    def on_average(self, e=None):
        if self.busy_averaging():
            return
        note = self.read_range_fields()
        if not self.avg:
            self.status("Enter a time range or drag across a chromatogram")
            return
        if self.avg[1] == self.avg[0]:
            self.chrom_pick(self.avg[0], None, None)
        else:
            self.pick_t = None
            self.plot_chroms(keep_view=True)
            self.compute_average()
        if note:
            self.status(note)

    def compute_average(self):
        if not self.data or not self.avg:
            return
        if getattr(self, "_averaging", False):
            # a range set from elsewhere (PDA view, Compare view, calibrant,
            # background, method) while averaging: averaged when this run ends
            self._average_next = True
            self.status("Still averaging the previous range; the new range follows")
            return
        t0, t1 = self.avg
        if t1 == t0:
            self.chrom_pick(t0, None, None)
            return
        bg_set = self.bg
        bg = self.bg if (self.use_bg.GetValue() and self.bg) else None
        binw = self.bin_width()
        keep = any(c["spec"] is not None for c in self.cols)
        try:
            out = self._average_all(t0, t1, bg, binw)
        except Exception as ex:
            self.__dict__.pop("_average_next", None)
            self._restore_spectrum_selection()
            _report_error(self, "The spectra of %.2f to %.2f min could not be averaged" % (t0, t1), ex)
            return
        if self.__dict__.pop("_average_next", False) and out is not None:
            # the range set meanwhile follows; until then the view keeps the
            # range of the spectra just averaged
            wx.CallAfter(self._average_later, self.avg, self.bg)
            self.avg, self.bg = (t0, t1), bg_set
        if out is None:
            self._restore_spectrum_selection()
            return
        if not any(n for d, n in out):
            rt = self.data.rt
            # the spectra shown stay, and so does the range they belong to
            self._restore_spectrum_selection()
            self.status("No scans between %.2f and %.2f min (the run goes from %.2f to %.2f min)"
                        % (t0, t1, float(rt.min()) if len(rt) else 0.0, float(rt.max()) if len(rt) else 0.0))
            return
        for col, (d, n) in zip(self.cols, out):
            col["spec"] = d
            self.spectrum_replaced(col)
            col["range"] = (t0, t1)
            col["measure"] = None
            col["m1"] = None
            note = ""
            if bg:
                # a background range without scans of this event subtracts nothing: say so
                idx = self.data.event_scans(col["e"])
                rt = self.data.rt[idx]
                has_bg = bool(np.any((rt >= min(bg)) & (rt <= max(bg))))
                note = (", background %.2f to %.2f min subtracted" % bg if has_bg else
                        ", no scans in the background range %.2f to %.2f min (nothing subtracted)" % bg)
            col["desc"] = "%.2f to %.2f min, %d scans%s" % (t0, t1, n, note)
            self.plot_spec(col, keep_view=keep)
        self._remember_spectrum_selection()
        self._status_binw("Averaged %.2f to %.2f min" % (t0, t1))
        self.spectra_changed()

    def _average_later(self, avg, bg):
        """A range set while the previous one was still being averaged."""
        if not self or not self.data or not avg:
            return
        self.avg, self.bg = avg, bg
        self.fill_range_fields()
        self.plot_chroms(keep_view=True)
        self.range_changed(False)  # averages it and moves the linked PDA view there again

    def _average_all(self, t0, t1, bg, binw):
        """The averaged spectrum of every scan event, or None if cancelled.
        Readers with a progress callback (HRMS: hundreds of scans of 600000
        points) run in a thread while a progress window with Cancel keeps
        the interface alive; the others run directly."""
        import inspect
        try:
            with_progress = "progress" in inspect.signature(self.data.average).parameters
        except (TypeError, ValueError):
            with_progress = False
        if not with_progress:
            wx.BeginBusyCursor()
            try:
                return [self.data.average(col["e"], t0, t1, bg=bg, binw=binw) for col in self.cols]
            finally:
                wx.EndBusyCursor()
        import threading
        state = {"cancel": False, "i": 0, "n": 0, "col": 0, "out": [], "err": None}

        def progress(i, n):
            state["i"], state["n"] = i, n
            return not state["cancel"]

        def work():
            try:
                for k, col in enumerate(self.cols):
                    state["col"] = k
                    state["out"].append(self.data.average(col["e"], t0, t1, bg=bg, binw=binw, progress=progress))
            except Exception as ex:
                state["err"] = ex
        self._averaging = True
        th = threading.Thread(target=work, name="average", daemon=True)
        th.start()
        win = None
        t_start = time.time()
        wx.BeginBusyCursor()
        try:
            while th.is_alive():
                th.join(0.05)
                if win is None and time.time() - t_start > 0.8:
                    import deconv_tab
                    win = deconv_tab.ProgressWindow(wx.GetTopLevelParent(self), "Averaging the spectrum",
                                                    lambda: state.__setitem__("cancel", True))
                if win is not None:
                    try:
                        n = state["n"]
                        txt = "%.2f to %.2f min: scan %d of %d" % (t0, t1, state["i"] + 1, n) if n else "reading"
                        if len(self.cols) > 1:
                            txt += " (%s)" % self.ev_short(self.cols[state["col"]]["e"])
                        win.update(state["i"], n, txt)
                    except RuntimeError:
                        pass
                wx.GetApp().Yield(True)
        finally:
            self._averaging = False
            wx.EndBusyCursor()
            if win is not None:
                win.finish()
        if state["cancel"]:
            self.status("Averaging cancelled")
            return None
        err = state["err"]
        if err is not None:
            if type(err).__name__ == "Cancelled":
                self.status("Averaging cancelled")
                return None
            raise err
        return state["out"]

    def plot_spec(self, col, keep_view=False):
        self.tree_changed()
        card = col["spec_card"]
        ax = card.ax
        view = ax.get_xlim() if (keep_view and card.full) else None
        card.reset()
        style_axes(ax, "m/z", "Intensity")
        e = col["e"]
        title = "Mass spectrum %s" % self.ev_short(e) if self.data else "Mass spectrum"
        card.side_name = (self.ev_short(e) + " MS") if self.data else "MS"  # (y title in a tile of both polarities)
        card.side_colour = self.ev_colour(e)
        spec = col["spec"]
        if spec is None or len(spec) == 0:
            msg = ("Drag across a chromatogram to average a time range" if self.data
                   else "The mass spectrum appears here")
            ax.text(0.5, 0.5, msg, transform=ax.transAxes, ha="center", va="center", fontsize=9, color=MUTED)
            card.full = None
            card.set_legend([])
            card.set_title(title, "")
            card.draw()
            return
        lo, hi = shared_mz_range(self.data, [c for c in self.cols if c["spec_card"].pair is card.pair] if
                                 card.pair is not None else [col])
        apex = self.APEX_LABELS and bool(getattr(self.data, "has_profile", False))
        sm = col.get("smark")  # the charge states of a found mass (Find mass: a found region clicked)
        # room above the peaks for the charge tags of a found mass or of the mass selected in a deconvolution
        # (inside the frame in a tile of both polarities, which has no room above it)
        dmarks = card.pair is not None and any(getattr(p, "marks_for", lambda c: False)(col) for p in self.panels)
        card.label_room_pt = 26.0 if (sm or dmarks) else 16.0
        card.top_tags = dmarks
        auto = draw_spectrum(card, spec, (lo, hi), self.spec_style.GetSelection() == 1, view, 0 if sm else 8,
                             self.MZ_FMT, self.apex_fn if apex else None, color=self.ev_colour(e))
        draw_spec_marks(ax, col, auto, self.PIN_FMT, lambda m: self._peak_height(col, m))
        if sm:
            try:
                draw_charge_marks(ax, spec, sm)
            except Exception as ex:
                _log("charge state marks: %r" % ex)
        if sm:  # (in the plot: a note of the title row would cover the labels of the tallest peaks)
            ax.plot([0.012], [0.935], transform=ax.transAxes, ls="none", marker="v", ms=5, mew=0,
                    color=sm["colour"], zorder=9)[0]._no_scale = True
            ax.text(0.022, 0.96, "charge states of %s" % (sm.get("name") or "%s Da" % _mass_text(sm["M"])),
                    transform=ax.transAxes, ha="left", va="top", fontsize=8, color=sm["colour"], zorder=9)
        sub = col["desc"]
        try:
            sp_ = col.get("spec")
            if sp_ is not None and len(sp_):
                ib = int(np.argmax(sp_[:, 1]))
                bx = self.apex_fn(sp_[:, 0], sp_[:, 1], ib) if apex else sp_[ib, 0]
                sub += "   \u00b7   base peak m/z " + (self.PIN_FMT % bx) + " (" + _g(sp_[ib, 1]) + ")"
        except Exception:
            pass
        card.set_title(title, sub)
        for p in self.panels:
            try:
                p.overlay(col, ax)
            except Exception as ex:
                _log("overlay failed: %s" % ex)
        card.draw()

    def _nearest_peak(self, col, x):
        return spec_nearest_peak(col["spec"], col["spec_card"].ax.get_xlim(), x, self.PICK_MIN)

    def spectrum_actions(self, col):
        """First entries of the spectrum menu (dialogs)."""
        return [("Deconvolute this spectrum\u2026", lambda: self.on_deconvolve(col["e"]))]

    def spec_menu(self, x, y, col):
        items = []
        e = col["e"]
        has = col["spec"] is not None and self.data is not None
        if has:
            items += self.spectrum_actions(col) + [(None, None)]
        pk = self._nearest_peak(col, x)
        if pk is not None and self.data:
            items.append((("Mass chromatogram of m/z " + self.PIN_FMT + " %s") % (pk, self.ev_short(e)),
                          lambda m=pk: self.add_xic(e, m)))
            items.append((("Mass chromatogram of m/z " + self.PIN_FMT + " with options\u2026") % pk,
                          lambda m=pk: self.ask_xic(e, mz=m)))
            if len(self.cols) > 1:
                items.append((("Mass chromatogram of m/z " + self.PIN_FMT + " in every scan event") % pk,
                              lambda m=pk: self._xic_all(m)))
            items.append(spec_pin_item(col, pk, self.PIN_FMT, lambda m: self.spec_tool(col, "label", m, None)))
            items.append((None, None))
        items += spec_clear_items(col, lambda: self.plot_spec(col, keep_view=True))
        items.append(spec_display_item(self.spec_style.GetSelection() == 1, self._set_style))
        if self.data:
            items.append(("Average a time range\u2026", self.ask_average))
        if has:
            items += [("Export spectrum\u2026", lambda: self.on_export_spectrum(e)),
                      ("Copy spectrum data", lambda: self.copy_spec(e))]
            if self.FULL_WINDOW:
                items.append(("Open this spectrum in the Deconvolute window", lambda: self.open_full(e)))
        if getattr(self, "dec", None) is not None and self.dec.results:
            items += [(None, None), ("Close the deconvolution results", self.dec.close)]
        items += [(None, None)] + mz_range_items(self, col["spec_card"])
        return items + self.view_menu()

    def _set_style(self, i):
        self.spec_style.SetSelection(i)
        self._style_changed()

    def _style_changed(self):
        """Profile or sticks chosen: the spectra drawn again, also the mass spectra of the PDA view."""
        for c in self.cols:
            self.plot_spec(c, keep_view=True)
        pda = getattr(self.frame, "pda", None)
        if pda is not None and getattr(pda, "mscols", None):
            for c in pda.mscols:
                pda.plot_ms(c, keep_view=True)

    def all_cards(self):
        return [v.card for v in self.views] + [c["spec_card"] for c in self.cols]

    def spec_tool(self, col, tool, x, y):
        if col["spec"] is None:
            return
        pk = self._nearest_peak(col, x)
        if pk is None:
            self.status("Click closer to a peak")
            return
        e = col["e"]
        if tool == "xic":
            self.add_xic(e, pk)
            self.status(("Mass chromatogram of m/z " + self.PIN_FMT + " %s added") % (pk, self.ev_short(e)))
        elif spec_click(col, tool, pk, self.PIN_FMT, lambda m1, m2: self._measure_text(e, m1, m2), self.status):
            self.plot_spec(col, keep_view=True)

    def _measure_text(self, e, m1, m2):
        return measure_text(m1, m2, self.data.adduct_sign(e) if self.data else 1)

    def _peak_height(self, col, m):
        return spec_peak_height(col["spec"], m, self.PICK_MIN)

    def _xic_all(self, m):
        for col in self.cols:
            self.add_xic(col["e"], m, update=False)
        self.on_update()

    def copy_spec(self, e):
        copy_spectrum(self, self.cols[e]["spec"])

    def _spec_name(self, e):
        col = self.cols[e]
        pol = {"+": "pos", "-": "neg"}.get(self.pol(e), "ev%d" % (e + 1))
        if pol in ("pos", "neg") and _pol_shared(self.data.events, e):
            pol += "_ev%d" % (e + 1)  # two events of one polarity: one file each
        if col["range"]:
            return "%s_%s_%.2f-%.2fmin" % (self.data.name, pol, col["range"][0], col["range"][1])
        return "%s_%s_scan" % (self.data.name, pol)

    def _chosen(self, e):
        if e is None:
            e = max(0, self.use_ev.GetSelection())
        return e if 0 <= e < len(self.cols) else 0

    def on_export_spectrum(self, e=None):
        if not self.data:
            self.status("Open a data file first")
            return
        e = self._chosen(e)
        col = self.cols[e]
        if col["spec"] is None:
            self.status("No spectrum to export")
            return
        path = self.frame.ask_save("Export mass spectrum %s" % self.ev_short(e), self._spec_name(e) + ".txt",
                                   SPEC_WILDCARD)
        if not path:
            return
        try:
            write_spectrum_file(path, col["spec"], self.data.name, self.pol(e), col["range"])
        except Exception as ex:
            _report_error(self, _write_failed(path), ex)
            return
        self.status("Saved " + path)

    def on_deconvolve(self, e=None):
        """Deconvolute the spectrum of scan event e in the Deconvolution tab."""
        if not self.data:
            self.status("Open a data file first")
            return
        e = self._chosen(e)
        if self.cols[e]["spec"] is None:
            self.status("Average a time range or click a time first")
            return
        self.dec.deconvolute(e)

    def open_full(self, e=None):
        """Open the spectrum in the full Deconvolute window (separate)."""
        if not self.data:
            return
        e = self._chosen(e)
        col = self.cols[e]
        if col["spec"] is None:
            return
        open_spectrum_in_deconvolute(self, self.frame, col["spec"], self._spec_name(e), self.data.adduct_sign(e),
                                     col.get("range"))

    def deconv_inputs(self):
        out = []
        if not self.data:
            return out
        for col in self.cols:
            if col["spec"] is None or not len(col["spec"]):
                continue
            e = col["e"]
            if col["range"]:
                what = "%.2f to %.2f min" % col["range"]
            else:
                what = col["desc"].split(",")[0]
            out.append({"label": "%s: %s" % (self.data.event_label(e), what), "spec": col["spec"],
                        "sign": self.data.adduct_sign(e), "name": self._spec_name(e), "e": e,
                        "time": ("%s averaged" % what) if col["range"] else what})
        return out

    def peak_selected(self, p):
        if not self.data or self.busy_averaging():
            return
        self.avg = (p["t0"], p["t1"])
        self.fill_range_fields()
        self.plot_chroms(keep_view=True)
        self.compute_average()


def open_in_unidec(path):
    """Open the file in the Deconvolute window (own process)."""
    here = os.path.dirname(os.path.abspath(__file__))
    launcher = os.path.join(here, "launch_unidec.py")
    exe = sys.executable
    if os.name == "nt" and exe.lower().endswith("python.exe"):
        w = exe[:-10] + "pythonw.exe"
        if os.path.isfile(w):
            exe = w
    try:
        subprocess.Popen([exe, launcher, path], cwd=os.path.dirname(here), close_fds=True)
    except Exception as e:
        wx.MessageBox("The Deconvolute window could not be started:\n%s" % e, TITLE, wx.ICON_ERROR)


# ==========================================================================
# PDA tab: heat map with sliders, chromatogram and UV spectrum follow them
# ==========================================================================
class MapCard(PlotCard):
    """Wavelength-time heat map. The handles on the top edge (time) and the
    right edge (wavelength) are sliders; the crosshair can be dragged too."""

    def __init__(self, parent, on_cross, on_step):
        PlotCard.__init__(self, parent, "Wavelength map", "", mode="map")
        self.on_cross, self.on_step = on_cross, on_step
        self.margins = [0.76, 0.48, 0.80, 0.30]
        self.t = self.wl = None
        self.full_y = None
        self.arts = {}
        self._drag = None
        self.canvas.mpl_connect("key_press_event", self._key)

    def setup(self, t, wl):
        from matplotlib.transforms import blended_transform_factory as blend, offset_copy
        ax, fig = self.ax, self.fig
        top = blend(ax.transData, ax.transAxes)
        right = blend(ax.transAxes, ax.transData)
        acc = C["accent"]
        tr_top = offset_copy(ax.transAxes, fig=fig, y=10, units="points")
        tr_right = offset_copy(ax.transAxes, fig=fig, x=10, units="points")
        for xs, ys, tr in (([0, 1], [1, 1], tr_top), ([1, 1], [0, 1], tr_right)):
            ln, = ax.plot(xs, ys, transform=tr, color="#DCE3EE", lw=5, solid_capstyle="round", clip_on=False)
            ln._no_scale = True
            ln._ui_only = True
        a = {}
        a["v"] = ax.axvline(t, color="white", lw=1.0, alpha=0.95)
        a["h"] = ax.axhline(wl, color="white", lw=1.0, alpha=0.95)
        a["v"]._ui_only = a["h"]._ui_only = True
        a["th"], = ax.plot([t], [1], transform=offset_copy(top, fig=fig, y=10, units="points"), marker="o", ms=9,
                           color=acc, mec="white", mew=1.5, clip_on=False, ls="none")
        a["wh"], = ax.plot([1], [wl], transform=offset_copy(right, fig=fig, x=10, units="points"), marker="o",
                           ms=9, color=acc, mec="white", mew=1.5, clip_on=False, ls="none")
        a["tl"] = ax.text(t, 1, "", transform=offset_copy(top, fig=fig, x=9, y=10, units="points"), ha="left",
                          va="center", fontsize=8, color=C["accent_text"], fontweight="bold", clip_on=False,
                          bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.9))
        a["wl"] = ax.text(1, wl, "", transform=offset_copy(right, fig=fig, x=18, units="points"), ha="left",
                          va="center", fontsize=8, color=C["accent_text"], fontweight="bold", clip_on=False,
                          bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.9))
        for art in a.values():
            art._no_scale = True
            art._ui_only = True
            self.add_live(art)
        self.arts = a
        self.set_cross(t, wl, redraw=False)

    def set_cross(self, t, wl, redraw=True):
        self.t, self.wl = t, wl
        a = self.arts
        if not a:
            return
        a["v"].set_xdata([t, t])
        a["h"].set_ydata([wl, wl])
        a["th"].set_xdata([t])
        a["wh"].set_ydata([wl])
        a["tl"].set_position((t, 1))
        a["tl"].set_text("%.2f min" % t)
        a["wl"].set_position((1, wl))
        a["wl"].set_text("%.0f nm" % wl)
        x0, x1 = sorted(self.ax.get_xlim())
        y0, y1 = sorted(self.ax.get_ylim())
        for k in ("v", "th", "tl"):
            a[k].set_visible(x0 <= t <= x1)
        for k in ("h", "wh", "wl"):
            a[k].set_visible(y0 <= wl <= y1)
        if redraw:
            self.refresh_cross()

    def refresh_cross(self):
        self.blit_live()

    # ------------------------------------------------------------- mouse
    def _mode_at(self, e):
        if not self.arts or self.t is None:
            return None
        bb = self.ax.bbox
        x, y = e.x, e.y
        if bb.x0 - 10 <= x <= bb.x1 + 10 and y > bb.y1 - 2:
            return "t"
        if x > bb.x1 - 2 and bb.y0 - 10 <= y <= bb.y1 + 10:
            return "w"
        if bb.x0 <= x <= bb.x1 and bb.y0 <= y <= bb.y1:
            px, py = self.ax.transData.transform((self.t, self.wl))
            near = self.FromDIP(6)
            if abs(x - px) <= near:
                return "t"
            if abs(y - py) <= near:
                return "w"
            return "tw"
        return None

    def _apply(self, e, final):
        tx, wy = self.ax.transData.inverted().transform((e.x, e.y))
        t = tx if "t" in self._drag else self.t
        wl = wy if "w" in self._drag else self.wl
        self.on_cross(t, wl, final)

    def _down(self, e):
        try:
            self.canvas.SetFocus()
        except Exception:
            pass
        if e.button == 1 and not e.dblclick and self._tool() in ("zoom", "pan"):
            return PlotCard._down(self, e)
        if e.button == 1:
            if e.dblclick:
                self.reset_view()
                return
            self._drag = self._mode_at(e)
            if self._drag:
                self._apply(e, False)
        elif e.button == 3 and self.on_menu and e.inaxes is self.ax:
            items = self.on_menu(e.xdata, e.ydata) or []
            if items:
                wx.CallAfter(self._popup, items)

    def _move(self, e):
        if self._press is not None or (self._drag is None and self._tool() in ("zoom", "pan")):
            return PlotCard._move(self, e)
        if self._drag:
            self._apply(e, False)
            return
        mode = self._mode_at(e)
        cur = {"t": wx.CURSOR_SIZEWE, "w": wx.CURSOR_SIZENS, "tw": wx.CURSOR_CROSS}.get(mode, wx.CURSOR_ARROW)
        self.canvas.SetCursor(wx.Cursor(cur))
        if self.readout_fmt is not None:
            ok = e.inaxes is self.ax and e.xdata is not None and self.full is not None
            self._set_readout(self.readout_fmt(e.xdata, e.ydata) if ok else "")

    def _up(self, e):
        if self._press is not None:
            return PlotCard._up(self, e)
        if self._drag:
            self._apply(e, True)
            self._drag = None

    def _scroll(self, e):
        if e.inaxes is not self.ax or e.xdata is None:
            return
        f = 0.8 if e.button == "up" else 1.25
        shift = bool(e.key and "shift" in e.key) or wx.GetKeyState(wx.WXK_SHIFT)
        if shift and self.full_y:
            y0, y1 = self.ax.get_ylim()
            a, b = e.ydata - (e.ydata - y0) * f, e.ydata + (y1 - e.ydata) * f
            a, b = max(a, self.full_y[0]), min(b, self.full_y[1])
            if b > a:
                self.ax.set_ylim(a, b)
        else:
            x0, x1 = self.ax.get_xlim()
            a, b = e.xdata - (e.xdata - x0) * f, e.xdata + (x1 - e.xdata) * f
            if self.full:
                a, b = max(a, self.full[0]), min(b, self.full[1])
            if b > a:
                self.ax.set_xlim(a, b)
        self.set_cross(self.t, self.wl, redraw=False)
        self.canvas.draw_idle()

    def reset_view(self):
        self.ylock = None
        if self.full:
            self.ax.set_xlim(*self.full)
        if self.full_y:
            self.ax.set_ylim(*self.full_y)
        self.set_cross(self.t, self.wl, redraw=False)
        self.canvas.draw_idle()

    def _key(self, e):
        k = (e.key or "").replace("shift+", "")
        big = 10 if (e.key or "").startswith("shift+") else 1
        step = {"left": (-big, 0), "right": (big, 0), "up": (0, big), "down": (0, -big)}.get(k)
        if step and self.on_step:
            self.on_step(*step)


class PDATab(TabBase):
    ylabel = "Absorbance (mAU)"
    units = "mAU"

    def __init__(self, parent, frame):
        TabBase.__init__(self, parent, frame)
        self.pda = None
        self.ms = None
        self.t_sel = None
        self.wl_sel = None
        self.uv = None
        self.uv_desc = ""
        self._live = False
        self._last = (None, None)
        self._uv_line = self._uv_marker = self._chrom_marker = None
        self._uv_labels = []
        self.scroller = RowScroller(self)
        self.rows = StackBox(self.scroller)
        self.rows.shape = tile_shape
        # by default the map, the chromatogram and UV spectrum and the mass spectra fill the window, the table
        # below (its own setting: the MS view keeps its choice)
        self.rows.fit_share = self._fit_share
        self.rows.fit_label = "Plots fill the window"
        self.rows.fit_floor = 60  # (plots lower than this: a very small window scrolls instead)
        self.rows.mode_key, self.rows.default_mode = "pda_tile_mode", "two"
        self.row2 = SplitBox(self.rows, wx.HORIZONTAL, min_size=200)
        self.map_card = MapCard(self.rows, self.on_cross, self.step_cross)
        self.map_card.on_menu = self.map_menu
        self.map_card.readout_fmt = self.map_readout
        self.chrom = PlotCard(self.row2, "Chromatogram", "", mode="range")
        self.view = ChromView(self, self.chrom, 0)
        self.view.file_tag = "_PDA"
        self.views = [self.view]
        self.spec_card = PlotCard(self.row2, "UV spectrum", "", mode="zoom")
        self.spec_card.on_pick = lambda x, y: self.on_cross(self.t_sel, x, True)
        self.spec_card.on_menu = self.spec_menu
        self.spec_card.on_view = lambda: self.plot_uv(keep_view=True)
        self.spec_card.readout_fmt = lambda x, y: "%.1f nm   %.1f mAU" % (x, y)
        self.spec_card.margins = [0.78, 0.52, 0.14, 0.14]
        self.table = PeakTable(self.rows, self.select_peak, self.delete_peak)
        # the mass spectra of the time chosen, one tile per scan event below the chromatogram (files with MS data)
        self.mscols = []
        self._ms_key = None
        self.map_card.image_name = lambda: self.frame.base_name() + "_PDA_map"
        self.chrom.image_name = lambda: self.frame.base_name() + "_PDA_chromatogram_%.0fnm" % (self.wl_sel or 0)
        self.spec_card.image_name = lambda: self.frame.base_name() + "_UV_%.2fmin" % (self.t_sel or 0)
        self.uv_pins = []
        self.map_card.tool_getter = lambda: self.tool
        self.map_card.on_view = lambda: self.map_card.set_cross(self.map_card.t, self.map_card.wl, redraw=False)
        self.spec_card.tool_getter = lambda: self.tool
        self.spec_card.tools_supported = {"label"}
        self.spec_card.on_tool = lambda tool, kind, x, y: self.uv_tool(x)

        sp = self.side
        sp.section("Map and sliders", T.GROUP["blue"])
        self.f_time = sp.text("", 70)
        sp.row("Time", self.f_time, "min")
        self.f_wl = sp.text("", 70)
        sp.row("Wavelength", self.f_wl, "nm")
        self.bw = sp.text("4", 70)
        sp.row("Bandwidth", self.bw, "nm")
        self.cscale = sp.choice(["Compressed", "Linear"])
        self.cscale.SetMinSize(wx.Size(self.FromDIP(120), -1))
        sp.row("Colour scale", self.cscale)

        sp.section("Chromatogram", T.GROUP["yellow"])
        self.wls = sp.text("", 120, tooltip="Separated by commas")
        sp.row("More wavelengths", self.wls, "nm")
        self.maxplot = sp.check("Max plot", False)
        self.smooth = sp.spin(0, 0, 51)
        sp.row("Smoothing", self.smooth, "points")
        self.norm = sp.check("Scale each trace to 100 %", False)
        self.overlay = sp.check("Overlay the MS trace", False)
        self.ms_event = sp.choice([])
        self.ms_event.SetMinSize(wx.Size(self.FromDIP(150), -1))
        sp.row("MS trace", self.ms_event)
        self.delay = sp.text("0.00")
        sp.row("MS detector delay", self.delay, "min")

        sp.section("UV spectrum", T.GROUP["yellow"])
        self.t0, self.t1 = sp.text("", 58), sp.text("", 58)
        sp.row("Average", self.t0, "to", self.t1, "min")
        self.b0, self.b1 = sp.text("", 58), sp.text("", 58)
        sp.row("Background", self.b0, "to", self.b1, "min")
        self.use_bg = sp.check("Subtract background", True)
        sp.buttons(flat(sp, "Average", handler=self.on_show_uv),
                   flat(sp, "Export\u2026", icon="export", handler=self.on_export_uv))
        self.build_integration()

        hs = wx.BoxSizer(wx.HORIZONTAL)
        self.row2.set_panes([self.chrom, self.spec_card], [3, 2])
        self._layout_tiles()
        self.scroller.set_box(self.rows)
        self.build_tools([("mode", "label", "Label", "label", self.TOOL_HINTS["label"])])
        left = wx.BoxSizer(wx.VERTICAL)
        left.Add(self.tools, 0, wx.EXPAND | wx.BOTTOM, self.FromDIP(10))
        left.Add(self.scroller, 1, wx.EXPAND)
        hs.Add(left, 1, wx.EXPAND | wx.ALL, self.FromDIP(12))
        hs.Add(self.side, 0, wx.EXPAND)
        self.SetSizer(hs)

        self.f_time.Bind(wx.EVT_TEXT_ENTER, lambda e: self.on_cross(_num(self.f_time, self.t_sel), self.wl_sel, True))
        self.f_wl.Bind(wx.EVT_TEXT_ENTER, lambda e: self.on_cross(self.t_sel, _num(self.f_wl, self.wl_sel), True))
        for c in (self.bw, self.wls, self.delay):
            c.Bind(wx.EVT_TEXT_ENTER, self.on_update)
        for c in (self.t0, self.t1, self.b0, self.b1):
            c.Bind(wx.EVT_TEXT_ENTER, self.on_show_uv)
        for c in (self.maxplot, self.overlay):
            c.Bind(wx.EVT_CHECKBOX, self.on_update)
        self.cscale.Bind(wx.EVT_CHOICE, lambda e: self.plot_map())
        self.norm.Bind(wx.EVT_CHECKBOX, lambda e: self.plot_chroms(keep_view=True))
        self.ms_event.Bind(wx.EVT_CHOICE, self.on_update)
        self.smooth.Bind(wx.EVT_SPINCTRL, self.on_update)
        self.plot_chroms()
        self.plot_uv()
        self.plot_map()

    def empty_text(self, v):
        if self.frame.path and not self.pda:
            return "This file has no PDA data"
        return "Open a data file with PDA data"

    def set_data(self, pda, ms):
        self.pda, self.ms = pda, ms
        self.peaks, self.sel_peak = [], None
        self.avg = self.bg = None
        self.uv = None
        for c in (self.t0, self.t1, self.b0, self.b1):
            c.SetValue("")
        self.ms_event.Set([ms.event_label(e) for e in range(ms.n_events)] if ms else [])
        if ms:
            self.ms_event.SetSelection(0)
        self.overlay.Enable(ms is not None)
        self._build_ms_tiles()
        self._last = (None, None)
        if pda is not None:
            w = pda.wavelengths
            want = getattr(self, "method_wl", None)  # (a method preset applied before the data were read)
            want = want if want is not None and w[0] <= want <= w[-1] else 280
            self.wl_sel = float(w[np.argmin(np.abs(w - want))]) if w[0] <= want <= w[-1] else float(np.median(w))
            t, y = pda.chromatogram(self.wl_sel, _num(self.bw, 4.0) or 4.0)
            self.t_sel = float(t[int(np.argmax(y))])
        else:
            self.t_sel = self.wl_sel = None
        self.plot_map()
        self._sync_fields()
        self.on_update()
        self.compute_uv()
        if pda is None:
            self.plot_uv()  # (compute_uv draws only with PDA data)

    # ---------------------------------------------------- slider logic
    def _snap(self, t, wl):
        p = self.pda
        t = float(p.times[int(np.argmin(np.abs(p.times - (t if t is not None else p.times[0]))))])
        wl = float(p.wavelengths[int(np.argmin(np.abs(p.wavelengths - (wl if wl is not None else 280.0))))])
        return t, wl

    def _sync_fields(self):
        if self.t_sel is not None:
            self.f_time.SetValue("%.3f" % self.t_sel)
            self.f_wl.SetValue("%.1f" % self.wl_sel)

    def on_cross(self, t, wl, final):
        if self.pda is None:
            return
        t, wl = self._snap(t, wl)
        changed = (t, wl) != (self.t_sel, self.wl_sel)
        t_moved = t != self.t_sel
        self.t_sel, self.wl_sel = t, wl
        self._sync_fields()
        self.map_card.set_cross(t, wl)
        if final:
            self._apply_live(force=True)
            if t_moved or getattr(self, "_link_pending", False):
                self._link_pending = False
                self.frame.link_time("pda", t=t)
        elif changed:
            if t_moved:
                self._link_pending = True  # the MS view follows on release
            if not self._live:  # the UV spectrum and chromatogram follow while dragging
                self._live = True
                wx.CallLater(45, self._apply_live)

    def step_cross(self, di, dj):
        if self.pda is None or self.t_sel is None:
            return
        p = self.pda
        i = int(np.argmin(np.abs(p.times - self.t_sel))) + di
        j = int(np.argmin(np.abs(p.wavelengths - self.wl_sel))) + dj
        i = min(max(i, 0), len(p.times) - 1)
        j = min(max(j, 0), len(p.wavelengths) - 1)
        self.on_cross(p.times[i], p.wavelengths[j], True)

    def _apply_live(self, force=False):
        self._live = False
        if self.pda is None:
            return
        _t0 = time.perf_counter()
        t_changed = self._last[0] != self.t_sel
        wl_changed = self._last[1] != self.wl_sel
        if not (t_changed or wl_changed or force):
            return
        self._last = (self.t_sel, self.wl_sel)
        had_range = bool(self.avg and self.avg[1] != self.avg[0])
        if t_changed:
            # a new time replaces an averaged time range; a new wavelength
            # alone keeps it (the UV average stays, only the chromatogram changes)
            self._clear_range()
        self._ms_hold = not force  # the mass spectra follow when the slider is let go
        try:
            if not force and not had_range and self._drag_update(t_changed, wl_changed):
                pass
            elif wl_changed or force:
                self.on_update(keep=True)
                self.compute_uv()
            else:
                self.plot_chroms(keep_view=True)
                self.compute_uv()
        finally:
            self._ms_hold = False
        if os.environ.get("UNILCMS_DEBUG"):
            _log("live update %.0f ms (time %s, wavelength %s)" % ((time.perf_counter() - _t0) * 1000,
                                                                  t_changed, wl_changed))

    def _drag_update(self, t_changed, wl_changed):
        """Fast update while a slider is dragged: only the live artists are
        redrawn (y ranges stay; the full replot follows on release)."""
        tr = self.view.traces[0] if self.view.traces else None
        if tr is None or tr.get("_line") is None or self.uv is None or getattr(self, "_uv_line", None) is None:
            return False
        if wl_changed:
            bw = max(_num(self.bw, 4.0) or 0.0, 0.0)
            t, y = self.pda.chromatogram(self.wl_sel, bw)
            y = LI.smooth(y, self.smooth.GetValue())
            tr["y"] = y
            tr["_line"].set_ydata(y * tr.get("_scale", 1.0))
            if self._uv_marker is not None:
                self._uv_marker.set_xdata([self.wl_sel, self.wl_sel])
            tr["label"] = "%.0f nm (slider)" % self.wl_sel
            self.chrom.legend_items = [(x["color"], x["label"]) for x in self.view.traces] + self.chrom.legend_items[len(self.view.traces):]
            self.chrom.set_title("Chromatogram", "bandwidth %g nm" % bw)
        if t_changed:
            if self._chrom_marker is not None:
                self._chrom_marker.set_xdata([self.t_sel, self.t_sel])
            bg = self.bg if (self.use_bg.GetValue() and self.bg) else None
            if bg and not self._has_spectra(*bg):
                bg = None  # (as compute_uv: nothing is subtracted)
            wl, a = self.pda.spectrum(self.t_sel, None, bg=bg)
            self.uv = (wl, a)
            self.uv_desc = "at %.3f min" % self.t_sel + (", background %.2f to %.2f min subtracted" % bg if bg else "")
            self._uv_line.set_ydata(a)
            self.spec_card.remove_live(self._uv_labels)
            self._uv_labels = [self.spec_card.add_live(x) for x in self._label_uv(self.spec_card.ax, wl, a)]
            self.spec_card.set_title("UV spectrum", self.uv_desc)
        self.chrom.blit_live()
        self.spec_card.blit_live()
        return True

    # ------------------------------------------------------ chromatogram
    def on_update(self, e=None, keep=True):
        self.tree_changed()
        if self.pda is None:
            self.view.traces = []
            self.plot_chroms()
            return
        bw = max(_num(self.bw, 4.0) or 0.0, 0.0)
        n = self.smooth.GetValue()
        traces = []
        main = self.wl_sel
        t, y = self.pda.chromatogram(main, bw)
        traces.append({"name": "%.0f nm" % main, "label": "%.0f nm (slider)" % main, "t": t, "y": LI.smooth(y, n),
                       "key": "slider"})
        outside = []
        for wl in _numlist(self.wls.GetValue()):
            if not (self.pda.wavelengths[0] - 2 <= wl <= self.pda.wavelengths[-1] + 2):
                outside.append("%g" % wl)
            elif abs(wl - main) > 0.6:
                t, y = self.pda.chromatogram(wl, bw)
                traces.append({"name": "%g nm" % wl, "label": "%g nm" % wl, "t": t, "y": LI.smooth(y, n),
                               "key": "wl:%g" % wl})
        if outside:
            self.status("No chromatogram at %s nm: the PDA data cover %.0f to %.0f nm"
                        % (", ".join(outside), self.pda.wavelengths[0], self.pda.wavelengths[-1]))
        if self.maxplot.GetValue():
            t, y = self.pda.max_plot()
            traces.append({"name": "Max plot", "label": "Max plot", "t": t, "y": LI.smooth(y, n), "key": "max"})
        for i, tr in enumerate(traces):
            tr["color"] = PAL[i % len(PAL)] if i else TRACE
        self.view.traces = traces
        self.update_trace_choice()
        sub = "bandwidth %g nm" % bw
        if n >= 3:
            sub += ", smoothed (%d points)" % n
        self.chrom.set_title("Chromatogram", sub)
        self.plot_chroms(keep_view=keep and bool(self.chrom.full))
        self.compute_ms()  # (a new MS detector delay)

    def chrom_extras(self, ax, scale_of, v):
        extra = []
        self._chrom_marker = None
        if self.t_sel is not None and not (self.avg and self.avg[1] != self.avg[0]):
            ln = ax.axvline(self.t_sel, color=C["accent"], lw=0.9, ls=(0, (4, 3)), alpha=0.9, zorder=4)
            ln._ui_only = True
            ln._no_scale = True
            self._chrom_marker = v.card.add_live(ln)
        if v.traces:
            v.card.add_live(v.traces[0].get("_line"))
        if self.overlay.GetValue() and self.ms is not None and v.traces:
            ev = max(0, self.ms_event.GetSelection())
            t, y = self.ms.chromatogram(ev, "tic")
            d = _num(self.delay, 0.0) or 0.0
            top = max(float(np.nanmax(tr["y"] * scale_of.get(tr["name"], 1.0))) for tr in v.traces)
            ym = float(np.nanmax(y)) if len(y) else 0.0
            if ym > 0:
                lab = "TIC %s%s" % (_pol_short(self.ms.events[ev]["polarity"], ev, self.ms.events),
                                    ", shifted by %g min" % -d if d else "")
                ln, = ax.plot(t - d, y / ym * top, color="#8A94A3", lw=LW_TRACE, alpha=0.8, zorder=2, label=lab)
                ln._no_scale = True
                extra.append(("#8A94A3", lab))
        return extra

    def chrom_pick(self, x, y, v):
        self._clear_range()  # a click shows the spectrum at one time, also at the time already set
        self.on_cross(x, self.wl_sel, True)

    def _clear_range(self):
        self.avg = None
        for c in (self.t0, self.t1):
            c.SetValue("")

    # ------------------------------------------------------ UV spectrum
    def range_changed(self, shift):
        self.compute_uv()
        if not shift and self.avg and self.avg[1] != self.avg[0]:
            self.frame.link_time("pda", rng=self.avg)

    def on_show_uv(self, e=None):
        note = self.read_range_fields()
        if not self.avg:
            self.status("Enter a time range, or drag across the chromatogram")
            return
        if self.avg[1] == self.avg[0]:
            # one time typed: the spectrum at that time, as a click on the chromatogram
            t = self.avg[0]
            self._clear_range()
            self.on_cross(t, self.wl_sel, True)
        else:
            self.plot_chroms(keep_view=True)
            self.compute_uv()
        if note:
            self.status(note)

    def _has_spectra(self, a, b):
        """The time range a..b meets the PDA run (a range narrower than the
        spectrum interval inside it uses the nearest spectrum)."""
        tt = self.pda.times
        if not len(tt):
            return False
        dt = float(getattr(self.pda, "interval_s", 0.0) or 0.0) / 60.0
        return max(a, b) >= tt[0] - dt and min(a, b) <= tt[-1] + dt

    def compute_uv(self):
        if self.pda is None:
            return
        bg = self.bg if (self.use_bg.GetValue() and self.bg) else None
        note = ""
        if bg and not self._has_spectra(*bg):
            # (lcms_pda would take the spectrum nearest to the range instead)
            note = ", no spectra in the background range %.2f to %.2f min (nothing subtracted)" % bg
            bg = None
        if self.avg and self.avg[1] != self.avg[0]:
            t0, t1 = self.avg
            if not self._has_spectra(t0, t1):
                tt = self.pda.times
                self.status("No PDA spectra between %.2f and %.2f min (the PDA data go from %.2f to %.2f min)"
                            % (t0, t1, tt[0], tt[-1]))
                self._clear_range()
                self.plot_chroms(keep_view=True)
                t0 = t1 = None
        if self.avg and self.avg[1] != self.avg[0]:
            wl, a = self.pda.spectrum(t0, t1, bg=bg)
            desc = "average %.2f to %.2f min" % (t0, t1)
        else:
            wl, a = self.pda.spectrum(self.t_sel, None, bg=bg)
            desc = "at %.3f min" % self.t_sel
        if bg:
            desc += ", background %.2f to %.2f min subtracted" % bg
        desc += note
        self.uv = (wl, a)
        self.uv_desc = desc
        self.plot_uv(keep_view=True)
        self.compute_ms()

    def plot_uv(self, keep_view=False):
        self.tree_changed()
        card = self.spec_card
        ax = card.ax
        view = ax.get_xlim() if (keep_view and card.full) else None
        card.reset()
        style_axes(ax, "Wavelength (nm)", "Absorbance (mAU)")
        self._uv_line = self._uv_marker = None
        self._uv_labels = []
        if self.uv is None:
            ax.text(0.5, 0.5, "Move the time slider on the map" if self.pda is not None else "No PDA data",
                    transform=ax.transAxes, ha="center", va="center", fontsize=9, color=MUTED)
            card.full = None
            card.set_title("UV spectrum", "")
            card.draw()
            return
        wl, a = self.uv
        self._uv_line = card.add_live(ax.plot(wl, a, color=TRACE, lw=LW_TRACE)[0])
        ax.axhline(0, color="#C8CFD9", lw=0.6, zorder=0)
        self._uv_marker = None
        if self.wl_sel is not None:
            ln = ax.axvline(self.wl_sel, color=C["accent"], lw=0.9, ls=(0, (4, 3)), alpha=0.9)
            ln._ui_only = True
            ln._no_scale = True
            self._uv_marker = card.add_live(ln)
        card.full = (float(wl[0]), float(wl[-1])) if wl[-1] > wl[0] else (wl[0] - 1.0, wl[0] + 1.0)
        if len(wl) == 1:  # one wavelength: a point
            ax.plot(wl, a, "o", color=TRACE, ms=3)
        ax.set_xlim(*(view if view else card.full))
        card.autoscale_y(pad=0.12)
        self._uv_labels = [card.add_live(t) for t in self._label_uv(ax, wl, a)]
        card.set_title("UV spectrum", self.uv_desc)
        card.draw()

    def _label_uv(self, ax, wl, a):
        x0, x1 = ax.get_xlim()
        arts = label_maxima(ax, wl, a, n=5, fmt="%.0f", min_sep=(x1 - x0) / 12.0, xlim=(x0, x1))
        for w in self.uv_pins:
            if not (x0 <= w <= x1):
                continue
            for q in list(arts):
                if abs(q.xy[0] - w) < (x1 - x0) / 40.0:
                    q.remove()
                    arts.remove(q)
            j = int(np.argmin(np.abs(wl - w)))
            arts.append(ax.annotate("%.0f nm\n%.0f mAU" % (wl[j], a[j]), (wl[j], a[j]), xytext=(0, 5),
                                    textcoords="offset points", ha="center", va="bottom", fontsize=7.5,
                                    color=C["accent_text"], fontweight="bold", zorder=5,
                                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec=C["accent"], lw=0.6,
                                              alpha=0.95)))
        return arts

    def uv_tool(self, x):
        if self.uv is None:
            return
        wl, a = self.uv
        j = int(np.argmin(np.abs(wl - x)))
        lo, hi = max(0, j - 3), min(len(wl) - 1, j + 3)
        j = lo + int(np.argmax(a[lo:hi + 1]))  # snap to the local band maximum
        w = float(wl[j])
        near = [p for p in self.uv_pins if abs(p - w) < 1.5]
        if near:
            self.uv_pins.remove(near[0])
        else:
            self.uv_pins.append(w)
        self.plot_uv(keep_view=True)

    def all_cards(self):
        return [self.map_card, self.chrom, self.spec_card] + [c["spec_card"] for c in self.mscols]

    # ------------------------------------------------------ mass spectra
    MS_WEIGHT = 4.0

    def _ms_tiles(self):
        tile = getattr(self, "ms_tile", None)
        return [tile] if tile is not None else [c["spec_card"] for c in self.mscols]

    def _layout_tiles(self):
        ms = self._ms_tiles()
        panes = [self.map_card, self.row2] + ms + [self.table]
        defaults = [7, 5] + [self.MS_WEIGHT] * len(ms) + [2.2]
        self.rows.set_panes(panes, [self.rows.weight_of(p, w) for p, w in zip(panes, defaults)])
        self.scroller.refit()

    # the tiles that fill the window ("Plots fill the window"): (plots, DIP for the title, sliders and axes);
    # each tile gets its DIP and a part of the rest for each plot in it, so the plots are about as high
    FIT_SHARES = {"map": (1.0, 125), "row2": (1.0, 110), "ms_pair": (2.0, 55), "ms": (1.0, 100)}

    def _fit_share(self, pane):
        if self.pda is None:
            return None  # (no PDA data: the view is not shown; the tiles as before)
        if pane is self.map_card:
            key = "map"
        elif pane is self.row2:
            key = "row2"
        elif pane in self._ms_tiles():
            key = "ms_pair" if isinstance(pane, PairCard) else "ms"
        else:
            return None
        share, dip = self.FIT_SHARES[key]
        return share, self.FromDIP(dip)

    def old_heights(self, old):
        """Tile layouts of an analysis saved before the plots filled the window (ms_project.set_views): the
        heights the tiles were dragged to then are left out (they would push the mass spectra out of the
        window again); closed tiles stay closed."""
        return [dict(st, frac={}) if i == 0 else st for i, st in enumerate(old)]

    def old_layouts(self, old, version=1):
        """Tile layouts of an analysis saved before (ms_project.set_views): 4.5 and older (version 1) had no
        mass spectra in this view (map, chromatogram and UV spectrum, table); version 2 a tile per mass
        spectrum."""
        m = len(self._ms_tiles())
        if not old:
            return old
        st = old[0]
        k = len(self.mscols) if version >= 2 else 0  # the mass spectrum tiles saved
        if len(st["weights"]) != 3 + k:
            return [] if version >= 2 else old
        # each tile now: the tiles saved that held it
        src = [[0], [1]] + ([list(range(2, 2 + k))] if m == 1 and k > 1 else [[2 + j] if k else [] for j in range(m)])
        src.append([2 + k])
        frac = {int(i): v for i, v in st["frac"].items()}
        new = dict(weights=[st["weights"][s_[0]] if s_ else self.MS_WEIGHT for s_ in src],
                   closed=[j for j, s_ in enumerate(src) if s_ and all(i in st["closed"] for i in s_)],
                   frac={j: min(0.95, sum(frac[i] for i in s_ if i in frac)) for j, s_ in enumerate(src)
                         if any(i in frac for i in s_)})
        return [new] + list(old[1:])

    def _ms_tab(self):
        tab = getattr(self.frame, "ms", None)
        return tab if isinstance(tab, MSTab) else None

    def _build_ms_tiles(self):
        """The mass spectrum of each scan event of the MS data (none without MS or PDA data): with several
        events one tile with a plot per event, one below the other."""
        for w in self._ms_tiles():
            try:
                w.Destroy()
            except RuntimeError:
                pass
        self.mscols = []
        self._ms_key = None
        ms = self.ms
        n = min(ms.n_events, MSTab.MAXCOLS) if (ms is not None and self.pda is not None) else 0
        self.ms_tile = PairCard(self.rows, "Mass spectra") if n > 1 else None
        for e in range(n):
            card = PlotCard(self.ms_tile or self.rows, "Mass spectrum", "", mode="zoom")
            col = {"e": e, "spec": None, "spec_card": card, "pins": [], "measure": None, "m1": None, "desc": "",
                   "range": None, "msg": ""}
            card.on_menu = lambda x, y, c=col: self.ms_menu(x, y, c)
            card.on_view = lambda c=col: self.plot_ms(c, keep_view=True)
            card.readout_fmt = lambda x, y: MSTab.READ_FMT % (x, _g(y))
            card.image_name = lambda c=col: self._ms_name(c)
            card.tool_getter = lambda: self.tool
            card.tools_supported = {"label"}
            card.on_tool = lambda tool, kind, x, y, c=col: self.ms_tool(c, tool, x)
            self.mscols.append(col)
        if self.ms_tile is not None:
            self.ms_tile.set_cards([c["spec_card"] for c in self.mscols])
        self._layout_tiles()
        for c in self.mscols:
            self.plot_ms(c)

    def _ms_times(self):
        """(MS time range a, b (a == b: one time), background range or None, delay) for the time or the range
        chosen here: MS time = PDA time + MS detector delay (as Link MS and PDA times)."""
        d = _num(self.delay, 0.0) or 0.0
        if self.avg and self.avg[1] != self.avg[0]:
            bg = self.bg if (self.use_bg.GetValue() and self.bg) else None
            return self.avg[0] + d, self.avg[1] + d, ((bg[0] + d, bg[1] + d) if bg else None), d
        return self.t_sel + d, self.t_sel + d, None, d

    def compute_ms(self):
        """The mass spectra of the time or time range chosen."""
        if not self.mscols or getattr(self, "_ms_hold", False) or self.pda is None or self.t_sel is None:
            return
        tab = self._ms_tab()
        if tab is not None and getattr(tab, "_averaging", False):
            return  # (the reader is busy in the averaging thread of the MS view)
        a, b, bg, d = self._ms_times()
        binw = tab.bin_width() if tab is not None else 0.05
        key = (a, b, bg, binw)
        if key == self._ms_key:
            return
        self._ms_key = key
        ms = self.ms
        later = " (delay %g min)" % d if d else ""
        wx.BeginBusyCursor()
        try:
            for col in self.mscols:
                e = col["e"]
                col.update(spec=None, range=None, desc="", msg="", measure=None, m1=None)
                try:
                    idx = ms.event_scans(e)
                    rt = np.asarray(ms.rt[idx], float)
                    gap = max(2.0 * float(np.median(np.diff(rt))) if len(rt) > 1 else 0.0, 0.02)
                    if not len(rt) or b < rt[0] - gap or a > rt[-1] + gap:
                        col["msg"] = "No MS scans at %s min%s" % (
                            ("%.2f" % a) if a == b else "%.2f to %.2f" % (a, b), later)
                        continue
                    n = 0
                    if a != b:
                        spec, n = ms.average(e, a, b, bg=bg, binw=binw)
                    if n:
                        col.update(spec=spec, range=(a, b), desc="%.2f to %.2f min, %d scans%s%s" % (
                            a, b, n, ", background %.2f to %.2f min subtracted" % bg if bg else "", later))
                    else:  # one time, or a range narrower than the scan interval: the nearest scan
                        t = 0.5 * (a + b)
                        spec, k = ms.scan_spectrum(e, t, binw=binw)
                        j = idx[int(np.argmin(np.abs(rt - t)))]
                        col.update(spec=spec, desc="scan %d at %.3f min%s" % (k + 1, ms.rt[j], later))
                except Exception as ex:
                    _log("PDA view mass spectrum: %s" % ex)
                    col["msg"] = "The mass spectrum could not be read"
        finally:
            wx.EndBusyCursor()
        for col in self.mscols:
            self.plot_ms(col, keep_view=True)

    def plot_ms(self, col, keep_view=False):
        """A mass spectrum tile, drawn as in the Mass spectrometry view."""
        card = col["spec_card"]
        ax = card.ax
        view = ax.get_xlim() if (keep_view and card.full) else None
        card.reset()
        style_axes(ax, "m/z", "Intensity")
        ms, e = self.ms, col["e"]
        title = "Mass spectrum %s" % _pol_short(ms.events[e]["polarity"], e, ms.events) if ms else "Mass spectrum"
        card.side_name = (_pol_short(ms.events[e]["polarity"], e, ms.events) + " MS") if ms else "MS"
        card.side_colour = NEG_COL if (ms and len(self.mscols) > 1 and ms.events[e].get("polarity") == "-") else TRACE
        spec = col["spec"]
        if spec is None or not len(spec):
            ax.text(0.5, 0.5, col.get("msg") or "Click or drag across the chromatogram", transform=ax.transAxes,
                    ha="center", va="center", fontsize=9, color=MUTED)
            card.full = None
            card.set_title(title, "")
            card.draw()
            return
        tab = self._ms_tab()
        sticks = tab is not None and tab.spec_style.GetSelection() == 1
        lo, hi = shared_mz_range(ms, self.mscols if card.pair is not None else [col])
        neg = len(self.mscols) > 1 and ms.events[e].get("polarity") == "-"
        auto = draw_spectrum(card, spec, (lo, hi), sticks, view, 8, MSTab.MZ_FMT,
                             color=NEG_COL if neg else TRACE)
        draw_spec_marks(ax, col, auto, MSTab.PIN_FMT, lambda m: spec_peak_height(spec, m, MSTab.PICK_MIN))
        ib = int(np.argmax(spec[:, 1]))
        sub = col["desc"] + "   \u00b7   base peak m/z " + (MSTab.PIN_FMT % spec[ib, 0]) + " (" + _g(spec[ib, 1]) + ")"
        card.set_title(title, sub)
        card.draw()

    def _ms_name(self, col):
        if self.ms is None or col["spec"] is None:
            return None
        pol = {"+": "pos", "-": "neg"}.get(self.ms.events[col["e"]]["polarity"], "ev%d" % (col["e"] + 1))
        rng = col["range"]
        return "%s_MS_%s_%s" % (self.frame.base_name(), pol, "%.2f-%.2fmin" % rng if rng else
                                "%.2fmin" % self._ms_times()[0])

    def ms_tool(self, col, tool, x):
        if col["spec"] is None:
            return
        pk = spec_nearest_peak(col["spec"], col["spec_card"].ax.get_xlim(), x, MSTab.PICK_MIN)
        if pk is None:
            self.status("Click closer to a peak")
            return
        if spec_click(col, tool, pk, MSTab.PIN_FMT, lambda m1, m2: measure_text(m1, m2, self.ms.adduct_sign(col["e"])),
                      self.status):
            self.plot_ms(col, keep_view=True)

    def _set_ms_style(self, i):
        tab = self._ms_tab()
        if tab is not None:
            tab._set_style(i)  # (as the Mass spectrometry view: its spectra and these)

    def ms_menu(self, x, y, col):
        items = []
        spec = col["spec"]
        pk = spec_nearest_peak(spec, col["spec_card"].ax.get_xlim(), x, MSTab.PICK_MIN) if spec is not None else None
        if pk is not None:
            items += [spec_pin_item(col, pk, MSTab.PIN_FMT, lambda m: self.ms_tool(col, "label", m)), (None, None)]
        items += spec_clear_items(col, lambda: self.plot_ms(col, keep_view=True))
        tab = self._ms_tab()
        if tab is not None:
            items.append(spec_display_item(tab.spec_style.GetSelection() == 1, self._set_ms_style))
        if spec is not None:
            if tab is not None:
                items.append(("Show in the Mass spectrometry view", self.show_in_ms))
            items += [("Export spectrum\u2026", lambda: self.export_ms(col)),
                      ("Copy spectrum data", lambda: copy_spectrum(self, spec))]
        items += [(None, None)] + mz_range_items(self, col["spec_card"])
        return items + self.view_menu()

    def show_in_ms(self):
        """The time or range of these spectra chosen in the Mass spectrometry view, which is shown."""
        tab = self._ms_tab()
        if tab is None or tab.data is None or tab.busy_averaging():
            return
        a, b, bg, d = self._ms_times()
        self.frame._linking = True  # (this view stays as it is)
        try:
            if a != b:
                tab.avg, tab.pick_t = (a, b), None
                tab.fill_range_fields()
                tab.plot_chroms(keep_view=True)
                tab.compute_average()
            else:
                tab.chrom_pick(a, None, None)
        finally:
            self.frame._linking = False
        self.frame.show_page(tab)

    def export_ms(self, col):
        if col["spec"] is None:
            return
        path = self.frame.ask_save("Export mass spectrum", self._ms_name(col) + ".txt", SPEC_WILDCARD)
        if not path:
            return
        try:
            write_spectrum_file(path, col["spec"], self.frame.base_name(), self.ms.events[col["e"]]["polarity"],
                                col["range"])
        except Exception as ex:
            _report_error(self, _write_failed(path), ex)
            return
        self.status("Saved " + path)

    def spec_menu(self, x, y):
        items = []
        if self.uv is not None:
            wl = float(self.uv[0][np.argmin(np.abs(self.uv[0] - x))])
            items.append(("Move the wavelength slider to %.0f nm" % wl, lambda w=wl: self.on_cross(self.t_sel, w, True)))
            items.append(("Add a chromatogram at %.0f nm" % wl, lambda w=wl: self.add_wavelength(w)))
            items.append((None, None))
        items += [("Average a time range\u2026", self.ask_average), (None, None),
                  ("Full view", self.spec_card.reset_view), ("Copy spectrum data", self.copy_uv),
                  ("Export spectrum\u2026", self.on_export_uv)]
        return items + self.view_menu()

    def add_wavelength(self, w):
        vals = _numlist(self.wls.GetValue())
        if not any(abs(v - w) < 0.6 for v in vals):
            vals.append(round(w))
        self.wls.SetValue(", ".join("%g" % v for v in vals))
        self.on_update()

    def copy_uv(self):
        if self.uv is None:
            return
        if _clip("\n".join("%.2f\t%.4f" % (a, b) for a, b in zip(*self.uv))):
            self.status("Spectrum copied")
        else:
            self.status("The clipboard is in use by another program; try again")

    def on_export_uv(self, e=None):
        if self.uv is None:
            self.status("No UV spectrum to export")
            return
        tag = ("%.2fmin" % self.t_sel) if not (self.avg and self.avg[1] != self.avg[0]) else \
            "%.2f-%.2fmin" % self.avg
        path = self.frame.ask_save("Export UV spectrum", "%s_UV_%s.txt" % (self.frame.base_name(), tag),
                                   "Text (*.txt)|*.txt|CSV (*.csv)|*.csv")
        if not path:
            return
        try:
            if path.lower().endswith(".csv"):
                with open(path, "w", encoding="utf-8") as fh:
                    fh.write("Wavelength (nm),Absorbance (mAU)\n")
                    for a, b in zip(*self.uv):
                        fh.write("%.2f,%.4f\n" % (a, b))
            else:
                lcms_pda.write_xy(path, self.uv[0], self.uv[1], "# %s, UV spectrum %s\n# nm\tmAU"
                                  % (self.frame.base_name(), self.uv_desc))
        except Exception as ex:
            _report_error(self, _write_failed(path), ex)
            return
        self.status("Saved " + path)

    def peak_selected(self, p):
        self._clear_range()  # the UV spectrum at the apex of the peak
        self.on_cross(p["rt"], self.wl_sel, True)

    # ---------------------------------------------------------------- map
    def plot_map(self):
        card = self.map_card
        ax = card.ax
        card.reset()
        card.arts = {}
        style_axes(ax, "Time (min)", "Wavelength (nm)")
        if self.pda is None:
            ax.text(0.5, 0.5, "No PDA data", transform=ax.transAxes, ha="center", va="center", fontsize=9,
                    color=MUTED)
            card.full = card.full_y = None
            card.set_title("Wavelength map", "")
            card.draw()
            return
        from matplotlib.colors import PowerNorm, Normalize
        from matplotlib.ticker import ScalarFormatter
        A = self.pda.A
        step = max(1, len(A) // 2400)
        img = np.clip(A[::step].T, 0, None)
        vmax = max(float(np.percentile(A, 99.7)) if A.size else 1.0, 1e-6)
        comp = self.cscale.GetSelection() != 1
        norm = PowerNorm(0.4, vmin=0, vmax=vmax) if comp else Normalize(vmin=0, vmax=vmax)
        t, wl = self.pda.times, self.pda.wavelengths
        w0, w1 = (float(wl[0]), float(wl[-1])) if wl[-1] > wl[0] else (wl[0] - 1.0, wl[0] + 1.0)  # (one wavelength)
        ax.imshow(img, origin="lower", aspect="auto", cmap="viridis", norm=norm,
                  extent=[t[0], t[-1], w0, w1], interpolation="nearest")
        ax.yaxis.set_major_formatter(ScalarFormatter())
        card.full = (float(t[0]), float(t[-1]))
        card.full_y = (w0, w1)
        ax.set_xlim(*card.full)
        ax.set_ylim(*card.full_y)
        card.set_title("Wavelength map", "colour 0 to %.0f mAU%s" % (vmax, ", compressed scale" if comp else ""))
        if self.t_sel is None:
            self.t_sel, self.wl_sel = self._snap(t[len(t) // 2], 280.0)
        card.setup(self.t_sel, self.wl_sel)
        card.draw()

    def map_readout(self, x, y):
        if self.pda is None or x is None or y is None:
            return ""
        i = int(np.argmin(np.abs(self.pda.times - x)))
        j = int(np.argmin(np.abs(self.pda.wavelengths - y)))
        return "%.3f min   %.0f nm   %.1f mAU" % (self.pda.times[i], self.pda.wavelengths[j], self.pda.A[i, j])

    def map_menu(self, x, y):
        if self.pda is None:
            return self.view_menu()
        wl = float(self.pda.wavelengths[np.argmin(np.abs(self.pda.wavelengths - y))])
        return [("Move the sliders here (%.2f min, %.0f nm)" % (x, wl), lambda: self.on_cross(x, y, True)),
                ("Go to a time and wavelength\u2026", self.ask_goto),
                ("Add a chromatogram at %.0f nm" % wl, lambda w=wl: self.add_wavelength(w)),
                (None, None), ("Map and chromatogram settings\u2026", self.ask_settings),
                (None, None), ("Full view", self.map_card.reset_view)] + self.view_menu()

    def extra_chrom_menu(self, x, y, v):
        if self.pda is None:
            return []
        return [(None, None), ("Go to a time and wavelength\u2026", self.ask_goto),
                ("Map and chromatogram settings\u2026", self.ask_settings)] + self.frame.link_menu()

    def average_fields(self):
        self.on_show_uv()

    def ask_goto(self):
        v = ask_form(self, "Go to", [dict(key="t", label="Time", value=self.f_time.GetValue(), unit="min"),
                                     dict(key="wl", label="Wavelength", value=self.f_wl.GetValue(), unit="nm")],
                     ok="Go")
        if v is None:
            return
        self.on_cross(_num_s(v["t"], self.t_sel), _num_s(v["wl"], self.wl_sel), True)

    def ask_settings(self):
        evs = list(self.ms_event.GetItems())
        fields = [dict(key="bw", label="Bandwidth", value=self.bw.GetValue(), unit="nm"),
                  dict(key="cscale", label="Colour scale", kind="choice", choices=["Compressed", "Linear"],
                       value=self.cscale.GetSelection(), width=130),
                  dict(key="wls", label="More wavelengths", value=self.wls.GetValue(), width=120, unit="nm",
                       tip="Separated by commas"),
                  dict(key="maxplot", label="Max plot", kind="check",
                       value=self.maxplot.GetValue()),
                  dict(key="smooth", label="Smoothing", value=str(self.smooth.GetValue()), unit="points"),
                  dict(key="norm", label="Scale each trace to 100 %", kind="check", value=self.norm.GetValue()),
                  dict(key="overlay", label="Overlay the MS trace", kind="check", value=self.overlay.GetValue())]
        if evs:
            fields.append(dict(key="ev", label="MS trace", kind="choice", choices=evs,
                               value=max(0, self.ms_event.GetSelection()), width=150))
        fields.append(dict(key="delay", label="MS detector delay", value=self.delay.GetValue(), unit="min"))
        v = ask_form(self, "Map and chromatogram", fields, ok="Apply", colour=T.GROUP["yellow"])
        if v is None:
            return
        self.bw.SetValue(v["bw"])
        self.cscale.SetSelection(v["cscale"])
        self.wls.SetValue(v["wls"])
        self.maxplot.SetValue(v["maxplot"])
        try:
            self.smooth.SetValue(int(float(v["smooth"] or 0)))
        except ValueError:
            pass
        self.norm.SetValue(v["norm"])
        self.overlay.SetValue(v["overlay"])
        if "ev" in v:
            self.ms_event.SetSelection(v["ev"])
        self.delay.SetValue(v["delay"])
        self.on_update()
        self.plot_map()

    # ------------------------------------------- method presets (method_presets.py)
    def method_get(self):
        """Settings of this view in a method preset (not the times, ranges or peaks of the file)."""
        return dict(wavelength=self.wl_sel if self.pda is not None else getattr(self, "method_wl", None),
                    bandwidth=_num(self.bw, None), colour_scale=("compressed", "linear")[max(0, self.cscale.GetSelection())],
                    more_wavelengths=_numlist(self.wls.GetValue()), max_plot=self.maxplot.GetValue(),
                    smoothing=self.smooth.GetValue(), scale_100=self.norm.GetValue(), overlay_ms=self.overlay.GetValue(),
                    ms_delay=_num(self.delay, None), subtract_background=self.use_bg.GetValue(), **self.integration_get())

    def method_set(self, v):
        """The values given (the others stay) as if chosen in the panel: computed and drawn again."""
        if "colour_scale" in v:
            self.cscale.SetSelection(("compressed", "linear").index(v["colour_scale"]))
        for key, ctrl in (("bandwidth", self.bw), ("ms_delay", self.delay)):
            if key in v:
                ctrl.ChangeValue("%.10g" % v[key])
        if "more_wavelengths" in v:
            self.wls.ChangeValue(", ".join("%g" % w for w in v["more_wavelengths"]))
        for key, cb in (("max_plot", self.maxplot), ("scale_100", self.norm), ("overlay_ms", self.overlay),
                        ("subtract_background", self.use_bg)):
            if key in v:
                cb.SetValue(bool(v[key]))
        if "smoothing" in v:
            self.smooth.SetValue(int(v["smoothing"]))
        self.integration_set(v)
        wl = v.get("wavelength")
        if wl is not None:
            self.method_wl = wl
        if wl is not None and self.pda is not None and self.t_sel is not None:
            self.on_cross(self.t_sel, wl, True)  # the chromatogram at that wavelength and the UV spectrum
        else:
            self.on_update()
            self.compute_uv()
        self.plot_map()


# ==========================================================================
# main window
# ==========================================================================
def show_text(parent, title, text, extra=None):
    """A long text (the help) in a window that scrolls: a message box cannot (it grows taller than the
    screen). Esc or OK closes it. extra: (label, function) of one more button, left of OK."""
    dlg = wx.Dialog(parent, title=title, style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
    tc = wx.TextCtrl(dlg, value=text, style=wx.TE_MULTILINE | wx.TE_READONLY | wx.TE_RICH2 | wx.BORDER_NONE)
    try:
        tc.SetFont(wx.Font(10, wx.FONTFAMILY_DEFAULT, wx.FONTSTYLE_NORMAL, wx.FONTWEIGHT_NORMAL))
    except Exception:
        pass
    ok = wx.Button(dlg, wx.ID_OK, "OK")
    box = wx.BoxSizer(wx.VERTICAL)
    box.Add(tc, 1, wx.EXPAND | wx.ALL, 12)
    row = wx.BoxSizer(wx.HORIZONTAL)
    if extra:
        eb = wx.Button(dlg, wx.ID_ANY, extra[0])
        eb.Bind(wx.EVT_BUTTON, lambda e: (dlg.EndModal(wx.ID_OK), wx.CallAfter(extra[1])))
        row.Add(eb, 0, wx.RIGHT, 8)
    row.Add(ok, 0)
    box.Add(row, 0, wx.ALIGN_RIGHT | wx.RIGHT | wx.BOTTOM, 12)
    dlg.SetSizer(box)
    area = T.display_area(dlg)
    w, h = min(dlg.FromDIP(720), int(area.width * 0.9)), int(area.height * 0.8)
    dlg.SetSize((w, h))
    dlg.CentreOnParent()
    ok.SetDefault()
    ok.SetFocus()
    tc.SetInsertionPoint(0)
    tc.ShowPosition(0)
    dlg.ShowModal()
    dlg.Destroy()


HELP_TEXT = """LCMS Analysis opens LC-MS and HPLC files: Shimadzu .lcd, Agilent .D,
Waters .raw, Thermo .raw, mzML, mzXML and ANDI .cdf. Sciex .wiff: convert
to mzML first (ProteoWizard MSConvert). Most tools are in the
right click menus of the plots; F9 shows or hides the side panel.

Mass spectrometry
  * A TIC or base peak chromatogram per scan event; right click it for
    mass chromatograms.
  * Click a chromatogram for the spectrum at that time, drag to average a
    range, Shift + drag to set the background range.
  * Saturated readings count as the highest normal reading and are marked by red ticks.
  * Right click a chromatogram > Find mass: where the charge states of a
    mass appear in the run. Click a shaded region for its spectrum.

PDA
  * Map: drag the handles on the top and right edges (or use the arrow
    keys) to set the time and the wavelength.
  * MS and PDA times are linked, with the detector delay. Right click a
    chromatogram > Link MS and PDA times to switch this off.

Compare (third view)
  * The chromatograms of every open file, stacked, offset or overlaid.
  * Right click a peak to align or scale the traces or to add a guide
    line; right click > Graph properties for fonts, axes and lines.
  * m/z (tool bar): click a peak of a trace for its mass spectra.

Deconvolution
  * Right click the spectrum > Deconvolute this spectrum (Bayesian or
    maximum entropy). The results appear below the spectrum.
  * The full deconvolution window: right click the spectrum > Open this
    spectrum in the Deconvolute window.
  * Right click the zero charge mass spectrum > Mass shifts.

Tool bar
  * Select, Zoom, Pan, Full view, Background range.
  * Integration: Auto, Drag, Click, Split, Delete, Clear. Areas are in
    signal x seconds.
  * Spectra: Mass chrom., Measure, Label.

Mouse and keys
  * Ctrl + drag zooms a box, double click gives the full view, the wheel
    scrolls the tiles (Ctrl + wheel zooms). Drag the gap below a tile to
    resize it.
  * Ctrl+Z and Ctrl+Y undo and redo; F1 lists every shortcut.
  * Ctrl+C and Ctrl+S copy or save the plot clicked last.

Files, methods and reports
  * Several files can be open: click one in the left panel to show it;
    x or Ctrl+W closes it.
  * Method (top bar) saves every setting under a name and applies it.
  * Report... (Ctrl+R): purity, deconvolution, Supporting Information or
    comparison report, as PDF or Word.
  * Everything saved goes into <file name>_analysis next to the data file.

MS reading uses OpenSZRaw (Apache-2.0)."""


class _Doc(object):
    """One open data file of an analysis window: its own views (tabs), data and
    results. Methods of the window run for a file through this object, so
    they see that file (path, data, views) while the window's widgets (tool
    bar, status bar) stay shared. Views get it as their "frame"."""

    def __init__(self, frame):
        object.__setattr__(self, "_frame", frame)
        object.__setattr__(self, "attrs", {})

    def __getattr__(self, name):
        a = object.__getattribute__(self, "attrs")
        if name in a:
            return a[name]
        fr = object.__getattribute__(self, "_frame")
        for klass in type(fr).__mro__:
            if klass.__module__.split(".")[0] == "wx":
                break
            if name in klass.__dict__:
                v = klass.__dict__[name]
                if isinstance(v, types.FunctionType):
                    return types.MethodType(v, self)
                if isinstance(v, staticmethod):
                    return v.__func__
                if isinstance(v, (classmethod, property)):
                    return getattr(fr, name)
                return v
        return getattr(fr, name)

    def __setattr__(self, name, value):
        object.__getattribute__(self, "attrs")[name] = value


class _VChip(object):
    """Value of a tool bar chip for one file; shown when that file is the
    one on screen."""

    def __init__(self, frame, doc, key):
        self.frame, self.doc, self.key = frame, doc, key

    def SetValue(self, v):
        self.doc.chip_vals[self.key] = v
        if self.frame.active is self.doc:
            self.frame.chips[self.key].SetValue(v)

    def GetValue(self):
        return self.doc.chip_vals.get(self.key, "")

    def __getattr__(self, name):
        return getattr(self.frame.chips[self.key], name)


class FilesPanel(wx.Panel):
    """Left panel of an analysis window, an analysis tree as in Bruker DataAnalysis: the open files (click to
    show one, x or Ctrl+W to close it, + to open more); under the file shown (or one unfolded with its arrow)
    its chromatograms (eye: hide a mass chromatogram or show it again, x: remove it), spectra, PDA traces and
    deconvolution results (eye, x; Up, Down, Space, Delete) and below them the data files of the folder,
    newest first (double click to open). Drag its right edge to make it wider or narrower (the width is kept).
    The views and results tell the panel when they change (tree_changed, watch_results); nothing is polled."""
    HEAD, ROW, LROW, FHEAD, FROW, GRIP = 30, 24, 21, 28, 21, 5  # heights and the grip at the right edge (DIP)
    WIDTH, WIDTH_MIN, WIDTH_MAX = 150, 110, 560
    GROUP_COL = {"Chromatograms": "#0BA064", "Spectra": "#C99A1E", "PDA": "#E0602F", "Deconvolution": "#8E5BD0"}
    BG, BAND, BAND_LINE, EDGE = "#FFFFFF", "#F1F4F8", "#E3E7ED", "#B7CBEE"
    SEL, HOT = "#DCE7FB", "#F1F4F9"

    def __init__(self, parent, frame):
        wx.Panel.__init__(self, parent, style=wx.BORDER_NONE | wx.WANTS_CHARS)
        self.frame = frame
        self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
        self.hover, self.hover_x, self.top, self.fsel = -1, False, 0, None
        self.hover_part = None  # "eye", "x" or "chev" of the item under the mouse
        self.kcur = None  # (deconvolution panel, result) with the keyboard cursor
        self.nsel = None  # key of the trace or spectrum clicked last
        self._acting = False  # a change made here (the cursor stays) rather than elsewhere
        self._drag = None  # (mouse x on the screen, width) while the right edge is dragged
        self._live = None  # LiveDrag while the right edge is dragged
        self._tree_due = False
        self.folder, self.fitems, self._fsig, self._info = "", [], None, {}
        st = settings()
        self.collapsed = bool(st.get("files_collapsed", False))
        try:
            self.width_dip = int(st.get("files_width") or self.WIDTH)
        except (TypeError, ValueError):
            self.width_dip = self.WIDTH
        self.width_dip = min(max(self.width_dip, self.WIDTH_MIN), self.WIDTH_MAX)
        self._apply_width()
        self.Bind(wx.EVT_PAINT, self._paint)
        self.Bind(wx.EVT_SIZE, lambda e: (self._clamp_top(), self.Refresh(), e.Skip()))
        self.Bind(wx.EVT_MOTION, self._motion)
        self.Bind(wx.EVT_LEAVE_WINDOW, self._leave)
        self.Bind(wx.EVT_LEFT_DOWN, self._click)
        self.Bind(wx.EVT_LEFT_UP, self._up)
        self.Bind(wx.EVT_MOUSE_CAPTURE_LOST, lambda e: self._end_drag())
        self.Bind(wx.EVT_LEFT_DCLICK, self._dclick)
        self.Bind(wx.EVT_RIGHT_UP, self._menu)
        self.Bind(wx.EVT_MOUSEWHEEL, self._wheel)
        self.Bind(wx.EVT_KEY_DOWN, self._key)
        self.Bind(wx.EVT_SET_FOCUS, lambda e: (self.Refresh(), e.Skip()))
        self.Bind(wx.EVT_KILL_FOCUS, lambda e: (self.Refresh(), e.Skip()))
        self.timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, lambda e: self.refresh_folder(), self.timer)
        self.timer.Start(15000)  # new runs of a sequence appear by themselves
        self.Bind(wx.EVT_WINDOW_DESTROY, self._destroyed)

    def _head_h(self):
        """Height of the header row (ANALYSES, + and <): as far down as the tool bar of the view shown beside the
        panel plus the same margin again, so that the row and the tool bar share their middle (a fixed height
        left the header higher than the tool bar). DIP HEAD when no tool bar is shown."""
        st = getattr(self, "_strip", None)
        try:
            ok = st is not None and st.IsShownOnScreen()
        except RuntimeError:
            ok = False
        from ms_plot_render import _LIVE
        if not ok and _LIVE["on"] and getattr(self, "_strip_drag", None) == _LIVE["n"]:
            ok, st = True, None   # (none found since this edge drag began: none while it lasts)
        if not ok:
            self._strip_drag = _LIVE["n"] if _LIVE["on"] else None
            st = None
            todo = list(self.frame.GetChildren())
            while todo and st is None:
                w = todo.pop(0)
                try:
                    if isinstance(w, ToolStrip) and w.IsShownOnScreen():
                        st = w
                        break
                    todo.extend(w.GetChildren())
                except RuntimeError:
                    continue
            self._strip = st
        if st is not None:
            try:
                y0 = self.ScreenToClient(st.GetScreenPosition()).y
                h = st.GetSize()[1]
                if 0 <= y0 < self.FromDIP(120) and h > 0:
                    self._band = (y0, h)  # the grey band: where the tool bar is, the same gap above it
                    return 2 * y0 + h
            except RuntimeError:
                pass
        self._band = (0, self.FromDIP(self.HEAD) - self.FromDIP(2))
        return self.FromDIP(self.HEAD)

    def tree_refresh(self):
        """A view changed its traces or spectra (TabBase.tree_changed)."""
        self._tree_due = False
        try:
            self._clamp_top()
            self.Refresh()
        except RuntimeError:  # (the window closed meanwhile)
            pass

    def _destroyed(self, e):
        if e.GetEventObject() is self:
            try:
                self.timer.Stop()
            except Exception:
                pass
        e.Skip()

    # results of the files (deconvolution), shown as a tree
    def watch_results(self, dec):
        """A deconvolution panel (deconv_tab.DeconvPanel) of one of the files
        calls back after its results were added, activated, hidden, shown
        or closed, and when a run starts or ends."""
        if self._results_changed not in dec.listeners:
            dec.listeners.append(self._results_changed)

    def _results_changed(self, dec):
        if self.kcur is not None and (not self._acting or self.kcur[1] not in self.kcur[0].results):
            self.kcur = None  # changed elsewhere (row menu, new run): the cursor goes back to the active one
        self._clamp_top()
        self.Refresh()

    @staticmethod
    def _decs(d):
        """Deconvolution panels of an open file (one per mass spectrometry view)."""
        out = []
        for _, pg in (getattr(d, "attrs", {}).get("pages") or []):
            dec = getattr(pg, "dec", None)
            if dec is not None and hasattr(dec, "results"):
                out.append(dec)
        return out

    def _branches(self, d):
        """(deconvolution panel, result) of a file, newest first."""
        return [(dec, ent) for dec in self._decs(d) for ent in list(dec.results)]

    def _expanded(self, d):
        """The results of a file are listed: the file shown, unless folded
        by hand; the others when unfolded by hand."""
        v = getattr(d, "attrs", {}).get("tree_open")
        return (d is self.frame.active) if v is None else bool(v)

    def toggle_tree(self, d, show=None):
        d.attrs["tree_open"] = (not self._expanded(d)) if show is None else bool(show)
        self._clamp_top()
        self.Refresh()

    def _focus_result(self, d, dec, ent):
        """A click on a result: its file shown, the result in the table and
        its row scrolled into view."""
        if d is not self.frame.active:
            self.frame.activate(d)
        self.kcur = (dec, ent)
        dec.focus(ent)
        self.Refresh()

    # geometry
    def _apply_width(self):
        """Width: the one dragged (at most 45 % of the window), 24 DIP when folded; True when it changed."""
        w = self.FromDIP(24) if self.collapsed else self._px_width(self.width_dip)
        if self.GetMinSize().width == w and self.GetSize()[0] == w:
            return False
        self.SetMinSize(wx.Size(w, -1))
        self.SetSize(wx.Size(w, self.GetSize()[1]))
        return True

    def _px_width(self, wd):
        """Pixels of a width in DIP, at most 45 % of the window."""
        w = self.FromDIP(wd)
        try:
            fw = wx.GetTopLevelParent(self).GetClientSize()[0]
            if fw > 200:
                w = min(w, int(fw * 0.45))
        except Exception:
            pass
        return w

    def rows(self):
        return [d for d in self.frame.docs if d.path or d.loading]

    @staticmethod
    def _sample(d):
        try:
            return ((d.attrs.get("sample_info") or {}).get("sample_name") or "").strip()
        except Exception:
            return ""

    def _lines(self, d):
        """The tree under a file: [("grp", (name, count)), ("node", node), ("res", (panel, result))], a node a dict
        (key, label, tip, page; shown and toggle for an eye, remove for an x)."""
        chrom, spec, pda = [], [], []
        for _, pg in getattr(d, "pages", None) or []:
            try:
                if hasattr(pg, "xics") and hasattr(pg, "cols") and getattr(pg, "data", None) is not None:
                    try:
                        k = pg.kind.GetSelection()
                    except Exception:
                        k = 0
                    for col in pg.cols:
                        e = col["e"]
                        sh = pg.ev_short(e)
                        try:
                            evl = pg.data.event_label(e)
                        except Exception:
                            evl = ""
                        if k in (0, 1):
                            chrom.append({"key": ("main", id(pg), e), "page": pg, "tip": evl,
                                          "label": ("%s %s" % ("TIC" if k == 0 else "BPC", sh)).strip()})
                        for x in pg.xics.get(e, []):
                            hid = (e, xic_key(x)) in pg.xic_hidden
                            chrom.append({"key": ("xic", id(pg), e, xic_key(x)), "page": pg,
                                          "label": (("m/z " + pg.PIN_FMT) % x["mz"] + (" " + sh if sh else "")),
                                          "tip": "%s, window %s" % (evl, xic_window_text(x)), "shown": not hid,
                                          "toggle": lambda pg=pg, e=e, x=x, h=hid: pg.set_xic_hidden(e, x, not h),
                                          "remove": lambda pg=pg, e=e, x=x: pg.remove_xic(e, x)})
                    for col in pg.cols:
                        if col.get("spec") is not None:
                            sh = pg.ev_short(col["e"])
                            desc = col.get("desc") or ""
                            spec.append({"key": ("spec", id(pg), col["e"]), "page": pg, "tip": desc,
                                         "label": ("MS %s %s" % (sh, desc.split(",")[0])).replace("  ", " ").strip()})
                elif hasattr(pg, "uv_desc") and getattr(pg, "pda", None) is not None:
                    if pg.wl_sel is not None:
                        pda.append({"key": ("pdawl", id(pg)), "page": pg, "label": "%.0f nm" % pg.wl_sel,
                                    "tip": "PDA chromatogram at %.1f nm" % pg.wl_sel})
                    if pg.uv is not None:
                        pda.append({"key": ("uv", id(pg)), "page": pg, "tip": pg.uv_desc,
                                    "label": ("UV " + pg.uv_desc.split(",")[0]).strip()})
            except Exception as ex:  # (a view being built or closed)
                _log("files tree: %s" % ex)
        out = []
        for name, nodes in (("Chromatograms", chrom), ("Spectra", spec), ("PDA", pda)):
            if nodes:
                out.append(("grp", (name, len(nodes))))
                out += [("node", n) for n in nodes]
        bs = self._branches(d)
        if bs:
            out.append(("grp", ("Deconvolution", len(bs))))
            out += [("res", b) for b in bs]
        return out

    def _items(self):
        """The open files and, under the expanded ones, their tree: ([(kind, file index, file, payload, rect)],
        y below the last one), kind "row" (a file), "grp", "node" or "res"; scrolled."""
        out = []
        r, lr = self.FromDIP(self.ROW), self.FromDIP(self.LROW)
        x, w = self.FromDIP(4), self.GetClientSize()[0] - self.FromDIP(8 + self.GRIP)
        ind = {"grp": self.FromDIP(20), "node": self.FromDIP(38), "res": self.FromDIP(38)}  # the levels of the tree
        y = self._head_h() - self.top
        rows = self.rows()
        for i, d in enumerate(rows):
            out.append(("row", i, d, None, wx.Rect(x, y, w, r)))
            y += r
            if self._expanded(d) and not d.loading:
                for kind, pl in self._lines(d):
                    out.append((kind, i, d, pl, wx.Rect(x + ind[kind], y, w - ind[kind], lr)))
                    y += lr
                y += self.FromDIP(3)
        if not rows:
            y += r  # the hint takes the place of one file
        return out, y

    def _fhead_y(self, end=None):
        return (self._items()[1] if end is None else end) + self.FromDIP(6)

    def _frow_rect(self, i, fhead=None):
        y0 = (self._fhead_y() if fhead is None else fhead) + self.FromDIP(self.FHEAD)
        r = self.FromDIP(self.FROW)
        return wx.Rect(self.FromDIP(4), y0 + i * r, self.GetClientSize()[0] - self.FromDIP(8 + self.GRIP), r)

    def _x_rect(self, rr):
        s = self.FromDIP(16)
        return wx.Rect(rr.x + rr.width - s - self.FromDIP(3), rr.y + (rr.height - s) // 2, s, s)

    def _chev_rect(self, rr):
        """Fold arrow of a file, at the left of its row."""
        s = self.FromDIP(14)
        return wx.Rect(rr.x + self.FromDIP(1), rr.y + (rr.height - s) // 2, s, s)

    def _eye_rect(self, rr):
        s = self.FromDIP(16)
        return wx.Rect(rr.x, rr.y + (rr.height - s) // 2, s, s)

    def _full_height(self):
        """Height of everything in the panel (unscrolled)."""
        end = self._items()[1] + self.top
        return end + self.FromDIP(6) + self.FromDIP(self.FHEAD) + len(self.fitems) * self.FromDIP(self.FROW) + \
            self.FromDIP(10)

    def _clamp_top(self):
        try:
            room = max(0, self._full_height() - self.GetClientSize()[1])
        except RuntimeError:
            return
        if self.top > room:
            self.top = int(room)

    def _head_buttons(self, fhead=None):
        s = self.FromDIP(20)
        y = (self._head_h() - s) // 2  # in the middle of the header row (the middle of the tool bar beside it)
        w = self.GetClientSize()[0] - (0 if self.collapsed else self.FromDIP(self.GRIP))
        if self.collapsed:
            w = self.GetClientSize()[0]
            return {"expand": wx.Rect((w - s) // 2, y, s, s),
                    "open": wx.Rect((w - s) // 2, y + s + self.FromDIP(6), s, s)}
        fy = (self._fhead_y() if fhead is None else fhead) + self.FromDIP(4)
        return {"collapse": wx.Rect(w - s - self.FromDIP(4), y - self.top, s, s),
                "open": wx.Rect(w - 2 * s - self.FromDIP(6), y - self.top, s, s),
                "folder": wx.Rect(w - s - self.FromDIP(4), fy, s, s)}

    def _on_grip(self, pt):
        return not self.collapsed and pt.x >= self.GetClientSize()[0] - self.FromDIP(self.GRIP)

    def _hit(self, pt):
        if self._on_grip(pt):
            return ("grip", None)
        items, end = self._items() if not self.collapsed else ([], 0)
        fhead = self._fhead_y(end) if not self.collapsed else None
        for k, r in self._head_buttons(fhead).items():
            if r.Contains(pt):
                return ("button", k)
        if self.collapsed:
            return (None, None)
        for kind, i, d, pl, rr in items:
            if not rr.Contains(pt):
                continue
            if kind == "row":
                if self._x_rect(rr).Contains(pt):
                    return ("x", i)
                if self._chev_rect(rr).Contains(pt):
                    return ("chev", i)
                return ("row", i)
            if kind == "grp":
                return ("grp", (i, d))
            if kind == "node":
                if pl.get("toggle") and self._eye_rect(rr).Contains(pt):
                    return ("neye", (i, d, pl))
                if pl.get("remove") and self._x_rect(rr).Contains(pt):
                    return ("nx", (i, d, pl))
                return ("node", (i, d, pl))
            if self._eye_rect(rr).Contains(pt):
                return ("eye", (i, d) + pl)
            if self._x_rect(rr).Contains(pt):
                return ("rx", (i, d) + pl)
            return ("res", (i, d) + pl)
        for i in range(len(self.fitems)):
            if self._frow_rect(i, fhead).Contains(pt):
                return ("frow", i)
        return (None, None)

    # folder of the file shown
    def refresh_folder(self, force=False):
        fr = self.frame
        d = getattr(self, "folder_override", None)
        if not d or not os.path.isdir(d):
            try:
                d = fr.folder() if fr.active is not None else ""
            except Exception:
                d = ""
        if not d or not os.path.isdir(d):
            if self.folder:
                self.folder, self.fitems = "", []
                self._clamp_top()
                self.Refresh()
            return
        try:
            sig = (d, os.path.getmtime(d))
        except OSError:
            sig = (d, 0)
        if sig == self._fsig and not force:
            return
        self._fsig, self.folder = sig, d
        items = []
        try:
            for ent in os.scandir(d):
                p = ent.path
                try:
                    if not fr.browse_accepts(p):
                        continue
                    stt = ent.stat()
                    size = stt.st_size if ent.is_file() else 0
                    items.append({"path": p, "name": ent.name, "mtime": stt.st_mtime, "size": size,
                                  "dir": ent.is_dir()})
                except OSError:
                    continue
        except OSError:
            items = []
        items.sort(key=lambda it: -it["mtime"])
        self.fitems = items[:300]
        self._clamp_top()  # a shorter list: no empty space scrolled into view
        self.Refresh()
        todo = [it["path"] for it in self.fitems if (it["path"], it["mtime"]) not in self._info]
        if todo:
            threading.Thread(target=self._read_infos, args=(todo, {it["path"]: it["mtime"] for it in self.fitems}),
                             daemon=True).start()

    def _read_infos(self, paths, mt):
        for p in paths:
            try:
                import lcms_sources
                info = lcms_sources.read_sample_info(p)
            except Exception:
                info = {}
            self._info[(p, mt.get(p))] = info
        try:
            wx.CallAfter(self.Refresh)
        except Exception:
            pass

    def _finfo(self, it):
        return self._info.get((it["path"], it["mtime"])) or {}

    def _is_open(self, path):
        key = os.path.normcase(os.path.normpath(path))
        return any(os.path.normcase(os.path.normpath(d.path or d.load_path or "")) == key for d in self.frame.docs
                   if d.path or d.load_path)

    # drawing
    def _paint(self, e):
        dc = own_buffered_dc(self)
        dc.SetBackground(wx.Brush(wx.Colour(self.BG)))
        dc.Clear()
        gc = T.crisp(dc)
        ui_font = kept_font
        w, h = self.GetClientSize()
        k = self.FromDIP(10) / 10.0
        items, end = self._items() if not self.collapsed else ([], 0)
        fy = self._fhead_y(end)
        btn = self._head_buttons(fy)

        def header_band(y, hh):
            gc.SetPen(wx.TRANSPARENT_PEN)
            gc.SetBrush(wx.Brush(wx.Colour(self.BAND)))
            gc.DrawRectangle(0, y, w, hh)
            gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(wx.Colour(self.BAND_LINE)).Width(1)))
            gc.StrokeLine(0, y + hh - 0.5, w, y + hh - 0.5)
            if y > 0:
                gc.StrokeLine(0, y + 0.5, w, y + 0.5)

        def edge():
            grip_hot = self.hover == "grip" or self._drag is not None
            gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(wx.Colour(C["accent"] if grip_hot else self.EDGE))
                                   .Width(2 * k if grip_hot else 1)))
            gc.StrokeLine(w - 0.5 - (k if grip_hot else 0), 0, w - 0.5 - (k if grip_hot else 0), h)

        def button(r, glyph, hot, colour=None):
            if hot:
                gc.SetPen(wx.TRANSPARENT_PEN)
                gc.SetBrush(wx.Brush(wx.Colour("#D6DEEB")))
                gc.DrawRoundedRectangle(r.x, r.y, r.width, r.height, 4 * k)
            gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(wx.Colour(colour or C["muted"])).Width(1.5 * k)))
            gc.SetBrush(wx.TRANSPARENT_BRUSH)
            cx, cy, q = r.x + r.width / 2.0, r.y + r.height / 2.0, r.width * 0.22
            if glyph == "+":
                gc.StrokeLine(cx - q, cy, cx + q, cy)
                gc.StrokeLine(cx, cy - q, cx, cy + q)
            elif glyph == "<":
                gc.StrokeLines([(cx + q * 0.5, cy - q), (cx - q * 0.5, cy), (cx + q * 0.5, cy + q)])
            elif glyph == ">":
                gc.StrokeLines([(cx - q * 0.5, cy - q), (cx + q * 0.5, cy), (cx - q * 0.5, cy + q)])
            elif glyph == "v":
                gc.StrokeLines([(cx - q, cy - q * 0.5), (cx, cy + q * 0.5), (cx + q, cy - q * 0.5)])
            elif glyph == "x":
                gc.StrokeLine(cx - q * 0.8, cy - q * 0.8, cx + q * 0.8, cy + q * 0.8)
                gc.StrokeLine(cx - q * 0.8, cy + q * 0.8, cx + q * 0.8, cy - q * 0.8)
            elif glyph == "folder":
                gc.DrawRoundedRectangle(cx - q * 1.1, cy - q * 0.6, q * 2.2, q * 1.5, 1.5 * k)
                gc.StrokeLine(cx - q * 1.1, cy - q * 0.6, cx - q * 0.3, cy - q * 0.6)

        def square(x, cy, colour, s=8.0):
            gc.SetPen(wx.TRANSPARENT_PEN)
            gc.SetBrush(wx.Brush(wx.Colour(colour)))
            gc.DrawRoundedRectangle(x, cy - s * k / 2.0, s * k, s * k, 2 * k)
            return x + s * k + self.FromDIP(5)

        def band(rr, colour):
            gc.SetPen(wx.TRANSPARENT_PEN)
            gc.SetBrush(wx.Brush(wx.Colour(colour)))
            gc.DrawRoundedRectangle(rr.x, rr.y + k, rr.width, rr.height - 2 * k, 4 * k)

        def text(t, x, rr, font, colour, room):
            gc.SetFont(font, wx.Colour(colour))
            t = _ellipsize(gc, t, room)
            tw, th = gc.GetTextExtent(t)
            gc.DrawText(t, x, rr.y + (rr.height - th) / 2.0)
            return tw

        hk = self.hover if isinstance(self.hover, str) else None
        if self.collapsed:
            edge()
            button(btn["expand"], ">", hk == "expand")
            button(btn["open"], "+", hk == "open")
            n = len(self.rows())
            if n:
                gc.SetFont(ui_font(8, 400), wx.Colour(C["accent_text"]))
                t = str(n)
                tw, th = gc.GetTextExtent(t)
                gc.DrawText(t, (w - tw) / 2.0, btn["open"].y + btn["open"].height + self.FromDIP(10))
            return
        hh = self._head_h()
        by, bh = self._band
        header_band(by - self.top, bh)
        gc.SetFont(ui_font(7.8, 400), wx.Colour(C["muted"]))
        th = gc.GetTextExtent("ANALYSES")[1]
        gc.DrawText("ANALYSES", self.FromDIP(9), by + (bh - th) / 2.0 - self.top)
        button(btn["open"], "+", hk == "open")
        button(btn["collapse"], "<", hk == "collapse")
        rows = self.rows()
        if not rows:
            gc.SetFont(ui_font(8, 400), wx.Colour(C["faint"]))
            gc.DrawText("No files open", self.FromDIP(9), hh + self.FromDIP(4) - self.top)
        self._paint_guides(gc, k, items, h)
        for kind, i, d, pl, rr in items:
            if rr.y + rr.height < 0 or rr.y > h:
                continue
            if kind == "row":
                active = d is self.frame.active
                hot = self.hover == ("r", i)
                if active:
                    band(rr, self.SEL)
                elif hot:
                    band(rr, self.HOT)
                cr = self._chev_rect(rr)
                button(cr, "v" if self._expanded(d) else ">", hot and self.hover_part == "chev")
                x = square(cr.x + cr.width + self.FromDIP(2), rr.y + rr.height / 2.0, C["accent"], 9.0)
                end = rr.x + rr.width - self.FromDIP(22 if (active or hot) else 4)
                tw = text(d.file_name, x, rr, ui_font(8.6, 400), C["text"], end - x)  # regular, as Explorer
                if d.loading:
                    sub, col = "reading …", C["accent_text"]
                elif any(getattr(pn, "busy", False) for _, pg in d.pages for pn in getattr(pg, "panels", [])):
                    sub, col = "deconvoluting …", C["accent_text"]
                else:
                    sub, col = self._sample(d), C["faint"]
                room = end - (x + tw + self.FromDIP(6))
                if sub and room > self.FromDIP(24):
                    text(sub, x + tw + self.FromDIP(6), rr, ui_font(7.6, 400), col, room)
                if active or hot:
                    button(self._x_rect(rr), "x", hot and self.hover_x)
            elif kind == "grp":
                name, n = pl
                x = square(rr.x, rr.y + rr.height / 2.0, self.GROUP_COL.get(name, C["muted"]))
                tw = text(name, x, rr, ui_font(8, 400), C["muted"], rr.x + rr.width - x)  # regular: nothing bold in the panel
                if rr.x + rr.width - (x + tw) > self.FromDIP(20):
                    text(str(n), x + tw + self.FromDIP(5), rr, ui_font(7.5, 400), C["faint"], self.FromDIP(30))
            elif kind == "node":
                hot = self.hover == ("n", pl["key"])
                sel = self.nsel == pl["key"] and d is self.frame.active
                if sel:
                    band(rr, self.SEL)
                elif hot:
                    band(rr, self.HOT)
                shown = pl.get("shown", True)
                if pl.get("toggle"):
                    er = self._eye_rect(rr)
                    if hot and self.hover_part == "eye":
                        gc.SetPen(wx.TRANSPARENT_PEN)
                        gc.SetBrush(wx.Brush(wx.Colour("#D6DEEB")))
                        gc.DrawRoundedRectangle(er.x, er.y, er.width, er.height, 4 * k)
                    self._paint_eye(gc, k, er, shown, C["muted"] if shown else "#A0A8B4")
                    x = er.x + er.width + self.FromDIP(4)
                else:
                    x = rr.x + self.FromDIP(20)
                end = rr.x + rr.width - self.FromDIP(22 if (hot and pl.get("remove")) else 4)
                text(pl["label"], x, rr, ui_font(8.3, 400), C["text"] if shown else "#8E97A3",
                     end - x)
                if hot and pl.get("remove"):
                    button(self._x_rect(rr), "x", self.hover_part == "x")
            else:
                self._paint_result(gc, k, rr, d, pl[0], pl[1], button, band, text)
        # folder
        header_band(fy, self.FromDIP(self.FHEAD) - self.FromDIP(2))
        gc.SetFont(ui_font(7.8, 400), wx.Colour(C["muted"]))
        gc.DrawText("FOLDER", self.FromDIP(9), fy + self.FromDIP(8))
        tw = gc.GetTextExtent("FOLDER")[0]
        button(btn["folder"], "folder", hk == "folder")
        gc.SetFont(ui_font(7.6, 400), wx.Colour(C["faint"]))
        fname = os.path.basename(self.folder.rstrip("\\/")) if self.folder else "no folder"
        gc.DrawText(_ellipsize(gc, fname, btn["folder"].x - self.FromDIP(16) - tw - self.FromDIP(9)),
                    self.FromDIP(9) + tw + self.FromDIP(6), fy + self.FromDIP(8))
        for i, it in enumerate(self.fitems):
            rr = self._frow_rect(i, fy)
            if rr.y + rr.height < 0:
                continue
            if rr.y > h:
                break
            opened = self._is_open(it["path"])
            hot = self.hover == ("f", i)
            sel = self.fsel == it["path"]
            if sel:
                band(rr, self.SEL)
            elif hot:
                band(rr, self.HOT)
            x = rr.x + self.FromDIP(6)
            gc.SetFont(ui_font(7.5, 400))
            t = datetime.datetime.fromtimestamp(it["mtime"]).strftime("%d.%m %H:%M")
            dw = gc.GetTextExtent(t)[0]
            show_date = rr.width > self.FromDIP(150)
            room = rr.width - self.FromDIP(12) - (dw + self.FromDIP(6) if show_date else 0)
            text(it["name"], x, rr, ui_font(8.3, 400), C["text"] if opened else C["muted"], room)
            if show_date:
                text(t, rr.x + rr.width - self.FromDIP(6) - dw, rr, ui_font(7.5, 400), C["faint"], dw + 2)
        edge()

    def _paint_guides(self, gc, k, items, h):
        """Lines of the tree: from each file (under its arrow) to its groups, from each group (under its square)
        to its traces, spectra and results."""
        gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(wx.Colour("#C3CCD8")).Width(1)))

        def draw(x, top, ends):
            if not ends:
                return
            y1 = ends[-1][1]
            if y1 < 0 or top > h:
                return
            gc.StrokeLine(x, top, x, y1)
            for xe, ym in ends:
                gc.StrokeLine(x, ym, xe, ym)

        file_x = grp_x = None
        file_top = grp_top = 0
        fends, gends = [], []
        for kind, i, d, pl, rr in items + [("end", None, None, None, None)]:
            if kind in ("row", "end", "grp"):
                draw(grp_x, grp_top, gends)
                gends = []
                grp_x = None
            if kind in ("row", "end"):
                draw(file_x, file_top, fends)
                fends = []
                if kind == "row":
                    cr = self._chev_rect(rr)
                    file_x = int(cr.x + cr.width / 2.0) + 0.5
                    file_top = rr.y + rr.height
            elif kind == "grp":
                ym = rr.y + rr.height / 2.0
                fends.append((rr.x - k, ym))
                grp_x = int(rr.x + 4 * k) + 0.5  # the middle of the group square
                grp_top = ym + 4 * k + k
            elif grp_x is not None:
                gends.append((rr.x - k, rr.y + rr.height / 2.0))

    def _paint_result(self, gc, k, rr, d, dec, ent, button, band, text):
        """One deconvolution result on one line: eye (shown or hidden), method and masses (the time and mass
        range in its tooltip); the one in the table is highlighted."""
        active = ent is dec.active
        hot = self.hover == ("b", id(ent))
        shown = ent.visible
        if active and d is self.frame.active:
            band(rr, self.SEL)
        elif hot:
            band(rr, self.HOT)
        cur = self.kcur[1] if self.kcur is not None else (dec.active if d is self.frame.active else None)
        if cur is ent and self.HasFocus():
            gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(wx.Colour(C["accent"])).Width(1).Style(wx.PENSTYLE_SHORT_DASH)))
            gc.SetBrush(wx.TRANSPARENT_BRUSH)
            gc.DrawRoundedRectangle(rr.x + 1, rr.y + 1, rr.width - 2, rr.height - 2, 4 * k)
        er = self._eye_rect(rr)
        eye_hot = hot and self.hover_part == "eye"
        if eye_hot:
            gc.SetPen(wx.TRANSPARENT_PEN)
            gc.SetBrush(wx.Brush(wx.Colour("#D6DEEB")))
            gc.DrawRoundedRectangle(er.x, er.y, er.width, er.height, 4 * k)
        self._paint_eye(gc, k, er, shown, C["accent_text"] if (shown and (active or eye_hot)) else
                        (C["muted"] if shown else "#A0A8B4"))
        try:
            line1 = dec.tree_text(ent)[0]
        except Exception:
            from ms_brand import method_name
            line1 = method_name(ent.res) or "Result"
        x = er.x + er.width + self.FromDIP(4)
        end = rr.x + rr.width - self.FromDIP(22 if hot else 4)  # the x appears under the mouse only
        poor = False
        try:
            qa = ent.res.get("quality")
            poor = isinstance(qa, dict) and str(qa.get("level", "")).lower() == "poor"
        except Exception:
            pass
        room = end - x - (self.FromDIP(14) if poor else 0)
        tw = text(line1 if shown else "hidden · " + line1, x, rr, ui_font(8.3, 400),
                  C["text"] if shown else "#8E97A3", room)
        if poor:  # warning triangle after the text
            s = self.FromDIP(10)
            self._paint_warning(gc, k, x + tw + self.FromDIP(4), rr.y + (rr.height - s) / 2.0, s)
        if hot:
            button(self._x_rect(rr), "x", self.hover_part == "x")

    @staticmethod
    def _paint_eye(gc, k, r, shown, colour):
        """Eye symbol (drawn): open when the result is shown, struck through
        when it is hidden."""
        col = wx.Colour(colour)
        cx, cy = r.x + r.width / 2.0, r.y + r.height / 2.0
        hw, hh = r.width * 0.36, r.height * 0.21
        path = gc.CreatePath()
        path.MoveToPoint(cx - hw, cy)
        path.AddQuadCurveToPoint(cx, cy - 2 * hh, cx + hw, cy)
        path.AddQuadCurveToPoint(cx, cy + 2 * hh, cx - hw, cy)
        path.CloseSubpath()
        gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(col).Width(1.3 * k).Join(wx.JOIN_ROUND)))
        gc.SetBrush(wx.TRANSPARENT_BRUSH)
        gc.StrokePath(path)
        rp = r.width * 0.1
        if shown:
            gc.SetPen(wx.TRANSPARENT_PEN)
            gc.SetBrush(wx.Brush(col))
            gc.DrawEllipse(cx - rp, cy - rp, 2 * rp, 2 * rp)
        else:
            gc.DrawEllipse(cx - rp, cy - rp, 2 * rp, 2 * rp)
            gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(col).Width(1.3 * k).Cap(wx.CAP_ROUND)))
            gc.StrokeLine(cx - hw * 0.85, cy + hh * 1.7, cx + hw * 0.85, cy - hh * 1.7)

    @staticmethod
    def _paint_warning(gc, k, x, y, s):
        """Small warning triangle (a poor result)."""
        col = wx.Colour(T.GROUP["red"])
        path = gc.CreatePath()
        path.MoveToPoint(x + s / 2.0, y)
        path.AddLineToPoint(x + s, y + s * 0.9)
        path.AddLineToPoint(x, y + s * 0.9)
        path.CloseSubpath()
        gc.SetPen(wx.TRANSPARENT_PEN)
        gc.SetBrush(wx.Brush(col))
        gc.FillPath(path)
        gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(wx.Colour("#FFFFFF")).Width(1.2 * k).Cap(wx.CAP_ROUND)))
        gc.StrokeLine(x + s / 2.0, y + s * 0.32, x + s / 2.0, y + s * 0.58)
        gc.StrokeLine(x + s / 2.0, y + s * 0.74, x + s / 2.0, y + s * 0.75)

    # mouse
    def _motion(self, e):
        if self._drag is not None:
            x0, w0 = self._drag
            sx = self.ClientToScreen(e.GetPosition()).x
            wd = int(round(self.ToDIP(w0 + sx - x0)))
            self._drag_w = min(max(wd, self.WIDTH_MIN), self.WIDTH_MAX)
            if self._live is not None:
                self._live.move(self._drag_w)
            return
        kind, v = self._hit(e.GetPosition())
        self.SetCursor(wx.Cursor(wx.CURSOR_SIZEWE if kind == "grip" else wx.CURSOR_ARROW))
        if kind in ("row", "x", "chev"):
            hv = ("r", v)
        elif kind in ("res", "eye", "rx"):
            hv = ("b", id(v[3]))
        elif kind in ("node", "neye", "nx"):
            hv = ("n", v[2]["key"])
        elif kind == "frow":
            hv = ("f", v)
        elif kind == "grip":
            hv = "grip"
        else:
            hv = v if kind == "button" else -1
        hx = kind == "x"
        part = {"eye": "eye", "rx": "x", "chev": "chev", "neye": "eye", "nx": "x"}.get(kind)
        if hv != self.hover or hx != self.hover_x or part != self.hover_part:
            self.hover, self.hover_x, self.hover_part = hv, hx, part
            tip = ""
            if kind in ("row", "x", "chev"):
                d = self.rows()[v]
                if hx:
                    tip = "Close this file (Ctrl+W)"
                elif kind == "chev":
                    tip = "Fold" if self._expanded(d) else "Unfold"
                else:
                    tip = "\n".join(t for t in (d.path or d.load_path or "", self._sample(d)) if t)
            elif kind in ("res", "eye", "rx"):
                dec, ent = v[2], v[3]
                if kind == "eye":
                    tip = ("Hide this result (Space)" if ent.visible else
                           "Show this result again (Space)")
                elif kind == "rx":
                    tip = "Close this result (Delete)"
                else:
                    try:
                        tip = dec.tree_text(ent)[2]
                    except Exception:
                        tip = ""
            elif kind in ("node", "neye", "nx"):
                nd = v[2]
                if kind == "neye":
                    tip = "Hide" if nd.get("shown", True) else "Show again"
                elif kind == "nx":
                    tip = "Remove"
                else:
                    tip = nd.get("tip") or nd["label"]
            elif kind == "frow":
                it = self.fitems[v]
                info = self._finfo(it)
                lines = [it["name"], datetime.datetime.fromtimestamp(it["mtime"]).strftime("%d.%m.%Y %H:%M") +
                         ("  ·  %.1f MB" % (it["size"] / 1048576.0) if it["size"] else "")]
                for key, lab in (("sample_name", "Sample"), ("sample_id", "Sample ID"), ("operator", "Acquired by"),
                                 ("method_file", "Method"), ("vial", "Vial")):
                    if info.get(key):
                        lines.append("%s: %s" % (lab, info[key]))
                if info.get("acquired"):
                    lines.append("Acquired: %s" % info["acquired"].strftime("%d.%m.%Y %H:%M"))
                tip = "\n".join(lines)
            elif kind == "button":
                tip = {"open": "Open a file (Ctrl+O)", "collapse": "Narrow the file list",
                       "expand": "Show the file list", "folder": "Show another folder"}[v]
            self.SetToolTip(tip)
            self.Refresh()

    def _leave(self, e):
        if self._drag is None:
            self.hover, self.hover_x, self.hover_part = -1, False, None
            self.Refresh()

    def _live_width(self, wd):
        if wd is not None and wd != self.width_dip:
            self.width_dip = wd
            if self._apply_width():
                if live_drag_on():
                    self._place_beside()
                else:
                    self.GetParent().Layout()

    def _place_beside(self):
        """While the edge is dragged: the views beside the panel moved and sized at once (as the sizer of the
        window would place them; it lays them out when the drag ends)."""
        body = self.GetParent()
        sz = body.GetSizer()
        items = [it for it in sz.GetChildren() if it.IsShown()] if sz is not None else []
        rest = items[1:]
        if (len(items) != 2 or items[0].GetWindow() is not self or rest[0].GetWindow() is None
                or rest[0].GetProportion() == 0 or any(it.GetBorder() for it in items)
                or sz.GetOrientation() != wx.HORIZONTAL):
            body.Layout()   # (not the usual panel and view beside it)
            return
        bw, bh = body.GetClientSize()
        fw = self.GetSize()[0]
        deferred = body.BeginRepositioningChildren()   # (moved together, at the end)
        try:
            self.SetSize(0, 0, fw, bh)
            rest[0].GetWindow().SetSize(fw, 0, max(bw - fw, 1), bh)   # the view shown (a file's views or Compare)
        finally:
            if deferred:
                body.EndRepositioningChildren()

    def _end_drag(self):
        if self._drag is None:
            return
        self._drag = None
        if self.HasCapture():
            self.ReleaseMouse()
        live, self._live = self._live, None

        def done():
            self._live_width(self._drag_w)
            if live is not None:
                live.end()
            self.GetParent().Layout()   # (while dragging the views beside were placed directly)
            self.GetParent().Refresh()   # (and only the parts that changed were painted)
        change_at_once(self.GetTopLevelParent(), done)
        T._save({"files_width": int(self.width_dip)})
        self.Refresh()

    def _up(self, e):
        self._end_drag()

    def _click(self, e):
        if self._on_grip(e.GetPosition()):  # the right edge: dragged to a new width
            self._drag = (self.ClientToScreen(e.GetPosition()).x, self.GetSize()[0])
            self._drag_w = self.width_dip
            self.CaptureMouse()
            self._live = LiveDrag(self, self._live_width)
            self.Refresh()
            return
        self._acting = True
        try:
            self._do_click(e)
        finally:
            self._acting = False

    def _show_node(self, d, nd):
        """A trace or spectrum of the tree: its file and its view shown."""
        if d is not self.frame.active:
            self.frame.activate(d)
        try:
            if nd.get("page") is not None:
                self.frame.show_page(nd["page"])
        except Exception as ex:
            _log("files tree: %s" % ex)
        self.nsel = nd["key"]
        self.Refresh()

    def _do_click(self, e):
        kind, v = self._hit(e.GetPosition())
        if kind in ("row", "x", "chev", "res", "eye", "rx", "node", "neye", "nx", "grp"):
            self.SetFocus()  # arrow keys, Space and Delete for the results
        if kind == "button":
            if v == "open":
                self.frame.on_open(None)
            elif v == "folder":
                dlg = wx.DirDialog(self, "Show the data files of a folder", defaultPath=self.folder or "")
                try:
                    if show_open_dialog(dlg) == wx.ID_OK:
                        T._save({self.frame.FOLDER_KEY: dlg.GetPath()})
                        self._list_folder(dlg.GetPath())
                finally:
                    dlg.Destroy()
            else:
                self.collapsed = v == "collapse"
                T._save({"files_collapsed": self.collapsed})
                self._apply_width()
                self.GetParent().Layout()
                self.Refresh()
        elif kind == "row":
            self.frame.activate(self.rows()[v])
        elif kind == "x":
            wx.CallAfter(self.frame.close_doc, self.rows()[v])
        elif kind == "chev":
            self.toggle_tree(self.rows()[v])
        elif kind == "grp":
            if v[1] is not self.frame.active:
                self.frame.activate(v[1])
        elif kind == "node":
            self._show_node(v[1], v[2])
        elif kind == "neye":
            v[2]["toggle"]()
        elif kind == "nx":
            wx.CallAfter(v[2]["remove"])
        elif kind == "res":
            self._focus_result(v[1], v[2], v[3])
        elif kind == "eye":
            self.kcur = (v[2], v[3])
            v[2].set_visible(v[3], not v[3].visible)
        elif kind == "rx":
            wx.CallAfter(v[2].close_result, v[3])
        elif kind == "frow":
            self.fsel = self.fitems[v]["path"]
            self.Refresh()

    def _list_folder(self, d):
        """Show another folder (until a file is opened)."""
        self.folder_override = d
        self.top = 0
        self.refresh_folder(force=True)

    def _dclick(self, e):
        kind, v = self._hit(e.GetPosition())
        if kind == "grip":  # double click on the edge: back to the usual width
            self.width_dip = self.WIDTH
            T._save({"files_width": self.WIDTH})
            if self._apply_width():
                self.GetParent().Layout()
            self.Refresh()
        elif kind == "frow":
            self.frame.load(self.fitems[v]["path"])
        elif kind == "row":
            self.frame.activate(self.rows()[v])
        elif kind in ("res", "chev", "eye", "neye", "node"):
            self._click(e)  # a quick second click: as a click (the arrow and the eye switch back)

    def _menu(self, e):
        kind, v = self._hit(e.GetPosition())
        rows = self.rows()
        items = []
        if kind in ("row", "x", "chev"):
            d = rows[v]
            items = [("Close %s (Ctrl+W)" % d.file_name, lambda: self.frame.close_doc(d))]
            if len(rows) > 1:
                items.append(("Close the other files", lambda: self.frame.close_others(d)))
            items += [("Close every file", self.frame.close_all), (None, None),
                      ("Fold" if self._expanded(d) else "Unfold", lambda: self.toggle_tree(d)), (None, None),
                      ("Show in folder", lambda: _show_in_folder(d.path or d.load_path)),
                      ("Copy the path", lambda: _copy_text(d.path or d.load_path or "")), (None, None)]
        elif kind in ("node", "neye", "nx"):
            d, nd = v[1], v[2]
            items = [("Show it", lambda: self._show_node(d, nd))]
            if nd.get("toggle"):
                items.append(("Hide" if nd.get("shown", True) else "Show again", nd["toggle"]))
            if nd.get("remove"):
                items.append(("Remove", nd["remove"]))
            items.append((None, None))
        elif kind in ("res", "eye", "rx"):
            d, dec, ent = v[1], v[2], v[3]
            self.kcur = (dec, ent)
            items = [("Show it in the table (Enter)", lambda: self._focus_result(d, dec, ent)),
                     ("Hide this result (Space)" if ent.visible else "Show this result again (Space)",
                      lambda: dec.set_visible(ent, not ent.visible)),
                     ("Close this result (Delete)", lambda: dec.close_result(ent)), (None, None)]
            if any(not r.visible for r in dec.results):
                items.append(("Show every result of this file", lambda: dec.set_all_visible(True)))
            if any(r.visible for r in dec.results):
                items.append(("Hide every result of this file", lambda: dec.set_all_visible(False)))
            items += [("Close every result of this file", dec.close), (None, None)]
        elif kind == "frow":
            it = self.fitems[v]
            items = [("Open %s" % it["name"], lambda: self.frame.load(it["path"])),
                     ("Show in folder", lambda: _show_in_folder(it["path"])),
                     ("Copy the path", lambda: _copy_text(it["path"])), (None, None)]
        items += [("Open a file… (Ctrl+O)", lambda: self.frame.on_open(None)),
                  ("Refresh the folder list", lambda: self.refresh_folder(force=True))]
        m = build_menu(self, items)
        self._acting = True
        try:
            popup_menu(self, m)
        finally:
            self._acting = False
        m.Destroy()

    def _wheel(self, e):
        room = max(0, self._full_height() - self.GetClientSize()[1])
        self.top = int(min(room, max(0, self.top - e.GetWheelRotation() / 120.0 * self.FromDIP(self.FROW) * 2)))
        self.Refresh()

    # keyboard (after a click in the list): the results of the files
    def _key(self, e):
        """Up and Down move through the results of the expanded files (each
        one shown in the table), Enter shows the one at the cursor, Space
        hides it or shows it again, Delete closes it, Left and Right fold
        and unfold the results of its file. F9 shows or hides the side
        panel, as in the views."""
        self._acting = True
        try:
            self._do_key(e)
        finally:
            self._acting = False

    def _do_key(self, e):
        code = e.GetKeyCode()
        if code == wx.WXK_F9 and hasattr(self.frame, "set_side"):
            self.frame.set_side(None)
            return
        items = [(d, b) for kind, i, d, b, rr in self._items()[0] if kind == "res"]
        cur = None
        if self.kcur is not None:
            cur = next((n for n, (d, b) in enumerate(items) if b[1] is self.kcur[1]), None)
        if cur is None:  # the result in the table of the file shown
            act = self.frame.active
            cur = next((n for n, (d, b) in enumerate(items) if d is act and b[1] is b[0].active), None)
        if code in (wx.WXK_UP, wx.WXK_DOWN, wx.WXK_NUMPAD_UP, wx.WXK_NUMPAD_DOWN):
            if not items:
                return
            step = -1 if code in (wx.WXK_UP, wx.WXK_NUMPAD_UP) else 1
            n = 0 if cur is None else max(0, min(len(items) - 1, cur + step))
            d, (dec, ent) = items[n]
            self.kcur = (dec, ent)
            if ent.visible:
                self._focus_result(d, dec, ent)
            else:
                self.Refresh()
            self._scroll_to_result(ent)
            return
        if cur is None:
            if code in (wx.WXK_LEFT, wx.WXK_RIGHT) and self.frame.active in self.rows():
                self.toggle_tree(self.frame.active, code == wx.WXK_RIGHT)
                return
            e.Skip()
            return
        d, (dec, ent) = items[cur]
        if code in (wx.WXK_RETURN, wx.WXK_NUMPAD_ENTER):
            self._focus_result(d, dec, ent)
        elif code == wx.WXK_SPACE:
            self.kcur = (dec, ent)
            dec.set_visible(ent, not ent.visible)
        elif code in (wx.WXK_DELETE, wx.WXK_NUMPAD_DELETE, wx.WXK_BACK):
            nxt = items[cur + 1] if cur + 1 < len(items) else (items[cur - 1] if cur > 0 else None)
            self.kcur = (nxt[1][0], nxt[1][1]) if nxt is not None and nxt[0] is d else None
            dec.close_result(ent)
        elif code in (wx.WXK_LEFT, wx.WXK_RIGHT):
            self.toggle_tree(d, code == wx.WXK_RIGHT)
        else:
            e.Skip()

    def _scroll_to_result(self, ent):
        """Scroll the panel so that the result's line is in view."""
        for kind, i, d, b, rr in self._items()[0]:
            if kind == "res" and b[1] is ent:
                h = self.GetClientSize()[1]
                if rr.y < 0:
                    self.top = max(0, self.top + rr.y - self.FromDIP(4))
                elif rr.y + rr.height > h:
                    self.top += rr.y + rr.height - h + self.FromDIP(4)
                self._clamp_top()
                self.Refresh()
                return


def _ellipsize(gc, text, room):
    if gc.GetTextExtent(text)[0] <= room:
        return text
    while text and gc.GetTextExtent(text + "…")[0] > room:
        text = text[:-1]
    return text + "…"


def _show_in_folder(path):
    if not path:
        return
    try:
        if os.name == "nt":
            subprocess.Popen(["explorer", "/select,", os.path.normpath(path)])
        else:
            subprocess.Popen(["xdg-open", os.path.dirname(path)])
    except Exception as ex:
        _log("show in folder: %s" % ex)


def show_open_dialog(dlg):
    """ShowModal of a file or folder dialog of Windows; the log says when it took long to become ready (the
    first one of a session can wait for Windows to list OneDrive or network places, while the window looks
    frozen). The background loading of libraries (T.start_preload) stops while the dialog is open: the
    dialog runs on the main thread, which needs Python's lock for every event of the program's windows
    while it opens, and an import holds that lock for long stretches (loading a library's DLLs), most of
    all on a busy computer."""
    t0 = time.perf_counter()

    def ready():
        dt = time.perf_counter() - t0
        if dt > 1.5:
            _log("the file dialog was ready after %.1f s (Windows was listing the folders)" % dt)
    wx.CallLater(1, ready)
    with T.hold_preload(wait=False):
        return dlg.ShowModal()


def _copy_text(text):
    if wx.TheClipboard.Open():
        try:
            wx.TheClipboard.SetData(wx.TextDataObject(text))
        finally:
            wx.TheClipboard.Close()


def _file_kind(path):
    try:
        name = data_formats.kind_name(data_formats.detect(path)[1])
        if os.path.isdir(path):
            return name + " folder" if name else "Folder"
        ext = os.path.splitext(path)[1].lstrip(".") or "file"
        mb = os.path.getsize(path) / 1048576.0
        return "%s  ·  %.0f MB" % (name or ext, mb)
    except OSError:
        return ""


class ReportDialog(wx.Dialog):
    """Which report (the five kinds; those that need something not done yet
    are greyed with the reason), the details to print, PDF or Word."""

    def __init__(self, parent, avail, cfg, formula_note="", fig_opts=None):
        import ms_report as R
        wx.Dialog.__init__(self, wx.GetTopLevelParent(parent), title="Create a report")
        self.fig_opts = fig_opts or {}
        self.SetBackgroundColour(wx.Colour(C["panel"]))
        # the form scrolls when the window would be taller than the screen (five kinds, long lists of figures)
        self.scroll = wx.ScrolledWindow(self, style=wx.VSCROLL | wx.BORDER_NONE)
        self.scroll.SetBackgroundColour(wx.Colour(C["panel"]))
        self.scroll.SetScrollRate(0, self.FromDIP(12))
        fm = self.form = Form(self.scroll, label_w=120, note_w=470)
        ss = wx.BoxSizer(wx.VERTICAL)
        ss.Add(fm, 1, wx.EXPAND)
        self.scroll.SetSizer(ss)
        fm.section("Report", T.GROUP["blue"])
        self.radios = []
        first = True
        choose = cfg.get("kind")
        for key, name, desc in R.KINDS:
            ok, why = avail.get(key, (False, ""))
            rb = wx.RadioButton(fm, label=name, style=wx.RB_GROUP if first else 0)
            first = False
            rb.SetFont(ui_font(9.5, 700))
            rb.Enable(ok)
            fm.full(rb)
            if not ok:  # a kind not available: one line with the reason
                st = wx.StaticText(fm, label="Not available: %s." % why)
                st.SetFont(ui_font(8.5, 400))
                st.SetForegroundColour(wx.Colour("#9AA3B2"))
                st.Wrap(self.FromDIP(470))
                fm.sizer.Add(st, 0, wx.LEFT, self.FromDIP(22))
            rb.Bind(wx.EVT_RADIOBUTTON, lambda e: self._kind_changed())
            self.radios.append((key, rb))
        oks = [k for k, rb in self.radios if rb.IsEnabled()]
        pick = choose if choose in oks else (oks[0] if oks else None)
        for k, rb in self.radios:
            rb.SetValue(k == pick)
        fm.section("Details", T.GROUP["green"])
        self.compound = fm.text(cfg.get("compound", ""), 300, "Compound name or code")
        fm.row("Compound", self.compound)
        self.operator = fm.text(cfg.get("operator", ""), 220, "Empty: as recorded in the file")
        fm.row("Operator", self.operator)
        self.notes = wx.TextCtrl(fm, value=str(cfg.get("notes") or ""), size=wx.Size(self.FromDIP(300),
                                                                                     self.FromDIP(40)),
                                 style=wx.TE_MULTILINE)
        self.notes.SetFont(ui_font(9, 400))
        fm.row("Notes", self.notes)
        self.assign = fm.text(cfg.get("assignment", ""), 300, "e.g. [M+H]+ 803.2, [M+HCOO]- 847.0")
        self.row_assign = fm.row("Assignment", self.assign)
        self.expected = wx.TextCtrl(fm, value=str(cfg.get("expected") or ""), size=wx.Size(self.FromDIP(300),
                                                                                          self.FromDIP(58)),
                                    style=wx.TE_MULTILINE)
        self.expected.SetFont(ui_font(9, 400))
        self.expected.SetToolTip("One per line: mass and name, e.g. 18363 lysozyme")
        self.row_exp = fm.row("Expected masses", self.expected)
        self.figno = fm.text(cfg.get("fig_no", "1"), 50)
        self.row_fig = fm.row("SI figure number", self.figno)
        if formula_note:
            fm.note(formula_note)
        # figures: one set of check boxes per kind of report, the chosen kind's shown
        fm.section("Figures in the report", T.GROUP["red"])
        saved = cfg.get("figs") if isinstance(cfg.get("figs"), dict) else {}
        self.fig_boxes, self.fig_panels = {}, {}
        for key, name, desc in R.KINDS:
            pn = wx.Panel(fm, style=wx.BORDER_NONE)
            pn.SetBackgroundColour(wx.Colour(C["panel"]))
            vs = wx.BoxSizer(wx.VERTICAL)
            boxes = {}
            got = saved.get(key) if isinstance(saved.get(key), dict) else {}
            for fk, flabel, fdef in self.fig_opts.get(key, []):
                cb = wx.CheckBox(pn, label=flabel)
                cb.SetFont(ui_font(9, 400))
                cb.SetValue(bool(got.get(fk, fdef)))
                vs.Add(cb, 0, wx.TOP, self.FromDIP(4))
                boxes[fk] = cb
            if not boxes:
                st = wx.StaticText(pn, label="No figures")
                st.SetFont(ui_font(8.5, 400))
                st.SetForegroundColour(wx.Colour(C["muted"]))
                vs.Add(st, 0, wx.TOP, self.FromDIP(4))
            pn.SetSizer(vs)
            fm.full(pn)
            self.fig_boxes[key], self.fig_panels[key] = boxes, pn
        fm.section("Output", T.GROUP["yellow"])
        self.fmt = fm.choice(["PDF", "Word (.docx)"], 200)
        self.fmt.SetSelection(1 if cfg.get("fmt") == "docx" else 0)
        fm.row("Format", self.fmt)
        self.open_after = fm.check("Open the report when it is saved", cfg.get("open", True))
        vs = wx.BoxSizer(wx.VERTICAL)
        vs.Add(self.scroll, 1, wx.EXPAND)
        self.ok = flat(self, "Create report…", "primary", handler=lambda e: self.EndModal(wx.ID_OK))
        self.ok.Enable(bool(oks))
        self._btns = dialog_buttons(self, flat(self, "Cancel", handler=lambda e: self.EndModal(wx.ID_CANCEL)),
                                    self.ok)
        vs.Add(self._btns, 0, wx.EXPAND)
        self.SetSizer(vs)
        self._kind_changed()
        self.CentreOnParent()
        self._on_screen()

    def kind(self):
        for k, rb in self.radios:
            if rb.GetValue() and rb.IsEnabled():
                return k
        return None

    def _show_row(self, row, show):
        for it in row.GetChildren():
            w = it.GetWindow()
            if w:
                w.Show(show)

    def _kind_changed(self):
        k = self.kind()
        self._show_row(self.row_exp, k == "deconv")
        self._show_row(self.row_assign, k in ("lcms",))
        self._show_row(self.row_fig, k == "si")
        for key, pn in self.fig_panels.items():
            pn.Show(key == k)
        self.form.Layout()
        self._fit_screen()

    def ShowModal(self):
        wx.CallAfter(self._refit)  # the sizes of the texts are known once the window is on screen
        return wx.Dialog.ShowModal(self)

    def Show(self, show=True):
        r = wx.Dialog.Show(self, show)
        if show:
            self._refit()
        return r

    def _refit(self):
        try:
            self._fit_screen()
            self.CentreOnParent()
            self._on_screen()
        except RuntimeError:  # closed already
            pass

    def _area(self):
        return T.display_area(self)

    def _fit_screen(self):
        """Natural size, but not taller than the screen: then the form scrolls."""
        self.form.InvalidateBestSize()
        best = self.form.GetBestSize()
        bh = self._btns.CalcMin().height if getattr(self, "_btns", None) is not None else self.FromDIP(60)
        h = best.height
        area = self._area()
        frame_h = self.GetSize().height - self.GetClientSize().height  # title bar and borders
        h = min(h, max(self.FromDIP(320), area.height - bh - frame_h - self.FromDIP(12)))
        w = best.width + (max(0, wx.SystemSettings.GetMetric(wx.SYS_VSCROLL_X, self)) if h < best.height else 0)
        self.scroll.SetMinSize(wx.Size(w, h))
        self.scroll.SetVirtualSize(best)
        self.scroll.FitInside()
        self.Layout()
        self.SetClientSize(wx.Size(w, h + bh))
        self._on_screen()

    def _on_screen(self):
        """The whole window (its buttons) on the screen."""
        T.fit_to_screen(self, grow=False, area=self._area())

    def values(self):
        return {"kind": self.kind(), "compound": self.compound.GetValue().strip(),
                "operator": self.operator.GetValue().strip(), "notes": self.notes.GetValue().strip(),
                "assignment": self.assign.GetValue().strip(), "expected": self.expected.GetValue().strip(),
                "fig_no": self.figno.GetValue().strip() or "1",
                "fmt": "docx" if self.fmt.GetSelection() == 1 else "pdf", "open": self.open_after.GetValue(),
                "figs": {fk: cb.GetValue() for fk, cb in self.fig_boxes.get(self.kind(), {}).items()},
                "figs_all": {k: {fk: cb.GetValue() for fk, cb in b.items()} for k, b in self.fig_boxes.items()}}


def draw_pending_plots(win):
    """The plots shown in win that are not drawn yet (a view just shown, a new size): drawn now. Again while
    drawing one asked for another to be drawn (linked axes)."""
    from ms_plot_render import canvas_class
    cls = canvas_class()
    for _ in range(3):
        drawn = 0
        stack = [win]
        while stack:
            w = stack.pop()
            for c in w.GetChildren():
                if not c.IsShown() or isinstance(c, wx.TopLevelWindow):
                    continue
                if isinstance(c, cls):
                    if not getattr(c, "_isDrawn", True):
                        c.draw()
                        drawn += 1
                        card = c.GetParent()   # its live layer now too (it would come one event later)
                        if getattr(card, "live", None) and getattr(card, "_bgcache", None) is not None:
                            card.blit_live()
                else:
                    stack.append(c)
        if not drawn:
            break


def change_at_once(fr, change):
    """change() shows another view or new tile sizes: the window keeps what it shows until the new layout is done,
    its plots drawn and the top bar shows its tools, then shows it at once (no parts drawn one by one, no
    controls seen at the corner of the window while they are made)."""
    if fr.IsFrozen():   # (inside another such change: that one shows it)
        change()
        return
    fr.Freeze()
    try:
        change()
        draw_pending_plots(fr)
        rb = fr.__dict__.get("ribbon")
        if rb is not None:
            try:
                rb.read_state()
            except RuntimeError:
                pass
    finally:
        fr.Thaw()
    try:   # painted now, the window and everything in it, in one pass (Update paints only the window itself)
        if wx.Platform == "__WXMSW__":
            import ctypes
            ctypes.windll.user32.RedrawWindow(ctypes.c_void_p(fr.GetHandle()), None, None,
                                              0x0080 | 0x0100)   # RDW_ALLCHILDREN | RDW_UPDATENOW
        else:
            fr.Update()
    except (RuntimeError, OSError, AttributeError):
        pass


def _view_at_once(method):
    """A method of the window that changes what it shows (another file, the side panel), shown at once."""
    def run(self, *a, **k):
        out = []
        change_at_once(self.window(), lambda: out.append(method(self, *a, **k)))
        return out[0] if out else None
    run.__name__, run.__doc__ = method.__name__, method.__doc__
    return run


class PostrunFrame(wx.Frame):
    """Window with a tool bar (open, tabs, file chips, Deconvolute, help), a
    list of the open files on the left, the views of the file shown and a
    status bar. Every open file keeps its own views and results (a _Doc);
    LCMS Analysis and HRMS Analysis derive from it."""
    TITLE = "Analysis"
    HELP = ""
    SPARE_CHIPS = ("pda", "ms")  # chips of the top bar left out first in a narrow window
    START_HINT = "Open a data file (Ctrl+O) or drop it here"

    def __init__(self, path=None):
        wx.Frame.__init__(self, None, title=self.TITLE)
        self._udp_scaled = True  # sizes here are DPI independent already
        self._home = self._user_display()
        self.docs, self.active = [], None
        self._load_queue = []
        self._load_scheduled = False
        self._side_shown = not T._load().get("side_hidden", True)
        self.SetBackgroundColour(wx.Colour(C["bg"]))
        self._set_icon()
        root = wx.Panel(self, style=wx.BORDER_NONE)
        root.SetBackgroundColour(wx.Colour(C["bg"]))
        tb = T._cls("ModernToolbar")(root, responsive=True)
        self.toolbar = tb
        tb.add(flat(tb, "Open ▾", "ghost", icon="open", tooltip="Open raw data or a saved analysis", height=34,
                    handler=self.on_open_menu), 4)
        self.save_btn=flat(tb, "Save", "ghost", icon="export", tooltip="Save the analysis (Ctrl+Shift+S)", height=34,
                           handler=lambda e: __import__('ms_project').save(self,e))
        tb.add(self.save_btn,4)
        tb.separator()
        import undo  # undo and redo of every edit of the window (undo.py); Undo and Redo in the top bar
        undo.install(self)
        undo.add_buttons(self, tb)
        body = wx.Panel(root, style=wx.BORDER_NONE)
        body.SetBackgroundColour(wx.Colour(C["bg"]))
        self.files = FilesPanel(body, self)
        self.docbook = wx.Simplebook(body)
        first = self._new_doc()
        self.seg = Segmented(tb, [p[0] for p in first.pages], self.on_tab)
        tb.add(self.seg, 4)
        if len(first.pages) < 2:
            self.seg.Hide()  # one view only
        tb.stretch()
        Chip = T._cls("Chip")
        self.chips = {}
        for key, label in self.chip_specs():
            self.chips[key] = Chip(tb, label)
            tb.add(self.chips[key], 3)
        import method_presets  # the Method chip: method presets of LCMS and HRMS Analysis
        method_presets.add_chip(self, tb)
        first.chips = {k: _VChip(self, first, k) for k in self.chips}
        tb.separator()
        tb.add(flat(tb, "Deconvolute…", "primary", icon="deconvolve", height=34,
                    tooltip="Deconvolute the spectrum shown (Ctrl+D)", handler=lambda e: self.deconvolute()), 4)
        tb.add(flat(tb, "Report…", "ghost", icon="report", height=34,
                    tooltip="PDF or Word report (Ctrl+R)", handler=lambda e: self.make_report()), 4)
        self.panel_btn = flat(tb, "Hide panel", "ghost", icon="panel", height=34,
                              tooltip="Show or hide the side panel (F9)", handler=lambda e: self.set_side(None))
        tb.add(self.panel_btn, 4)
        tb.add(flat(tb, "", "ghost", icon="help", tooltip="About " + self.TITLE, height=34, padx=8,
                    handler=self.on_help), 6)
        hs = wx.BoxSizer(wx.HORIZONTAL)
        hs.Add(self.files, 0, wx.EXPAND)
        hs.Add(self.docbook, 1, wx.EXPAND)
        body.SetSizer(hs)
        vs = wx.BoxSizer(wx.VERTICAL)
        vs.Add(tb, 0, wx.EXPAND)
        try:  # a newer MSpektra: the information bar under the top bar (also without the start screen)
            import updater
            vs.Add(updater.bar(root), 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, self.FromDIP(8))
            wx.CallLater(500, updater.start)
        except Exception as ex:
            _log("update: %s" % ex)
        vs.Add(body, 1, wx.EXPAND)
        import ms_ribbon  # the ribbon in place of the top bar (ms_ribbon.py; "ribbon": false for the top bar)
        if ms_ribbon.enabled(self):
            try:
                ms_ribbon.install(self, root, tb, vs)
            except Exception:
                import traceback
                traceback.print_exc()
        root.SetSizer(vs)
        fs = wx.BoxSizer(wx.VERTICAL)
        fs.Add(root, 1, wx.EXPAND)
        self.SetSizer(fs)
        self.CreateStatusBar(2)
        self.SetStatusWidths([-1, screen_dip(self, 330, 240)])
        self.SetStatusText(self.START_HINT, 0)
        # icon buttons of the top bar: their labels go when the bar would need a second row
        self._top_labels = {b: b.GetLabel() for b in tb.GetChildren()
                            if isinstance(b, T._cls("FlatButton")) and b._kind == "ghost" and b._icon and b.GetLabel()}
        fit_rows = tb._fit_rows
        tb._fit_rows = lambda: (self._fit_top_bar(), fit_rows())
        self._screen_fit_pending = False
        self.Bind(wx.EVT_SIZE, self._on_frame_size)
        self.SetDropTarget(_Drop(self))
        self.activate(first)
        self.set_side(self._side_shown, remember=False)
        import shortcuts  # every keyboard shortcut (accelerators, keys of the views, F1 list): shortcuts.py
        shortcuts.install(self)
        import ms_project
        ms_project.install(self)
        self._size_to_display()
        if path:
            wx.CallAfter(self.load, path)

    def Show(self, show=True):
        if not show or getattr(self, "_first_show_done", False):
            return wx.Frame.Show(self, show)
        self._first_show_done = True
        # Lay out controls and tiles at their final sizes before the first paint.
        # Windows otherwise briefly paints new child controls at (0, 0).
        self.Freeze()
        try:
            self._fit_screen()
            self.Layout()
            result = wx.Frame.Show(self, True)
            self._layout_opening_views()
        finally:
            self.Thaw()
        self.Refresh()
        return result

    def _layout_opening_views(self):
        self.Layout()
        for doc in self.docs:
            doc.book.Layout()
            for _, page in doc.pages:
                page.Layout()
                scroller = getattr(page, "scroller", None)
                if scroller is not None:
                    scroller.refit()
                # A canvas resized during startup must paint its current bitmap.
                for card in getattr(page, "all_cards", lambda: [])():
                    card._bgcache = None
                    card.canvas.draw_idle()

    # ------------------------------------------------------ screen size
    def _on_frame_size(self, e):
        e.Skip()
        if not self._screen_fit_pending:
            self._screen_fit_pending = True
            wx.CallAfter(self._fit_screen)

    def _fit_screen(self):
        """Fixed widths (side panels, file list, status bar) for the size of
        the window: compact on a small screen (screen_dip)."""
        if not self:
            return
        self._screen_fit_pending = False
        changed = False
        pages = [p for d in self.docs for _, p in d.pages] + [self.__dict__.get("_compare_tab")]
        for p in pages:
            side = getattr(p, "side", None)
            if isinstance(side, SidePanel) and side.fit_screen():
                changed = True
                p.Layout()
        changed = self.files._apply_width() or changed
        self.SetStatusWidths([-1, screen_dip(self, 330, 240)])
        if changed:
            self.Layout()

    def _fit_top_bar(self):
        """The top bar in one row: in a narrow window the icon buttons lose
        their labels, the last ones first (their tooltips stay), then the
        information chips of SPARE_CHIPS are left out (the status bar
        has the same)."""
        tb = self.toolbar
        if not tb or not self._top_labels or self.__dict__.get("ribbon"):  # (the ribbon: no top bar to fit)
            return
        room = tb.GetClientSize().width
        btns = [b for b in self._top_labels if b.IsShown()]
        spare = [self.chips[k] for k in self.SPARE_CHIPS if k in getattr(self, "chips", {})]
        need = sum(it.CalcMin().width for it in tb.sizer.GetChildren() if it.IsShown())
        need += sum(self._label_width(b, self._top_labels[b]) for b in btns if not b.GetLabel())
        for c in spare:  # left out here before: the room it would take
            if c._value and not c.IsShown():
                need += tb.sizer.GetItem(c).CalcMin().width
        short, drop = set(), set()
        for b in reversed(btns):
            if need <= room:
                break
            need -= self._label_width(b, self._top_labels[b])
            short.add(b)
        for c in spare:
            if need <= room:
                break
            if c._value:
                need -= tb.sizer.GetItem(c).CalcMin().width
                drop.add(c)
        changed = False
        for b, label in self._top_labels.items():
            want = "" if b in short else label
            if b.GetLabel() != want:
                b.SetLabel(want)
                changed = True
        for c in spare:
            show = bool(c._value) and c not in drop
            if c.IsShown() != show:
                c.Show(show)
                changed = True
        if changed:
            tb.sizer.Layout()

    @staticmethod
    def _label_width(b, label):
        dc = wx.ClientDC(b)
        dc.SetFont(b.GetFont())
        return dc.GetTextExtent(label)[0] + b.FromDIP(7)

    def _set_top_label(self, b, label):
        self._top_labels[b] = label
        if b.GetLabel():  # (shown only as an icon: the bar decides when it is laid out)
            b.SetLabel(label)

    def __getattr__(self, name):
        # data, views and results belong to the file shown
        d = self.__dict__.get("active")
        if d is not None:
            a = object.__getattribute__(d, "attrs")
            if name in a:
                return a[name]
        raise AttributeError(name)

    # ------------------------------------------------------------ override
    def build_tabs(self, book):
        return []

    def chip_specs(self):
        return [("file", "File")]

    def accepts(self, path):
        return False

    def browse_accepts(self, path):
        """Files listed in the folder part of the file panel."""
        return self.accepts(path)

    def read_data(self, path, progress):
        """Runs in a worker thread: returns (dict of data, list of errors)."""
        return {}, ["not implemented"]

    def data_loaded(self, path, data, errors):
        pass

    def help_text(self):
        return self.HELP

    def title_text(self):
        return self.TITLE

    def write_spectrum(self, path, spec):
        lcms_data.write_spectrum_txt(path, spec)

    # ------------------------------------------------------------- files
    def window(self):
        """The window itself (also when called for one file)."""
        return object.__getattribute__(self, "_frame") if isinstance(self, _Doc) else self

    def _doc(self):
        return self if isinstance(self, _Doc) else self.active

    def _is_active(self):
        return self.window().active is self._doc()

    def _new_doc(self):
        fr = self.window()
        doc = _Doc(fr)
        doc.path, doc.load_path, doc.loading = None, None, False
        doc.ms_data = None
        doc.file_name, doc.file_kind = "", ""
        doc.doc_title, doc.status1, doc.chip_vals = fr.TITLE, "", {}
        doc.chips = {k: _VChip(fr, doc, k) for k in getattr(fr, "chips", {})}
        doc.book = wx.Simplebook(fr.docbook)
        doc.pages = doc.build_tabs(doc.book)
        for label, page in doc.pages:
            doc.book.AddPage(page, label)
        doc.seg_enabled = [True] * len(doc.pages)
        fr.docbook.AddPage(doc.book, "")
        fr.docs.append(doc)
        for _, pg in doc.pages:
            if hasattr(pg, "side"):
                pg.side.Show(bool(fr._side_shown))
                pg.Layout()
        return doc

    @_view_at_once
    def activate(self, doc):
        fr = self.window()
        if doc is None or doc not in fr.docs:
            return
        fr.active = doc
        i = fr.docbook.FindPage(doc.book)
        if i >= 0 and fr.docbook.GetSelection() != i:
            fr.docbook.ChangeSelection(i)
        for k, chip in getattr(fr, "chips", {}).items():
            chip.SetValue(doc.chip_vals.get(k, ""))
        seg = getattr(fr, "seg", None)
        if seg is not None:
            for i, en in enumerate(doc.seg_enabled):
                seg.enable_segment(i, en)
            seg.set_active(max(0, doc.book.GetSelection()))
        fr.SetTitle(doc.doc_title)
        try:
            fr.files.refresh_folder()
        except Exception:
            pass
        if fr.GetStatusBar():
            fr.SetStatusText(doc.status1, 1)
            if doc.path:
                fr.SetStatusText(doc.file_name, 0)
        fr.files.Refresh()
        fr.Layout()

    def next_doc(self, step=1):
        fr = self.window()
        rows = fr.files.rows()
        if len(rows) < 2 or fr.active not in rows:
            return
        fr.activate(rows[(rows.index(fr.active) + step) % len(rows)])

    def close_doc(self, doc):
        """Close one file (its views, results and data); the window stays."""
        fr = self.window()
        if doc is None or doc not in fr.docs:
            return
        if not doc.path and not doc.loading:
            return  # the empty start view
        if any(getattr(pg, "_averaging", False) for _, pg in doc.pages):
            # the averaging thread still reads the file: closing it now would free
            # the reader under it
            fr.status_msg("%s is still being averaged; close it afterwards" % doc.file_name)
            return
        for _, pg in doc.pages:
            for pn in getattr(pg, "panels", []):
                if getattr(pn, "busy", False) and hasattr(pn, "cancel"):
                    try:
                        pn.cancel()
                    except Exception:
                        pass
        for key in ("ms_data", "pda_data"):
            dd = doc.attrs.get(key)
            if dd is not None and hasattr(dd, "close"):
                try:
                    dd.close()
                except Exception:
                    pass
        i = fr.docs.index(doc)
        fr.docs.remove(doc)
        import undo
        undo.file_closed(fr, doc)  # its undo steps go with it
        if fr.active is doc:
            if not any(d.path or d.loading for d in fr.docs):
                nd = next((d for d in fr.docs if not d.path and not d.loading), None) or fr._new_doc()
            else:
                shown = [d for d in fr.docs if d.path or d.loading]
                nd = shown[min(i, len(shown) - 1)]
            fr.activate(nd)
            if not nd.path:
                fr.SetStatusText(fr.START_HINT, 0)
        k = fr.docbook.FindPage(doc.book)
        if k >= 0:
            fr.docbook.DeletePage(k)
        name = doc.file_name
        doc.attrs.clear()
        doc.attrs.update(path=None, loading=False, pages=[], file_name=name)
        fr.files.Refresh()
        fr.SetStatusText("Closed " + name if name else fr.START_HINT, 0)
        import gc
        gc.collect()

    def close_others(self, keep):
        fr = self.window()
        for d in [d for d in fr.docs if d is not keep and (d.path or d.loading)]:
            fr.close_doc(d)
        fr.activate(keep)

    def close_all(self):
        fr = self.window()
        for d in [d for d in fr.docs if d.path or d.loading]:
            fr.close_doc(d)

    def enable_segment(self, i, flag):
        """View i (e.g. PDA) available for this file or not."""
        doc = self._doc()
        doc.seg_enabled[i] = bool(flag)
        if self._is_active():
            self.window().seg.enable_segment(i, bool(flag))

    def set_file_status(self, text):
        doc = self._doc()
        doc.status1 = text
        if self._is_active():
            self.window().SetStatusText(text, 1)

    # ----------------------------------------------------------------- misc
    def _set_icon(self):
        try:
            T.set_app_icon(self)
            return
        except Exception:
            pass
        try:
            import unidec
            ico = os.path.join(os.path.dirname(unidec.__file__), "bin", "logo.ico")
            if os.path.isfile(ico):
                self.SetIcon(wx.Icon(ico, wx.BITMAP_TYPE_ICO))
        except Exception:
            pass

    def _user_display(self):
        """Work area of the display the user works on: that of the window in front (the start screen whose tile
        was clicked, or the other analysis window), else the one under the mouse. The new window is moved there
        before its contents are made, so they are made at the scale of that display: Windows puts a new window
        without a parent on the main display (a laptop at 200 %), and moving it to the other display (a monitor
        at 100 %) afterwards left sizes of the old scale."""
        area = None
        try:
            act = wx.GetActiveWindow()
            top = wx.GetTopLevelParent(act) if act is not None else None
            if top is not None and top is not self and top.IsShown() and not top.IsIconized():
                area = T.display_area(top, own_first=True)
        except (RuntimeError, AttributeError):
            area = None
        if area is None:
            n = wx.Display.GetFromPoint(wx.GetMousePosition())
            area = wx.Display(n if 0 <= n < wx.Display.GetCount() else 0).GetClientArea()
        self.SetPosition(wx.Point(area.x + area.width // 8, area.y + area.height // 8))
        return area

    def _size_to_display(self):
        area = getattr(self, "_home", None) or T.display_area(self)
        w = min(self.FromDIP(1480), int(area.width * 0.94))
        h = min(self.FromDIP(940), int(area.height * 0.94))
        self.SetSize(wx.Rect(area.x + (area.width - w) // 2, area.y + (area.height - h) // 2, w, h))
        self.SetMinSize(wx.Size(min(self.FromDIP(1000), area.width), min(self.FromDIP(680), area.height)))
        # opens maximised (the size above is the one it gets when restored);
        # "start_maximized": false in the settings file turns this off
        if settings().get("start_maximized", True) or area.width < self.FromDIP(1200) or \
                area.height < self.FromDIP(760):
            self.Maximize()

    def base_name(self):
        if not self.path:
            return "analysis"
        b = os.path.basename(self.path.rstrip("\\/"))
        return os.path.splitext(b)[0]

    def folder(self):
        if self.path:
            return os.path.dirname(self.path.rstrip("\\/"))
        return T._load().get(self.FOLDER_KEY) or T._load().get("last_folder") or ""

    FOLDER_KEY = "lcms_folder"

    def out_dir(self):
        """One folder per data file for everything saved from it:
        <data folder>\\<file name>_analysis (created when first needed)."""
        doc=self._doc()
        saved_folder=getattr(doc,'attrs',{}).get('analysis_folder') if doc else None
        if saved_folder:
            os.makedirs(saved_folder,exist_ok=True)
            return saved_folder
        if not self.path:
            d = os.path.join(os.path.expanduser("~"), APP)
            try:
                os.makedirs(d, exist_ok=True)
                return d
            except OSError:
                return self.folder()
        cands = [os.path.join(self.folder(), self.base_name() + OUT_SUFFIX),
                 os.path.join(os.path.expanduser("~"), APP, self.base_name())]
        for d in cands:
            try:
                os.makedirs(d, exist_ok=True)
                test = os.path.join(d, ".write_test")
                open(test, "w").close()
                os.remove(test)
                return d
            except OSError:
                continue
        return self.folder()

    def ask_save(self, title, name, wildcard):
        dlg = wx.FileDialog(self.window(), title, defaultDir=self.out_dir(), defaultFile=name, wildcard=wildcard,
                            style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
        try:
            if dlg.ShowModal() != wx.ID_OK:
                return None
            path = dlg.GetPath()
            ext = wildcard.split("|")[1 + 2 * dlg.GetFilterIndex()].replace("*", "").split(";")[0]
            if ext and not path.lower().endswith(ext.lower()):
                path += ext
            return path
        finally:
            dlg.Destroy()

    def page_index(self, page):
        for i, (_, p) in enumerate(self.pages):
            if p is page:
                return i
        return -1

    def show_page(self, page):
        i = self.page_index(page)
        if i >= 0:
            self.book.SetSelection(i)
            if self._is_active():
                self.window().seg.set_active(i)

    def on_tab(self, i):
        change_at_once(self.window(), lambda: self.book.SetSelection(i))

    def on_help(self, e=None):
        ver = getattr(T, "APP_VERSION", "")
        text = "MSpektra %s\n\n%s" % (ver, self.help_text())
        from ms_brand import CREDITS
        extra = None
        try:
            import updater
            prev = updater.previous()
            if prev:  # the version before the last update is kept: back to it
                extra = ("Back to MSpektra %s" % prev, updater.back_to_previous)
        except Exception:
            pass
        show_text(self.window(), "About %s (MSpektra %s)" % (self.title_text(), ver), text+"\n\n"+CREDITS,
                  extra=extra)


    # ---------------------------------------------------------------- report
    def make_report(self, kind=None):
        """Report window: choose the kind of report and the details, then save
        it as PDF or Word next to the data (in the _analysis folder)."""
        import ms_report as R
        fr = self.window()
        if not fr.active or not fr.path:
            fr.status_msg("Open a data file first")
            return
        try:
            avail = R.availability(fr)
        except Exception as ex:
            _report_error(fr, "The report could not be prepared", ex)
            return
        cfg = T._load().get("report_cfg")
        cfg = dict(cfg) if isinstance(cfg, dict) else {}  # (a settings file edited by hand)
        si = {}
        try:
            si = R.sample_info(fr)
        except Exception:
            pass
        doc = fr.active
        cfg["compound"] = doc.attrs.get("report_compound", "")
        cfg["notes"] = doc.attrs.get("report_notes", "")
        cfg["assignment"] = doc.attrs.get("report_assignment", "")
        if not cfg.get("operator"):
            cfg["operator"] = si.get("operator", "") if si.get("operator") != "System Administrator" else ""
        if kind:
            cfg["kind"] = kind
        elif R.compare_shown(fr) and avail.get("compare", (False,))[0]:
            cfg["kind"] = "compare"  # the Compare view on screen: its report
        elif cfg.get("kind") == "compare":
            cfg["kind"] = "lcms"  # a file's view on screen: a report of that file
        note = ""
        tab = getattr(fr, "ms", None)
        if tab is not None and not getattr(fr, "pda", None):
            col = R._spec_col(tab) if getattr(fr, "ms_data", None) is not None else None
            fi = R._formula_info(col) if col else None
            note = ("Exact mass in the report: %s." % fi["name"]) if fi else \
                "No formula checked (right click the spectrum > Check a formula)."
        fig_opts = {}
        for k in avail:
            try:
                fig_opts[k] = R.figure_options(fr, k)
            except Exception as ex:
                _log("figure options (%s): %s" % (k, ex))
        dlg = ReportDialog(fr, avail, cfg, note, fig_opts)
        try:
            if dlg.ShowModal() != wx.ID_OK:
                return
            v = dlg.values()
        finally:
            dlg.Destroy()
        if not v["kind"]:
            return
        doc.report_compound, doc.report_notes, doc.report_assignment = v["compound"], v["notes"], v["assignment"]
        figs = dict(cfg.get("figs") if isinstance(cfg.get("figs"), dict) else {})
        for k, d in v.get("figs_all", {}).items():
            merged = dict(figs.get(k) if isinstance(figs.get(k), dict) else {})
            merged.update(d)  # choices of figures this file does not have are kept
            figs[k] = merged
        T._save({"report_cfg": {"kind": v["kind"], "operator": v["operator"], "expected": v["expected"],
                                "fmt": v["fmt"], "open": v["open"], "fig_no": v["fig_no"], "figs": figs}})
        ext = ".docx" if v["fmt"] == "docx" else ".pdf"
        tag = {"hrms": "HRMS_report", "lcms": "LCMS_purity_report", "deconv": "deconvolution_report",
               "si": "SI_page", "compare": "comparison_report"}[v["kind"]]
        wild = "Word document (*.docx)|*.docx" if ext == ".docx" else "PDF (*.pdf)|*.pdf"
        if v["kind"] == "compare":
            # several runs: saved next to them (the folder of the first run), as the exports of the Compare view
            ctab, drawn = R._compare_drawn(fr)
            first = drawn[0]["entry"]["doc"].attrs.get("path") if drawn else fr.path
            try:
                name = ctab.image_name()  # e.g. comparison_254nm
            except Exception:
                name = "comparison"
            dlg = wx.FileDialog(fr, "Save the report", defaultDir=os.path.dirname((first or "").rstrip("\\/")),
                                defaultFile=safe_file_name(name + "_report") + ext, wildcard=wild,
                                style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)
            try:
                path = dlg.GetPath() if dlg.ShowModal() == wx.ID_OK else None
            finally:
                dlg.Destroy()
            if path and not path.lower().endswith(ext):
                path += ext
        else:
            path = fr.ask_save("Save the report", safe_file_name("%s_%s" % (fr.base_name(), tag)) + ext, wild)
        if not path:
            return
        wx.BeginBusyCursor()
        fr.status_msg("Creating the report \u2026")
        try:
            wx.Yield()
            R.create(fr, v["kind"], v, path, v["fmt"])
        except Exception as ex:
            _report_error(fr, _write_failed(path) if isinstance(ex, PermissionError) else
                          "The report could not be created", ex)
            return
        finally:
            try:
                wx.EndBusyCursor()
            except Exception:
                pass
        fr.status_msg("Report saved: " + path)
        toast(fr, "Report saved: %s\nin %s" % (os.path.basename(path), os.path.dirname(path)))
        if v["open"]:
            try:
                if os.name == "nt":
                    os.startfile(path)
                else:
                    subprocess.Popen(["xdg-open", path])
            except Exception as ex:
                _log("open report: %s" % ex)

    # --------------------------------------------------------- deconvolution
    def deconv_inputs(self):
        return self.ms.deconv_inputs() if getattr(self, "ms", None) else []

    def deconvolute(self, e=None):
        """Deconvolute the spectrum of scan event e (or the one chosen in the
        panel) in the mass spectrometry view."""
        if not getattr(self, "ms", None) or not self.ms.data:
            self.status_msg("Open a data file first")
            return
        self.show_page(self.ms)
        self.ms.on_deconvolve(e)

    def status_msg(self, text):
        self.SetStatusText(text, 0)

    def link_time(self, src, t=None, rng=None):
        """Time or time range chosen in one view, shown in the others (LCMS
        Analysis: MS and PDA). Nothing here."""
        pass

    def link_menu(self):
        return []

    @_view_at_once
    def set_side(self, show=None, remember=True):
        """Side panel of every view (of every open file) shown or hidden
        (None: switch)."""
        fr = self.window()
        pages = [p for d in fr.docs for _, p in d.pages if hasattr(p, "side")]
        if not pages:
            return
        if show is None:
            show = not fr._side_shown
        fr._side_shown = bool(show)
        for p in pages:
            p.side.Show(bool(show))
            p.Layout()
        fr._set_top_label(fr.panel_btn, "Hide panel" if show else "Show panel")
        fr.toolbar.Layout()
        if remember:
            T._save({"side_hidden": not show})
            fr.status_msg("Side panel shown" if show else
                          "Side panel hidden (F9 shows it)")

    # ------------------------------------------------------------- open
    def on_open(self, e=None):
        pass

    def raw_open_items(self):
        """(label, function) of the raw data entries of the Open menu."""
        if hasattr(self, "open_mzml"):  # HRMS Analysis: .d folders and mzML files, in this menu (no second one)
            return [("Open a data file (mzML, Thermo .raw)…", self.open_mzml),
                    ("Open a data folder (Bruker or Agilent .d)…", self.open_d),
                    ("Open a data folder (Waters .raw)…", self.open_waters)]
        return [("Open raw data…", lambda: self.on_open(None))]

    def on_open_menu(self, e=None):
        """Open ▾ in the top bar: raw data, an analysis folder or a saved analysis (Ctrl+O: raw data only)."""
        import ms_project
        fr = self.window()
        items = fr.raw_open_items() + [(None, None), ("Open a saved analysis…", lambda: ms_project.open_file(fr))]
        m = build_menu(fr, [(label, (lambda f=f: wx.CallAfter(f)) if f else None) for label, f in items])
        btn = e.GetEventObject() if e is not None else None
        if isinstance(btn, wx.Window) and not isinstance(btn, wx.Frame):
            popup_menu(btn, m, wx.Point(0, btn.GetSize()[1]))
        else:
            popup_menu(fr, m)
        m.Destroy()

    def load(self, path):
        """Queue file setup so a multi-file drop yields between workspaces."""
        fr = self.window()
        path = os.path.abspath(path)
        key = os.path.normcase(os.path.normpath(path.rstrip("\\/")))
        if not any(os.path.normcase(os.path.normpath(p.rstrip("\\/"))) == key for p in fr._load_queue):
            fr._load_queue.append(path)
        if not fr._load_scheduled:
            fr._load_scheduled = True
            wx.CallLater(20, fr._drain_load_queue)

    def _drain_load_queue(self):
        if not self:
            return
        if not self._load_queue:
            self._load_scheduled = False
            return
        if getattr(self, "_asking", False):  # a "Saved analysis" question (ms_project) or "Could not read" is open
            wx.CallLater(100, self._drain_load_queue)
            return
        path = self._load_queue.pop(0)
        try:
            self._start_load(path)
        finally:
            if self and self._load_queue:
                wx.CallLater(40, self._drain_load_queue)
            elif self:
                self._load_scheduled = False

    def _start_load(self, path):
        """Open a file: in a new file slot (the empty start view is used for
        the first one); a file that is open already is just shown."""
        fr = self.window()
        path = os.path.abspath(path)
        key = os.path.normcase(os.path.normpath(path.rstrip("\\/")))
        for d in fr.docs:
            p = d.path or d.load_path
            if p and os.path.normcase(os.path.normpath(p.rstrip("\\/"))) == key:
                fr.activate(d)
                return
        cur = fr.active
        doc = cur if (cur is not None and not cur.path and not cur.loading) else fr._new_doc()
        name = os.path.basename(path.rstrip("\\/"))
        doc.loading, doc.load_path = True, path
        doc.file_name, doc.file_kind = name, _file_kind(path)
        fr.activate(doc)
        fr.SetStatusText("Reading %s …" % name, 0)
        import method_presets
        method_presets.file_opening(fr, doc)  # the default method of the window, if one is set
        wx.BeginBusyCursor()

        def progress(i, n):
            if fr.active is not doc:
                return
            if n:
                wx.CallAfter(fr.SetStatusText, "Reading %s: scan %d of %d" % (name, i, n), 0)
            else:
                wx.CallAfter(fr.SetStatusText, "Reading %s: scan %d" % (name, i), 0)

        t0 = time.perf_counter()
        who = type(fr).__name__.replace("Frame", " Analysis")

        def work():
            try:
                data, errors = doc.read_data(path, progress)
            except Exception as ex:
                data, errors = {}, [str(ex)]
            print("[%s] read %s in %.1f s" % (who, name, time.perf_counter() - t0))  # timing in the log file
            wx.CallAfter(fr._loaded, doc, path, data, errors)

        threading.Thread(target=work, daemon=True).start()

    def _loaded(self, doc, path, data, errors):
        fr = self.window()
        if fr and getattr(fr, "_asking", False):
            # a message about another file is open: its message loop runs this CallAfter, so a second "Could
            # not read" box opened on top of it (and the window changed behind it); this file comes after it
            wx.CallLater(100, fr._loaded, doc, path, data, errors)
            return
        try:
            wx.EndBusyCursor()
        except Exception:
            pass
        name = os.path.basename(path.rstrip("\\/"))
        if not fr or doc not in fr.docs:  # closed (or the window closed) while it was read
            for v in data.values():
                if v is not None and hasattr(v, "close"):
                    try:
                        v.close()
                    except Exception:
                        pass
            return
        doc.loading = False
        if not any(v is not None for v in data.values()):
            fr._asking = True  # (the next files wait: _loaded, _drain_load_queue)
            try:
                wx.MessageBox("Could not read %s\n\n%s" % (path, "\n".join(errors)), fr.TITLE, wx.ICON_ERROR)
            finally:
                fr._asking = False
            doc.load_path = None
            fr.SetStatusText("Could not read " + name, 0)
            if any(d.path for d in fr.docs):
                doc.path = "x"  # let close_doc remove this slot
                fr.close_doc(doc)
            else:
                doc.file_name = doc.file_kind = ""
                fr.files.Refresh()
            return
        doc.path = path
        fr.files.folder_override = None
        si = doc.attrs.pop("source_info", None)  # read with the data (lcms_sources)
        if si is None and path.lower().endswith(".lcd"):
            try:
                si = lcms_data.read_sample_info(path)
            except Exception:
                si = {}
        if si is not None:
            doc.sample_info = si
            bits = [si.get("sample_name") or ""] + ([si["acquired"].strftime("%d.%m.%Y %H:%M")]
                                                   if si.get("acquired") else [])
            if any(bits):
                doc.file_kind = "  \u00b7  ".join(b for b in bits if b)
        T._save({fr.FOLDER_KEY: os.path.dirname(path.rstrip("\\/"))})
        doc.doc_title = "%s: %s" % (fr.TITLE, name)
        if fr.active is doc:
            fr.SetTitle(doc.doc_title)
        doc.chips["file"].SetValue(name)
        try:
            doc.data_loaded(path, data, errors)
        except Exception as ex:
            _report_error(fr, "%s was read, but it could not be shown completely" % name, ex)
        import undo
        undo.file_opened(fr, doc)  # nothing to undo yet for this file
        fr.files.Refresh()
        fr.Layout()


def _read_problem(path, detail):
    """Why a data file could not be read, in plain words, with the reader's
    own message as the detail."""
    t = (detail or "").lower()
    if not os.path.exists(path):
        why = "the file is no longer there (moved, renamed or deleted?)"
    elif "permission" in t or "denied" in t or "sharing" in t or "cannot open" in t or "being used" in t:
        why = "the file is in use by another program or cannot be read"
    elif "tlm raw data" in t:
        why = "this is a triple quadrupole file; its MS data cannot be read yet"
    elif "storage found" in t or "no ms scans" in t:
        why = "the file holds no MS scans"
    elif any(k in t for k in ("ole2", "cfb", "sector", "fat chain", "compound file", "container", "truncat")):
        why = "the file is incomplete or damaged (still being recorded or copied?)"
    else:
        return detail or "unknown error"
    return "%s (%s)" % (why, detail) if detail else why


class LCMSFrame(PostrunFrame):
    TITLE = TITLE
    HELP = HELP_TEXT
    START_HINT = "Open a data file (Ctrl+O) or drop it here"
    FOLDER_KEY = "lcms_folder"
    _linking = False

    def link_time(self, src, t=None, rng=None):
        """A time (click, slider) or a time range (drag) chosen in the MS or
        the PDA view is chosen in the other one too, shifted by the MS
        detector delay (PDA settings): MS time = PDA time + delay."""
        if self._linking or not T._load().get("link_ms_pda", True):
            return
        ms, pda = getattr(self, "ms", None), getattr(self, "pda", None)
        if ms is None or pda is None or ms.data is None or pda.pda is None:
            return
        delay = _num(pda.delay, 0.0) or 0.0
        self._linking = True
        try:
            if src == "pda":
                if rng:
                    ms.avg = (rng[0] + delay, rng[1] + delay)
                    ms.pick_t = None
                    ms.fill_range_fields()
                    ms.plot_chroms(keep_view=True)
                    ms.compute_average()
                elif t is not None:
                    ms.chrom_pick(t + delay, None, None)
            else:
                if rng:
                    pda.avg = (rng[0] - delay, rng[1] - delay)
                    pda.fill_range_fields()
                    pda.plot_chroms(keep_view=True)
                    pda.compute_uv()
                elif t is not None and pda.wl_sel is not None:
                    pda._clear_range()  # one time in the MS view: one time in the PDA view
                    pda.on_cross(t - delay, pda.wl_sel, True)
        except Exception as ex:
            _log("linking MS and PDA: %s" % ex)
        finally:
            self._linking = False

    def link_menu(self):
        on = bool(T._load().get("link_ms_pda", True))
        return [(None, None), ("Link MS and PDA times",
                               lambda: T._save({"link_ms_pda": not on}), on)]

    def build_tabs(self, book):
        self.pda_data = None
        self.ms = MSTab(book, self)
        self.pda = PDATab(book, self)
        import lcms_compare  # the Compare view is shared by the window: this page only stands for it
        return [("Mass spectrometry", self.ms), ("PDA (UV/Vis)", self.pda),
                ("Compare", lcms_compare.ComparePage(book, self))]

    # ------------------------------------------------------- Compare view (all open files)
    def compare_tab(self, create=True):
        """The Compare view of the window (made when first shown)."""
        fr = self.window()
        tab = fr.__dict__.get("_compare_tab")
        if tab is None and create:
            import lcms_compare
            tab = lcms_compare.CompareTab(fr.docbook.GetParent(), fr)
            fr.docbook.GetContainingSizer().Add(tab, 1, wx.EXPAND)
            tab.Hide()
            fr.__dict__["_compare_tab"] = tab
            tab.files_changed()
            pend = fr.__dict__.pop("_method_compare", None)  # a method applied before the view was made
            if pend:
                tab.method_set(pend)
        return tab

    def _compare_chosen(self):
        doc = self.window().active
        try:
            i = doc.book.GetSelection()
            return i >= 0 and bool(getattr(doc.pages[i][1], "is_compare", False))
        except Exception:
            return False

    def _sync_compare(self):
        """The Compare view instead of the views of the file while its page is chosen."""
        fr = self.window()
        show = fr._compare_chosen()
        tab = fr.compare_tab(create=show)
        if tab is None:
            return
        try:
            if show and not tab.IsShown():
                fr.docbook.Hide()
                # the settings of the comparison are in its side panel: shown unless hidden in this view
                side = not T._load().get("compare_side_hidden", False)
                tab.side.Show(side)
                fr._set_top_label(fr.panel_btn, "Hide panel" if side else "Show panel")
                tab.Show()
                tab.Layout()
                tab.select_doc(fr.active)
            elif not show and tab.IsShown():
                tab.Hide()
                fr._leave_compare_pages()
                fr.docbook.Show()
                fr._set_top_label(fr.panel_btn, "Hide panel" if fr._side_shown else "Show panel")
            fr.toolbar.Layout()
            fr.docbook.GetParent().Layout()
        except RuntimeError:
            pass

    def _leave_compare_pages(self):
        """The Compare view was left: the other files, put on their Compare page while they were
        chosen during the comparison, show the view chosen now (a click on one of them in the list
        of files brought the Compare view back)."""
        fr = self.window()
        act = fr.active
        try:
            i = act.book.GetSelection()
        except (RuntimeError, AttributeError):
            return
        for d in fr.docs:
            try:
                k = d.book.GetSelection()
                if d is act or k < 0 or not getattr(d.pages[k][1], "is_compare", False):
                    continue
                en = list(getattr(d, "seg_enabled", None) or [True] * len(d.pages))
                j = i if 0 <= i < len(en) and en[i] else next((n for n, ok in enumerate(en) if ok), 0)
                d.book.ChangeSelection(j)
            except (RuntimeError, IndexError, AttributeError):
                pass

    def on_tab(self, i):
        def change():
            self.book.SetSelection(i)
            self._sync_compare()
        change_at_once(self.window(), change)

    @_view_at_once
    def activate(self, doc):
        fr = self.window()
        tab = fr.__dict__.get("_compare_tab")
        if tab is not None and tab.IsShown() and doc is not None and doc in fr.docs and doc.pages:
            # another file chosen while comparing: the Compare view stays (that file is chosen in its list)
            k = next((i for i, (_, pg) in enumerate(doc.pages) if getattr(pg, "is_compare", False)), -1)
            if k >= 0 and doc.book.GetSelection() != k:
                doc.book.ChangeSelection(k)
        PostrunFrame.activate(self, doc)
        fr._sync_compare()
        if tab is not None and tab.IsShown():
            tab.select_doc(doc)

    def show_page(self, page):
        tab = self.window().__dict__.get("_compare_tab")
        if tab is not None and tab.IsShown() and self._is_active() and page is getattr(self, "pda", None) and \
                getattr(self, "ms_data", None) is None:
            # a file without MS data opened while comparing (data_loaded shows such a file in its PDA
            # view): the Compare view stays, and lists it
            return
        PostrunFrame.show_page(self, page)
        if self._is_active():
            self.window()._sync_compare()

    def close_doc(self, doc):
        PostrunFrame.close_doc(self, doc)
        tab = self.window().__dict__.get("_compare_tab")
        if tab is not None:
            tab.files_changed()

    @_view_at_once
    def set_side(self, show=None, remember=True):
        fr = self.window()
        tab = fr.__dict__.get("_compare_tab")
        if tab is not None and tab.IsShown():  # the Compare view keeps its own choice
            show = (not tab.side.IsShown()) if show is None else bool(show)
            tab.side.Show(show)
            tab.Layout()
            fr._set_top_label(fr.panel_btn, "Hide panel" if show else "Show panel")
            fr.toolbar.Layout()
            if remember:
                T._save({"compare_side_hidden": not show})
            return
        PostrunFrame.set_side(self, show, remember)

    def chip_specs(self):
        return [("file", "File"), ("ms", "MS"), ("pda", "PDA")]

    def accepts(self, path):
        return data_formats.detect(path)[1] in data_formats.LCMS_KINDS

    def browse_accepts(self, path):
        """Data sets listed in the file panel (a Thermo .raw file without
        reading its header)."""
        low = path.lower().rstrip("\\/")
        if low.endswith((".raw", ".lcd", ".mzml", ".mzxml", ".cdf")):
            return os.path.isfile(path) or data_formats.detect(path)[1] == "waters_raw"
        return low.endswith(".d") and data_formats.detect(path)[1] == "agilent_chemstation"

    def load(self, path):
        path = data_formats.data_path(path)
        if data_formats.detect(path)[1] == "sciex_wiff":
            import lcms_sources
            wx.MessageBox(lcms_sources.WIFF_TEXT, self.TITLE, wx.ICON_INFORMATION)
            return
        PostrunFrame.load(self, path)

    def raw_open_items(self):
        return [("Open a data file…", self.open_files),
                ("Open a data folder (Agilent .D, Waters .raw)…", self.open_folder)]

    def on_open(self, e=None):
        self.open_files()

    def open_files(self):
        dlg = wx.FileDialog(self.window(), "Open LC-MS or HPLC data", defaultDir=self.folder(),
                            wildcard=data_formats.LCMS_WILDCARD,
                            style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST | wx.FD_MULTIPLE)
        try:
            if show_open_dialog(dlg) == wx.ID_OK:
                for p in dlg.GetPaths():
                    self.load(p)
        finally:
            dlg.Destroy()

    def open_folder(self):
        dlg = wx.DirDialog(self.window(), "Open an Agilent .D or Waters .raw folder", defaultPath=self.folder(),
                           style=wx.DD_DEFAULT_STYLE | wx.DD_DIR_MUST_EXIST)
        try:
            if show_open_dialog(dlg) != wx.ID_OK:
                return
            path = dlg.GetPath()
        finally:
            dlg.Destroy()
        if not self.accepts(path):
            wx.MessageBox("%s is not a data folder LCMS Analysis can read." % path, self.TITLE, wx.ICON_WARNING)
            return
        self.load(path)

    def read_data(self, path, progress):
        if data_formats.detect(path)[1] != "shimadzu_lcd":  # other vendors: lcms_sources
            import lcms_sources
            data, notes = lcms_sources.read(path, progress)
            self.source_info = data.pop("info", {})
            if data.get("ms") is None and data.get("pda") is None:
                raise RuntimeError("no MS or UV data in this file")
            return data, notes
        ms = pda = None
        errors = []
        engine = None
        engine_err = ""
        try:  # the C++ core when available (same data, a few seconds instead of a minute)
            import msengine_py
            engine = msengine_py.open_file(path, progress)
        except Exception as ex:
            print("[engine] Python reader used: %s" % ex)
            engine, engine_err = None, str(ex)
        if engine is not None:
            ms = engine
            try:
                pda = engine.pda()
            except Exception as ex:
                errors.append("PDA data: %s" % _read_problem(path, str(ex)))
            return {"ms": ms, "pda": pda}, errors
        try:
            ms = lcms_data.LCDFile(path, progress=progress)
        except Exception as ex:
            # the reason in plain words; the message of the C++ reader is the clearer detail
            errors.append("MS data: %s" % _read_problem(path, engine_err or str(ex)))
        try:
            if lcms_pda.has_pda(path):
                wx.CallAfter(self.SetStatusText, "Reading %s: PDA data" % os.path.basename(path), 0)
                pda = lcms_pda.PDAData(path)
        except Exception as ex:
            errors.append("PDA data: %s" % _read_problem(path, str(ex)))
        return {"ms": ms, "pda": pda}, errors

    def data_loaded(self, path, data, errors):
        ms, pda = data.get("ms"), data.get("pda")
        self.ms_data, self.pda_data = ms, pda
        if ms:
            pols = [e["polarity"] for e in ms.events]
            txt = " / ".join({"+": "+", "-": "−"}.get(p, "?") for p in pols) if all(pols) else \
                "%d event%s" % (ms.n_events, "s" if ms.n_events > 1 else "")
            self.chips["ms"].SetValue("ESI " + txt if all(pols) else txt)
        else:
            self.chips["ms"].SetValue("none")
        self.chips["pda"].SetValue(("%.0f to %.0f nm" % (pda.wavelengths[0], pda.wavelengths[-1]) if
                                    len(pda.wavelengths) > 1 else "%.0f nm" % pda.wavelengths[0]) if pda else "none")
        self.ms.set_data(ms)
        self.pda.set_data(pda, ms)
        self.enable_segment(0, ms is not None)
        self.enable_segment(1, pda is not None)
        cur = self.book.GetSelection()
        if ms is None and pda is not None:
            self.show_page(self.pda)
        elif cur == 1 and pda is None:
            self.show_page(self.ms)
        info = []
        if ms:
            info.append(ms.summary())
            if ms.saturated.any():
                info.append("%d scans with saturated signal" % int(ms.saturated.sum()))
        if pda:
            info.append("PDA " + pda.summary())
        if errors:
            info.append("; ".join(errors))
        self.SetStatusText(" | ".join(info), 0)
        self.set_file_status(os.path.basename(path))
        tab = self.window().__dict__.get("_compare_tab")
        if tab is not None:  # the Compare view lists the new file
            tab.request_files_changed()


UniLCMSFrame = LCMSFrame  # older name


class _Drop(wx.FileDropTarget):
    def __init__(self, frame):
        wx.FileDropTarget.__init__(self)
        self.frame = frame

    def OnDropFiles(self, x, y, names):
        ok = [n for n in names if self.frame.accepts(n)]
        for n in ok:  # every file dropped is opened
            wx.CallAfter(self.frame.load, n)
        bad = [os.path.basename(n.rstrip("\\/")) for n in names if n not in ok]
        if bad:  # say why nothing (or not everything) opened
            wx.CallAfter(self.frame.SetStatusText, "Not opened (not a data file of this window): %s"
                         % ", ".join(bad[:5]) + (" and %d more" % (len(bad) - 5) if len(bad) > 5 else ""), 0)
        return bool(ok)


def open_window(path=None):
    """The LCMS Analysis window; path: a data file or a list of them (each
    opens as a file of the window)."""
    t0 = time.perf_counter()
    paths = [p for p in (path if isinstance(path, (list, tuple)) else [path]) if p]
    f = LCMSFrame(paths[0] if paths else None)
    for p in paths[1:]:
        wx.CallAfter(f.load, p)
    f.Show()
    f.Raise()
    _log("window built in %.1f s" % (time.perf_counter() - t0))
    return f


def run(path=None):
    """Stand-alone start (own wx.App); path: a data file or a list of them."""
    app = wx.App(False)
    app.SetAppName(APP)
    open_window(path)
    app.MainLoop()

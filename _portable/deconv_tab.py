"""
Deconvolution inside the mass spectrometry view of LCMS Postrun and HRMS
Postrun (MS Analysis): no separate window or tab.

Right click a spectrum > Deconvolute this spectrum (or Deconvolute in the
top bar): a window asks for the settings (DeconvDialog); OK closes it and
runs the deconvolution on the spectrum shown for that scan event (averaged
range or single scan). The fit and the charge states of the selected mass
are drawn on that spectrum, and a result row (zero charge mass spectrum,
charge states) plus a mass table appear below it. Every run adds a row
(the newest at the top, at most 8 are kept); the table and the fit on the
spectrum follow the active result, the newest or the one clicked last.
The settings are kept for the next time.
  Bayesian deconvolution  Bayesian charge state deconvolution (Marty et
                          al., Anal. Chem. 2015, 87, 4370)
  Maximum entropy         charge state deconvolution by maximum entropy
                          (ms_deconv.py, written for MS Analysis)
  IsoDec                  monoisotopic masses from isotope resolved data
                          (IsoDec; HRMS spectra only)
Results (mass spectrum, peak list, input spectrum) are saved in the output
folder of the data file.
"""
import os
import sys
import time
import pickle
import struct
import threading
import subprocess
import weakref

import numpy as np
import wx

import unidec_theme as T
from unidec_theme import C, ui_font
import unilcms as U
import ms_deconv as D

METHODS = ["Bayesian deconvolution", "Maximum entropy", "IsoDec (isotope resolved)"]
METHOD_TAGS = ["unidec", "maxent", "isodec"]
PRESETS = [
    ("Peptides and small proteins, 0.5 to 20 kDa", dict(z=(1, 30), mass=(500, 20000), step=0.5, width=0, win=5, thr=5)),
    ("Proteins, denatured, 5 to 80 kDa", dict(z=(5, 60), mass=(5000, 80000), step=1, width=0, win=20, thr=5)),
    ("Antibodies and large proteins, 20 to 200 kDa", dict(z=(10, 100), mass=(20000, 200000), step=2, width=0, win=50,
                                                          thr=5)),
    ("Native proteins and complexes, 10 kDa to 1 MDa", dict(z=(5, 100), mass=(10000, 1000000), step=10, width=0,
                                                            win=500, thr=5)),
    ("Isotope resolved, peptides (HRMS)", dict(z=(1, 10), mass=(1000, 6000), step=0.01, width=0, win=0.4, thr=2,
                                               method=1, iso=1, hr_only=True)),
    ("Monoisotopic masses (IsoDec, HRMS)", dict(z=(1, 30), mass=(100, 50000), step=0.01, width=0, win=0.5, thr=1,
                                                method=2, hr_only=True)),
    ("Custom", None),
]
DEFAULT_PRESET = 1


def suggest_ranges(spec, hr, sign=1, carrier=None):
    """Charge and mass ranges from the spectrum: the tallest peak is given
    the charge whose neighbouring charge states (and, for resolved data,
    isotope spacing) are best supported. Returns a dict or None."""
    spec = np.asarray(spec, float)
    spec = spec[spec[:, 1] > 0]
    if len(spec) < 20:
        return None
    x, y = spec[:, 0], spec[:, 1]
    i = int(np.argmax(y))
    mz, top = float(x[i]), float(y[i])
    H = float(carrier) if carrier else D.PROTON * (-1 if sign < 0 else 1)

    def height(m, tol):
        lo, hi = np.searchsorted(x, [m - tol, m + tol])
        return float(y[lo:hi].max()) / top if hi > lo else 0.0

    tol_of = (lambda m: max(m * 15e-6, 0.004)) if hr else (lambda m: 0.35)
    best, best_z = -1.0, None
    for z in range(1, 101):
        M = z * (mz - H)
        if M <= 0:
            continue
        score = 0.0
        for k in (-3, -2, -1, 1, 2, 3):
            zz = z + k
            if zz < 1:
                continue
            m2 = (M + zz * H) / zz
            if x[0] <= m2 <= x[-1]:
                score += height(m2, tol_of(m2)) * (1.0 if abs(k) == 1 else 0.5)
        if hr:
            d = 1.00235 / z
            iso = height(mz + d, tol_of(mz)) + height(mz - d, tol_of(mz))
            half = height(mz + d / 2, tol_of(mz)) + height(mz - d / 2, tol_of(mz))  # z too small by 2
            score += 2.0 * iso - 1.0 * half
        if score > best + 1e-9:
            best, best_z = score, z
    if best_z is None:
        return None
    z = best_z
    M = z * (mz - H)
    res_power = D.estimate_resolution(spec) if hr else 0
    resolved = bool(hr and res_power > 1.3 * M)  # isotopes of the species resolved

    def present(zz):
        m = (M + zz * H) / zz
        if not (x[0] <= m <= x[-1]) or height(m, tol_of(m)) < 0.05:
            return False
        if resolved:  # its isotope peak must be there too, at 1/z spacing
            d = 1.00235 / zz
            return height(m + d, tol_of(m)) > 0.02 or height(m - d, tol_of(m)) > 0.02
        return True
    # charge states present (only these are used: charges without peaks let
    # the deconvolution explain other compounds or noise with this mass)
    zs = [zz for zz in range(max(1, z - 40), z + 41) if present(zz)]
    zlo, zhi = (min(zs), max(zs)) if zs else (z, z)
    out = {"mass_guess": M, "resolved": resolved, "resolution": res_power}
    if resolved:
        out["z"] = (zlo, zhi)
        out["mass"] = (max(1.0, np.floor(M - 30)), np.ceil(M + 80))
        out["step"] = 0.01 if M < 20000 else 0.02
    else:
        out["z"] = (max(1, zlo - 1), zhi + 1)
        span = max(0.1 * M, 200.0)
        out["mass"] = (max(1.0, np.floor(M - span)), np.ceil(M + span))
        out["step"] = 0.5 if M < 20000 else (1.0 if M < 100000 else 5.0)
    # m/z range: the charge states used, with a margin
    lo_m, hi_m = out["mass"]
    a = (min((lo_m + zz * H) / zz for zz in range(out["z"][0], out["z"][1] + 1)))
    b = (max((hi_m + zz * H) / zz for zz in range(out["z"][0], out["z"][1] + 1)))
    pad = 0.01 * (b - a) + 2.0
    out["mz"] = (max(float(x[0]), np.floor(a - pad)), min(float(x[-1]), np.ceil(b + pad)))
    return out


class TableCard(wx.Panel):
    """White card with a title, an information line and a list."""

    def __init__(self, parent, title, cols, on_select=None, on_activate=None, empty="", multi=False):
        wx.Panel.__init__(self, parent, style=wx.BORDER_NONE)
        self.SetBackgroundStyle(wx.BG_STYLE_PAINT)
        self.SetBackgroundColour(C["bg"])
        self.title, self.info, self.empty = title, "", empty
        self.badge = None  # (text, colour) of a warning pill after the title (e.g. a poor result)
        self.on_select, self.on_activate = on_select, on_activate
        self.extra_menu = None  # callback() -> menu entries shown before Copy table
        self.on_delete = None  # callback(list of rows) for the Delete key
        self.list = wx.ListCtrl(self, style=wx.LC_REPORT | wx.BORDER_NONE | (0 if multi else wx.LC_SINGLE_SEL))
        self.list.SetFont(ui_font(9, 400))
        self.set_columns(cols)
        s = wx.BoxSizer(wx.VERTICAL)
        s.AddSpacer(self.FromDIP(30))
        s.Add(self.list, 1, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, self.FromDIP(10))
        self.SetSizer(s)
        self.Bind(wx.EVT_PAINT, self._paint)
        self.Bind(wx.EVT_SIZE, lambda e: (self.Layout(), self.Refresh(), e.Skip()))
        self.list.Bind(wx.EVT_LIST_ITEM_SELECTED, lambda e: self.on_select and self.on_select(e.GetIndex()))
        self.list.Bind(wx.EVT_LIST_ITEM_ACTIVATED, lambda e: self.on_activate and self.on_activate(e.GetIndex()))
        self.list.Bind(wx.EVT_CONTEXT_MENU, self._menu)
        self.list.Bind(wx.EVT_KEY_DOWN, self._key)
        self.rows = []

    def selected(self):
        out, i = [], self.list.GetFirstSelected()
        while i >= 0:
            out.append(i)
            i = self.list.GetNextSelected(i)
        return out

    def _key(self, e):
        if e.GetKeyCode() in (wx.WXK_DELETE, wx.WXK_BACK) and self.on_delete and self.selected():
            self.on_delete(self.selected())
        elif e.GetKeyCode() == ord("A") and e.ControlDown():
            for i in range(self.list.GetItemCount()):
                self.list.Select(i)
        else:
            e.Skip()

    def set_columns(self, cols):
        self.list.ClearAll()
        self.cols = cols
        for n, (name, w) in enumerate(cols):
            self.list.InsertColumn(n, name, wx.LIST_FORMAT_RIGHT if n else wx.LIST_FORMAT_LEFT, width=self.FromDIP(w))

    def set_rows(self, rows, info=""):
        self.rows = rows
        self.list.DeleteAllItems()
        for n, r in enumerate(rows):
            i = self.list.InsertItem(n, str(r[0]))
            for c, v in enumerate(r[1:], 1):
                self.list.SetItem(i, c, str(v))
        self.info = info
        self.Refresh()

    def set_badge(self, text=None, colour=None, tip=""):
        """Warning pill after the title (None removes it); tip: tooltip of
        the title line."""
        self.badge = (text, colour or T.GROUP["red"]) if text else None
        tip = tip or ""
        if tip != getattr(self, "_tip", None):
            self._tip = tip
            self.SetToolTip(tip)
        self.Refresh()

    def select(self, i, only=True):
        if 0 <= i < self.list.GetItemCount():
            if only:
                for k in self.selected():
                    if k != i:
                        self.list.Select(k, False)
            self.list.Select(i)
            self.list.EnsureVisible(i)

    def text(self):
        head = "\t".join(c[0] for c in self.cols)
        return "\n".join([head] + ["\t".join(str(v) for v in r) for r in self.rows])

    def _menu(self, e):
        items = list(self.extra_menu()) if self.extra_menu else []
        if self.rows:
            items.append(("Copy table", lambda: (U._clip(self.text()), U._status(self, "Table copied"))))
        if U._stack_of(self) is not None:
            items += U.tile_menu(self)
        if not items:
            return
        m = U.build_menu(self, items)
        U.popup_menu(self, m)
        m.Destroy()

    def _paint(self, e):
        dc = wx.AutoBufferedPaintDC(self)
        dc.SetBackground(wx.Brush(wx.Colour(C["bg"])))
        dc.Clear()
        gc = T.crisp(dc)
        w, h = self.GetClientSize()
        r = self.FromDIP(12)
        gc.SetPen(wx.TRANSPARENT_PEN)
        gc.SetBrush(wx.Brush(wx.Colour(14, 28, 48, 12)))
        gc.DrawRoundedRectangle(1, 3, w - 2, h - 4, r)
        gc.SetBrush(wx.Brush(wx.Colour("#FFFFFF")))
        gc.SetPen(gc.CreatePen(wx.GraphicsPenInfo(wx.Colour(C["line"])).Width(1.0)))
        gc.DrawRoundedRectangle(0.5, 0.5, w - 1.5, h - 3, r)
        gc.SetFont(ui_font(10, 700), wx.Colour(C["text"]))
        tw, th = gc.GetTextExtent(self.title)
        gc.DrawText(self.title, self.FromDIP(14), self.FromDIP(9))
        if self.badge:
            bw = U.draw_badge(gc, self, self.badge[0], self.badge[1], self.FromDIP(14) + tw + self.FromDIP(10),
                              self.FromDIP(9), th)
            tw += self.FromDIP(10) + bw
        gc.SetFont(ui_font(9, 400), wx.Colour(C["muted"]))
        txt = self.info or self.empty
        room = w - tw - self.FromDIP(40)
        if gc.GetTextExtent(txt)[0] > room:
            while txt and gc.GetTextExtent(txt + "…")[0] > room:
                txt = txt[:-1]
            txt = txt + "…" if txt else ""
        sw, sh = gc.GetTextExtent(txt)
        gc.DrawText(txt, self.FromDIP(14) + tw + self.FromDIP(10), self.FromDIP(9) + (th - sh) / 2.0)


def _mass_fmt(m, step):
    if step < 0.05:
        return "%.4f" % m
    if step < 0.5:
        return "%.2f" % m
    return "%.1f" % m


UNIDEC, MAXENT, ISODEC = 0, 1, 2
# charge carriers: mass added per charge, with its sign (m/z = (M + z a) / z)
CARRIERS_POS = [("H\u207a (protonation)", 1.007276467), ("Na\u207a", 22.989218), ("K\u207a", 38.963158),
                ("NH\u2084\u207a", 18.033823), ("Custom mass", None)]
CARRIERS_NEG = [("H\u207a loss (deprotonation)", -1.007276467), ("Cl\u207b (chloride adduct)", 34.969401),
                ("Formate (HCOO\u207b)", 44.998203), ("Acetate (CH\u2083COO\u207b)", 59.013853),
                ("Custom mass", None)]
ADDUCTS = CARRIERS_POS  # older name


def carriers(sign):
    """Charge carriers offered for positive or negative ions."""
    return CARRIERS_NEG if (sign or 1) < 0 else CARRIERS_POS

# every setting of the window: key, section, label, kind, default, unit, methods (None = all), hr only, tip
# kinds: text, choice (extra = list), check
SETTINGS = [
    ("adduct", "spectrum", "Charge carrier", "choice", 0, None, None, False,
     "Ion carrying the charge; the list follows the polarity"),
    ("adduct_custom", "spectrum", "Carrier mass", "text", "", "Da", None, False,
     "Mass added per charge, with its sign (Na+: 22.9892)"),
    ("baseline", "spectrum", "Subtract background", "check", None, None, None, False,
     "Removes broad humps of chemical noise first"),
    ("bl_width", "spectrum", "Baseline window", "text", "15", "peak widths", None, False,
     "Anything wider than this is taken as baseline"),
    ("rp", "shape", "Resolving power", "text", "0", "0 = measured", (MAXENT,), True, "m/FWHM of the instrument"),
    ("width", "shape", "Peak width (half height)", "text", "0", "m/z, 0 = measured", (UNIDEC, MAXENT), False, ""),
    ("psfun", "shape", "Peak shape", "choice", 0, None, (UNIDEC,), False, ""),
    ("iso", "shape", "Isotopes", "choice", 0, None, (MAXENT,), True,
     "Resolved needs isotopes resolved in the data"),
    ("isotopemode", "shape", "Isotope mode", "choice", 0, None, (UNIDEC,), False, "Not used by this engine"),
    ("disp_rp", "shape", "Show at resolving power", "text", "", "empty = as fitted", (MAXENT,), True,
     "Lower: smoother peaks, heights kept"),
    ("numit", "unidec", "Iterations", "text", "100", None, (UNIDEC,), False, ""),
    ("zzsig", "unidec", "Charge smoothing", "text", "1", None, (UNIDEC,), False,
     "1 or more; larger: smoother charge envelopes"),
    ("psig", "unidec", "Point smoothing", "text", "1", None, (UNIDEC,), False, "0 = off"),
    ("beta", "unidec", "Suppression (beta)", "text", "0", None, (UNIDEC,), False,
     "Suppresses weak charge states; 0 = off"),
    ("msig", "unidec", "Mass smoothing", "text", "0", "Da", (UNIDEC,), False, "0 = off"),
    ("poolflag", "unidec", "m/z to mass", "choice", 2, None, (UNIDEC,), False, ""),
    ("binning", "unidec", "Data reduction", "choice", 0, None, (UNIDEC,), False, ""),
    ("mzbins", "unidec", "Bin size", "text", "", "m/z, empty = automatic", (UNIDEC,), False,
     "Constant resolution: m/z of the first bin"),
    ("smooth", "unidec", "Data smoothing", "text", "0", "points", (UNIDEC,), False, "Gaussian; 0 = off"),
    ("rounds", "maxent", "Rounds", "text", "10", None, (MAXENT,), False,
     "Maximum rounds; each gives a sharper result"),
    ("min_rounds", "maxent", "Minimum rounds", "text", "4", None, (MAXENT,), False,
     "Rounds done before the result may be accepted"),
    ("iterations", "maxent", "Refinement steps per round", "text", "50", None, (MAXENT,), False, ""),
    ("phaseres", "isodec", "Model", "choice", 0, None, (ISODEC,), False, "IsoDec neural network model"),
    ("iso_matchtol", "isodec", "Matching tolerance", "text", "5", "ppm", (ISODEC,), False,
     "Tolerance for matching isotope peaks"),
    ("iso_minpeaks", "isodec", "Minimum isotope peaks", "text", "3", None, (ISODEC,), False,
     "Isotope peaks needed for a match"),
    ("iso_css_thresh", "isodec", "Pattern match threshold", "text", "0.7", "0 to 1", (ISODEC,), False,
     "Minimum match of measured and calculated patterns"),
    ("iso_knockdown_rounds", "isodec", "Overlap passes", "text", "5", None, (ISODEC,), False,
     "Repeats to separate overlapping patterns"),
    ("iso_datathreshold", "isodec", "Data threshold", "text", "0.05", None, (ISODEC,), False,
     "Relative intensity below which data are ignored"),
    ("iso_window", "isodec", "Peak window", "text", "5", "points", (ISODEC,), False,
     "Window for picking the centroids"),
    ("iso_thresh", "isodec", "Peak threshold", "text", "0.0001", None, (ISODEC,), False,
     "Relative threshold for picking the centroids"),
    ("win", "peaks", "Peak window", "text", "", "Da", None, False, "Peaks closer than this are merged into one"),
    ("thr", "peaks", "Threshold", "text", "", "% of the largest", None, False, "Smaller masses are not listed"),
    ("labels", "peaks", "Labels", "text", "12", "peaks at most", None, False, "Peaks labelled with their mass"),
]
CHOICES = {
    "adduct": [a[0] for a in ADDUCTS],
    "psfun": ["Gaussian", "Lorentzian", "Split Gaussian / Lorentzian"],
    "iso": ["Envelope (average masses)", "Resolved (isotope peaks)"],
    "isotopemode": ["Off", "Monoisotopic masses", "Average masses"],
    "poolflag": ["Integrate", "Interpolate", "Smart"],
    "binning": ["Automatic", "None", "Linear bins", "Constant resolution"],
    "phaseres": ["8 (default)", "4"],
}
SECTIONS = [("spectrum", "Spectrum and data", "blue"), ("method", "Method", "yellow"), ("ranges", "Ranges", "green"),
            ("shape", "Peak shape and isotopes", "yellow"), ("unidec", "Bayesian deconvolution settings", "blue"),
            ("maxent", "Maximum entropy settings", "blue"), ("isodec", "IsoDec settings", "blue"),
            ("peaks", "Peaks in the result", "red")]
RANGE_KEYS = ("z0", "z1", "m0", "m1", "step", "mz0", "mz1")
FIELDS = RANGE_KEYS + ("min_int",) + tuple(s[0] for s in SETTINGS if s[3] == "text")
MIN_INT_UNITS = ["counts (y axis)", "% of the base peak"]


def _default_cfg(hr):
    presets = [p for p in PRESETS if hr or not (p[1] or {}).get("hr_only")]
    cfg = {"method": 0, "preset": DEFAULT_PRESET, "pol": 0, "suggest_open": False, "mz0": "", "mz1": "", "adduct_neg": 0,
           "min_int": "0", "min_int_unit": 0}
    for key, sec, lab, kind, dflt, unit, meths, hr_only, tip in SETTINGS:
        cfg[key] = (bool(hr) if key == "baseline" else dflt) if kind != "text" else dflt
    p = presets[DEFAULT_PRESET][1]
    cfg.update(z0="%d" % p["z"][0], z1="%d" % p["z"][1], m0="", m1="",
               step="%g" % p["step"], width="%g" % p["width"], win="%g" % p["win"], thr="%g" % p["thr"])
    return cfg


def params_from(cfg, hr):
    """Deconvolution parameters (ms_deconv) from the stored settings."""
    def num(key, default=None):
        value = str(cfg.get(key, "")).strip().replace(",", ".")
        if not value:
            return default
        try:
            result = float(value)
        except ValueError:
            raise ValueError("Enter a valid number for %s" % key)
        if not np.isfinite(result):
            raise ValueError("Enter a finite number for %s" % key)
        return result
    z0, z1, m0, m1 = num("z0"), num("z1"), num("m0"), num("m1")
    if None in (z0, z1, m0, m1) or z1 < z0 or m1 <= m0 or m0 <= 0 or z0 < 1:
        raise ValueError("Check the charge and mass ranges")
    if z0 != int(z0) or z1 != int(z1):
        raise ValueError("Charge limits must be whole numbers")
    step = num("step", 1.0)
    if step <= 0:
        raise ValueError("The mass step must be larger than 0")
    for key in ("width", "rp", "disp_rp", "smooth", "zzsig", "psig", "beta", "msig", "min_int"):
        if num(key, 0) < 0:
            raise ValueError("%s must be 0 or larger" % key)
    for key in ("mzbins", "bl_width", "win"):
        value = num(key)
        if value is not None and value <= 0:
            raise ValueError("%s must be larger than 0" % key)
    for key in ("numit", "rounds", "iterations", "min_rounds", "iso_minpeaks", "iso_knockdown_rounds", "iso_window"):
        value = num(key)
        if value is not None and (value < 1 or value != int(value)):
            raise ValueError("%s must be a positive whole number" % key)
    if not 0 <= num("thr", 5) <= 100:
        raise ValueError("Peak threshold must be between 0 and 100 %")
    sign = -1 if cfg.get("pol") == 1 else 1
    lst = carriers(sign)
    k = int(cfg.get("adduct_neg" if sign < 0 else "adduct", 0) or 0)
    k = max(0, min(k, len(lst) - 1))
    adduct = lst[k][1]
    if adduct is None:
        adduct = num("adduct_custom", None)
        if adduct is None or adduct == 0:
            raise ValueError("Enter the mass of the charge carrier (Custom mass): the mass added per charge, e.g. "
                             "22.9892 for Na+, -1.00728 for the loss of H+ or 34.9694 for Cl-")
    p = {"z_range": (int(z0), int(z1)), "mass_range": (m0, m1), "mass_step": step,
         "peak_width": num("width", 0) or 0, "peak_window": num("win", 10),
         "peak_thresh": (num("thr", 5) or 0) / 100.0, "sign": sign,
         "adduct_mass": float(adduct),
         "baseline": bool(cfg.get("baseline")), "baseline_width": num("bl_width", 15.0) or 15.0,
         "labels": int(num("labels", 12))}
    a, b = num("mz0"), num("mz1")
    if a is not None and b is not None and b > a:
        p["mz_range"] = (a, b)
    elif (a is None) != (b is None) or (a is not None and b <= a):
        raise ValueError("Check the m/z range (both empty = the whole spectrum)")
    mi = num("min_int", 0)
    if mi is None or mi < 0:
        raise ValueError("The minimum intensity must be 0 or larger")
    pct = int(cfg.get("min_int_unit", 0) or 0) == 1
    if pct and mi >= 100:
        raise ValueError("The minimum intensity must be below 100 % of the base peak")
    p["min_int"], p["min_int_unit"] = float(mi), ("pct" if pct else "abs")
    if hr:
        p["tof"] = True
        p["isotopes"] = "resolved" if cfg.get("iso") == 1 else "envelope"
        p["resolution"] = num("rp", 0) or 0
        p["display_resolution"] = num("disp_rp", 0) or 0
    # UniDec
    for key, cast in (("numit", int), ("zzsig", float), ("psig", float), ("beta", float), ("msig", float),
                      ("smooth", float)):
        v = num(key, None)
        if v is not None:
            p[key] = cast(v)
    p["psfun"] = int(cfg.get("psfun", 0) or 0)
    p["isotopemode"] = int(cfg.get("isotopemode", 0) or 0)
    p["poolflag"] = int(cfg.get("poolflag", 2) if cfg.get("poolflag") is not None else 2)
    p["binning"] = ["auto", "none", "linear", "resolution"][int(cfg.get("binning", 0) or 0)]
    v = num("mzbins", None)
    if v:
        p["mzbins"] = v
    # maximum entropy
    for key in ("rounds", "min_rounds", "iterations"):
        v = num(key, None)
        if v:
            p[key] = int(v)
    # IsoDec
    p["phaseres"] = 4 if int(cfg.get("phaseres", 0) or 0) == 1 else 8
    for key in ("iso_matchtol", "iso_minpeaks", "iso_css_thresh", "iso_knockdown_rounds", "iso_datathreshold",
                "iso_window", "iso_thresh"):
        v = num(key, None)
        if v is not None:
            p[key] = v
    return p


class DeconvDialog(wx.Dialog):
    """Every setting of one deconvolution in one window: spectrum and data,
    method, ranges, peak shape and isotopes, the settings of the method, and
    the peaks shown in the result. Deconvolute closes the window and starts
    the run."""

    def __init__(self, panel, inputs, select=None):
        wx.Dialog.__init__(self, wx.GetTopLevelParent(panel.tab), title="Deconvolute",
                           style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER | wx.MAXIMIZE_BOX)
        self.panel, self.hr, self.inputs = panel, panel.hr, inputs
        self._setting = False
        self.SetBackgroundColour(wx.Colour(C["panel"]))
        self.presets = [p for p in PRESETS if self.hr or not (p[1] or {}).get("hr_only")]
        cfg = dict(_default_cfg(self.hr))
        cfg.update(panel.cfg)
        self._carrier_sel = {1: int(cfg.get("adduct", 0) or 0), -1: int(cfg.get("adduct_neg", 0) or 0)}
        self._sign, self._ready = None, False
        self.scroll = wx.ScrolledWindow(self, style=wx.VSCROLL | wx.BORDER_NONE)
        self.scroll.SetBackgroundColour(wx.Colour(C["panel"]))
        self.scroll.SetScrollRate(0, self.FromDIP(14))
        left = self.left = U.Form(self.scroll, label_w=128, note_w=320)
        right = self.right = U.Form(self.scroll, label_w=138, note_w=320)
        self.ctrl, self.rows, self.headers = {}, {}, {}
        col_of = {"spectrum": left, "method": left, "shape": left, "peaks": left, "ranges": right, "unidec": right,
                  "maxent": right, "isodec": right}
        for sec, title, colour in SECTIONS:
            fm = col_of[sec]
            self.headers[sec] = fm.section(title, T.GROUP[colour])
            if sec == "spectrum":
                self.src = fm.choice([d["label"] for d in inputs], 320, tooltip="Spectrum to deconvolute")
                fm.full(self.src)
                self.pol = fm.choice(["Positive ions", "Negative ions"], 160)
                fm.row("Polarity", self.pol)
                self.mz0 = fm.text("", 72, tooltip="Both empty: the whole spectrum")
                self.mz1 = fm.text("", 72, tooltip="Both empty: the whole spectrum")
                fm.row("m/z range", self.mz0, "to", self.mz1)
                self.min_int = fm.text(str(cfg.get("min_int", "0")), 90,
                                       tooltip="Weaker data points are ignored; 0: use all")
                self.min_unit = fm.choice(MIN_INT_UNITS, 150)
                self.min_unit.SetSelection(1 if int(cfg.get("min_int_unit", 0) or 0) == 1 else 0)
                fm.row("Minimum intensity", self.min_int, self.min_unit)
                self.note_data = fm.note("")
            elif sec == "method":
                self.method = fm.choice([m for k, m in enumerate(METHODS) if self.hr or k < 2], 200)
                fm.row("Method", self.method)
            elif sec == "ranges":
                self.preset = fm.choice([p[0] for p in self.presets], 320)
                fm.full(self.preset)
                fm.buttons(U.flat(fm, "Suggest ranges from this spectrum", icon="replot", height=28,
                                  handler=lambda e: self.suggest(),
                                  tooltip="Charge, mass and m/z ranges from the spectrum"))
                self.z0, self.z1 = fm.text("", 60), fm.text("", 60)
                fm.row("Charge", self.z0, "to", self.z1)
                self.m0, self.m1 = fm.text("", 80), fm.text("", 80)
                fm.row("Mass", self.m0, "to", self.m1, "Da")
                self.step = fm.text("", 64, tooltip="Isotope peaks: 0.01 to 0.02 Da; envelopes: 0.5 to 2 Da")
                fm.row("Mass step", self.step, "Da")
            for key, s2, lab, kind, dflt, unit, meths, hr_only, tip in SETTINGS:
                if s2 != sec or (hr_only and not self.hr):
                    continue
                if kind == "check":
                    c = fm.check(lab, bool(cfg.get(key)), tooltip=tip)
                    self.ctrl[key], self.rows[key] = c, c
                    continue
                if kind == "choice":
                    c = fm.choice(CHOICES[key], 200, tooltip=tip)
                    c.SetSelection(max(0, min(int(cfg.get(key, dflt) or 0), len(CHOICES[key]) - 1)))
                else:
                    c = fm.text(str(cfg.get(key, dflt)), 72, tooltip=tip)
                self.ctrl[key] = c
                self.rows[key] = fm.row(lab, c, *([unit] if unit else []))
        for k in ("z0", "z1", "m0", "m1", "step", "mz0", "mz1", "min_int"):
            self.ctrl[k] = getattr(self, k)
        body = wx.BoxSizer(wx.HORIZONTAL)
        body.Add(left, 0, wx.EXPAND)
        sep = wx.Panel(self.scroll, size=wx.Size(1, -1))
        sep.SetBackgroundColour(wx.Colour(C["line"]))
        body.Add(sep, 0, wx.EXPAND | wx.TOP | wx.BOTTOM, self.FromDIP(14))
        body.Add(right, 0, wx.EXPAND)
        self.scroll.SetSizer(body)

        defaults = U.flat(self, "Method defaults", handler=lambda e: self.method_defaults(),
                          tooltip="Reset the settings of this method")
        cancel = U.flat(self, "Cancel", handler=lambda e: self.EndModal(wx.ID_CANCEL))
        ok = U.flat(self, "Deconvolute", "primary", icon="deconvolve", handler=lambda e: self.accept())
        vs = wx.BoxSizer(wx.VERTICAL)
        vs.Add(self.scroll, 1, wx.EXPAND)
        vs.Add(U.dialog_buttons(self, cancel, ok, left=(defaults,)), 0, wx.EXPAND)
        self.SetSizer(vs)

        self._setting = True
        try:
            for k in RANGE_KEYS:
                getattr(self, k).SetValue(str(cfg.get(k, "")))
            if 0 <= cfg.get("method", 0) < self.method.GetCount():
                self.method.SetSelection(cfg.get("method", 0))
            self.preset.SetSelection(max(0, min(cfg.get("preset", DEFAULT_PRESET), len(self.presets) - 1)))
        finally:
            self._setting = False
        i = 0
        if select is not None:
            i = next((k for k, d in enumerate(inputs) if d.get("e") == select), 0)
        if inputs:
            self.src.SetSelection(i)
        self.source_changed()

        self.preset.Bind(wx.EVT_CHOICE, lambda e: self.apply_preset())
        self.method.Bind(wx.EVT_CHOICE, lambda e: self.method_changed())
        self.src.Bind(wx.EVT_CHOICE, lambda e: self.source_changed())
        self.pol.Bind(wx.EVT_CHOICE, lambda e: self.pol_changed())
        for c in (self.mz0, self.mz1, self.min_int):
            c.Bind(wx.EVT_TEXT, lambda e: self.data_changed())
        self.min_unit.Bind(wx.EVT_CHOICE, lambda e: self.unit_changed())
        if "iso" in self.ctrl:
            self.ctrl["iso"].Bind(wx.EVT_CHOICE, lambda e: self.iso_changed())
        self.ctrl["adduct"].Bind(wx.EVT_CHOICE, lambda e: self.method_changed())
        for c in (self.z0, self.z1, self.m0, self.m1, self.step, self.ctrl["width"], self.ctrl["win"],
                  self.ctrl["thr"]):
            c.Bind(wx.EVT_TEXT, lambda e: self._custom())
        for k in FIELDS:
            if k in self.ctrl:
                self.ctrl[k].Bind(wx.EVT_TEXT_ENTER, lambda e: self.accept())
        self.Bind(wx.EVT_CHAR_HOOK, self._keys)
        self.method_changed(fit=False)
        self._ready = True
        self._fit()

    def _keys(self, e):
        if e.GetKeyCode() == wx.WXK_ESCAPE:
            self.EndModal(wx.ID_CANCEL)
        else:
            e.Skip()

    def _fit(self):
        """Size of the window from its contents. The first time it is centred
        on the main window; later (another method shows other settings) it
        stays where the user put it, and keeps the size the user gave it."""
        self.scroll.Layout()
        self.scroll.FitInside()
        best = self.scroll.GetSizer().GetMinSize()
        sb = wx.SystemSettings.GetMetric(wx.SYS_VSCROLL_X, self)
        try:
            area = wx.Display(max(0, wx.Display.GetFromWindow(self.GetParent()))).GetClientArea()
            hmax, wmax = int(area.height * 0.9), int(area.width * 0.95)
        except Exception:
            area = None
            hmax, wmax = self.FromDIP(900), self.FromDIP(1400)
        w = min(best.width + max(sb, self.FromDIP(16)) + self.FromDIP(4), wmax)
        h = min(best.height + self.FromDIP(70), hmax)
        placed = getattr(self, "_fit_size", None)
        if placed is None:
            self.SetClientSize(wx.Size(w, h))
            self.SetMinSize(wx.Size(min(w, self.FromDIP(600)), min(h, self.FromDIP(420))))
            self.CentreOnParent()
            self._fit_size = tuple(self.GetSize())
            return
        if self.IsMaximized() or tuple(self.GetSize()) != placed:
            return  # the user sized the window: leave it
        pos = self.GetPosition()
        self.SetMinSize(wx.Size(min(w, self.FromDIP(600)), min(h, self.FromDIP(420))))
        self.SetClientSize(wx.Size(w, h))
        self._fit_size = tuple(self.GetSize())
        # same top left corner as before; moved only as far as needed to stay
        # on the screen when it grew
        x, y = pos
        if area is not None:
            sw, sh = self.GetSize()
            x = max(area.x, min(x, area.x + area.width - sw))
            y = max(area.y, min(y, area.y + area.height - sh))
        self.SetPosition(wx.Point(x, y))

    def status(self, text):
        self.panel.status(text)
        self.SetTitle("Deconvolute – " + text if len(text) < 110 else "Deconvolute")

    # ----------------------------------------------------------- settings
    def _custom(self):
        if not self._setting and self.preset.GetSelection() != len(self.presets) - 1:
            self.preset.SetSelection(len(self.presets) - 1)

    def apply_preset(self):
        p = self.presets[self.preset.GetSelection()][1]
        if not p:
            return
        self._setting = True
        try:
            self.z0.SetValue("%d" % p["z"][0])
            self.z1.SetValue("%d" % p["z"][1])
            # the mass range of the preset too (its name gives it); before, the
            # mass range typed for another kind of sample was kept, e.g. 7000 to
            # 16000 Da with the charges 10 to 100 of the antibody preset
            self.m0.SetValue("%g" % p["mass"][0])
            self.m1.SetValue("%g" % p["mass"][1])
            self.step.SetValue("%g" % p["step"])
            self.ctrl["width"].SetValue("%g" % p["width"])
            self.ctrl["win"].SetValue("%g" % p["win"])
            self.ctrl["thr"].SetValue("%g" % p["thr"])
            if "method" in p and p["method"] < self.method.GetCount():
                self.method.SetSelection(p["method"])
            if "iso" in self.ctrl:
                self.ctrl["iso"].SetSelection(p.get("iso", 0))
        finally:
            self._setting = False
        self.method_changed()

    def method_changed(self, fit=True):
        k = self.method.GetSelection()
        for key, sec, lab, kind, dflt, unit, meths, hr_only, tip in SETTINGS:
            if key not in self.rows:
                continue
            show = meths is None or k in meths
            if key == "width" and self.hr and k == MAXENT:
                show = False  # the resolving power sets the width
            if key == "adduct_custom":
                show = self.ctrl["adduct"].GetSelection() == self.ctrl["adduct"].GetCount() - 1
            row = self.rows[key]
            if isinstance(row, wx.Window):
                row.Show(show)
            else:
                row.ShowItems(show)
        for sec, meth in (("unidec", UNIDEC), ("maxent", MAXENT), ("isodec", ISODEC)):
            self.headers[sec].GetContainingSizer().ShowItems(k == meth)
        self.left.Layout()
        self.right.Layout()
        if fit:
            self._fit()

    def iso_changed(self):
        if self.ctrl["iso"].GetSelection() == 1 and (U._num(self.step, 1.0) or 1.0) > 0.1:
            self.step.SetValue("0.01")
            self.status("Mass step set to 0.01 Da for resolved isotopes")

    def method_defaults(self):
        k = self.method.GetSelection()
        d = _default_cfg(self.hr)
        for key, sec, lab, kind, dflt, unit, meths, hr_only, tip in SETTINGS:
            if key not in self.ctrl or meths is None or k not in meths:
                continue
            c = self.ctrl[key]
            if kind == "text":
                c.SetValue(str(d[key]))
            elif kind == "choice":
                c.SetSelection(int(d[key]))
            else:
                c.SetValue(bool(d[key]))
        self.status("%s settings back to their defaults" % METHODS[k].split(" (")[0])

    def current_input(self):
        i = self.src.GetSelection()
        return self.inputs[i] if 0 <= i < len(self.inputs) else None

    def source_changed(self):
        d = self.current_input()
        if d:
            self.pol.SetSelection(1 if d.get("sign", 1) < 0 else 0)
            try:
                x = np.asarray(d["spec"], float)[:, 0]
                self.mz0.SetHint("%.1f" % x.min())
                self.mz1.SetHint("%.1f" % x.max())
            except Exception:
                pass
        self.pol_changed()
        self.data_changed()

    def _data_params(self):
        a, b = U._num(self.mz0, None), U._num(self.mz1, None)
        return {"mz_range": (a, b) if a is not None and b is not None and b > a else None,
                "min_int": max(0.0, U._num(self.min_int, 0) or 0),
                "min_int_unit": "pct" if self.min_unit.GetSelection() == 1 else "abs"}

    def data_changed(self):
        """The line under the minimum intensity: time and m/z range used, base
        peak, and the points and peaks above the minimum."""
        d = self.current_input()
        if d is None or not hasattr(self, "note_data"):
            return
        try:
            s = D.data_summary(np.asarray(d["spec"], float), self._data_params())
        except Exception:
            return
        if not s.get("mz"):
            txt = "No data in this m/z range"
        else:
            what = d.get("time") or d["label"].split(": ", 1)[-1]
            txt = "%s · m/z %.1f–%.1f · %d points" % (what, s["mz"][0], s["mz"][1], s["points_used"])
        self.left.set_note(self.note_data, txt)
        if self._ready:
            self.scroll.FitInside()

    def unit_changed(self):
        """Keep the same threshold when the unit changes (counts <-> % of
        the base peak of the spectrum chosen)."""
        d = self.current_input()
        v = U._num(self.min_int, 0) or 0
        if d is not None and v > 0:
            try:
                prm = self._data_params()
                s = D.data_summary(np.asarray(d["spec"], float), dict(prm, min_int=0))
                base = s["base"]
                if base > 0:
                    if prm["min_int_unit"] == "pct":  # was counts
                        self.min_int.ChangeValue("%.3g" % min(99.0, 100.0 * v / base))
                    else:
                        self.min_int.ChangeValue(D._fmt_int(v / 100.0 * base))
            except Exception:
                pass
        self.data_changed()

    def pol_changed(self):
        """Charge carriers of the polarity chosen (negative ions: loss of H+,
        Cl-, formate, acetate); each polarity keeps its own choice."""
        sign = -1 if self.pol.GetSelection() == 1 else 1
        c = self.ctrl.get("adduct")
        if c is None:
            return
        if self._sign is not None and self._sign != sign:
            self._carrier_sel[self._sign] = c.GetSelection()
        if self._sign != sign:
            c.SetItems([a[0] for a in carriers(sign)])
            c.SetSelection(max(0, min(int(self._carrier_sel.get(sign, 0) or 0), c.GetCount() - 1)))
            self._sign = sign
        if self._ready:
            self.method_changed()

    def carrier(self):
        """Mass added per charge for the carrier chosen (signed)."""
        sign = -1 if self.pol.GetSelection() == 1 else 1
        lst = carriers(sign)
        k = self.ctrl["adduct"].GetSelection()
        a = lst[k][1] if 0 <= k < len(lst) else None
        if a is None:
            a = U._num(self.ctrl["adduct_custom"], None) if "adduct_custom" in self.ctrl else None
        return a or D.PROTON * sign

    def suggest(self):
        d = self.current_input()
        if d is None:
            self.status("No spectrum")
            return
        spec = np.asarray(d["spec"], float)
        g = suggest_ranges(spec, self.hr, -1 if self.pol.GetSelection() == 1 else 1, self.carrier())
        if g is None:
            self.status("No clear charge state series found")
            return
        z0, z1 = int(g["z"][0]), int(g["z"][1])
        if z1 <= z0 and self.method.GetSelection() == UNIDEC:
            z1 = z0 + 1  # UniDec's engine needs two charge states at least
        self._setting = True
        try:
            self.z0.SetValue("%d" % z0)
            self.z1.SetValue("%d" % z1)
            self.m0.SetValue("%g" % g["mass"][0])
            self.m1.SetValue("%g" % g["mass"][1])
            self.step.SetValue("%g" % g["step"])
            if g.get("mz"):
                self.mz0.SetValue("%g" % g["mz"][0])
                self.mz1.SetValue("%g" % g["mz"][1])
            self.ctrl["win"].SetValue("0.4" if g["resolved"] else "%g" % max(5, round(g["mass_guess"] / 1000)))
            self.preset.SetSelection(len(self.presets) - 1)
            if "iso" in self.ctrl:
                self.ctrl["iso"].SetSelection(1 if g["resolved"] else 0)
        finally:
            self._setting = False
        self.method_changed()
        msg = "Main species about %.1f Da, charges %d to %d, m/z %g to %g" % (
            g["mass_guess"], g["z"][0], g["z"][1], g["mz"][0], g["mz"][1])
        if self.hr:
            msg += ("; isotopes resolved (resolving power %.0f): mass step %g Da" % (g["resolution"], g["step"])
                    if g["resolved"] else "; isotopes not resolved at resolving power %.0f: average masses"
                    % g["resolution"])
        self.status(msg)

    def cfg(self):
        c = {k: self.ctrl[k].GetValue().strip() for k in FIELDS if k in self.ctrl}
        for key, sec, lab, kind, dflt, unit, meths, hr_only, tip in SETTINGS:
            if key not in self.ctrl:
                continue
            if kind == "choice":
                c[key] = self.ctrl[key].GetSelection()
            elif kind == "check":
                c[key] = self.ctrl[key].GetValue()
        c.update(method=self.method.GetSelection(), preset=self.preset.GetSelection(), pol=self.pol.GetSelection(),
                 suggest_open=False, min_int_unit=self.min_unit.GetSelection())
        if "adduct" in self.ctrl and self._sign is not None:
            self._carrier_sel[self._sign] = self.ctrl["adduct"].GetSelection()
        c["adduct"] = self._carrier_sel.get(1, 0)
        c["adduct_neg"] = self._carrier_sel.get(-1, 0)
        return c

    def accept(self):
        if self.current_input() is None:
            wx.MessageBox("No spectrum to deconvolute: average a time range on a chromatogram first.",
                          "Deconvolute", wx.ICON_INFORMATION, self)
            return
        c = self.cfg()
        try:
            p = params_from(c, self.hr)
            D.check_grid(p, METHOD_TAGS[c["method"]])
        except ValueError as ex:
            from ms_brand import display
            wx.MessageBox(display(str(ex)), "Deconvolute", wx.ICON_WARNING, self)
            return
        self.result_cfg, self.result_params, self.result_input = c, p, self.current_input()
        self.EndModal(wx.ID_OK)


# ==========================================================================
# the worker process (deconv_worker.py) and the progress window
# ==========================================================================
def _read_exact(stream, n):
    buf = b""
    while len(buf) < n:
        chunk = stream.read(n - len(buf))
        if not chunk:
            return None
        buf += chunk
    return buf


def _label_room(ax, screen=False):
    """Width and height of the plot in points of on-screen text, the smaller
    of the tile and the copied image (whose text is scaled); screen=True: of
    the tile only."""
    try:
        bb = ax.get_window_extent()
        dpi = ax.figure.dpi
        sw, sh = bb.width / dpi * 72.0, bb.height / dpi * 72.0
    except Exception:
        sw, sh = 600.0, 300.0
    try:
        cfg = {} if screen else U.export_settings()
        if cfg and cfg.get("mode") != "screen":
            k = max(0.3, float(cfg["text"]) / 10.0)
            sw = min(sw, (float(cfg["width_cm"]) / 2.54 * 72.0 - 50.0) / k)
            sh = min(sh, (float(cfg["height_cm"]) / 2.54 * 72.0 - 45.0) / k)
    except Exception:
        pass
    return max(sw, 60.0), max(sh, 40.0)


def _labels_inside(ax):
    """Peak labels near the left or right edge are aligned towards the plot,
    so they are not cut off on screen or in a saved image."""
    from matplotlib.text import Annotation
    x0, x1 = ax.get_xlim()
    w = (x1 - x0) or 1.0
    for t in ax.texts:
        if isinstance(t, Annotation):
            f = (t.xy[0] - x0) / w
            t.set_ha("right" if f > 0.94 else "left" if f < 0.06 else "center")


def _savetxt(path, a, fmt, header=None):
    """The same file as np.savetxt(path, a, fmt=fmt, header=header), several times faster: rows
    are formatted 65536 at a time by one string operation instead of one write per row."""
    a = np.asarray(a, float)
    if a.ndim == 1:
        a = a[:, None]
    row = (fmt if "\t" in fmt or a.shape[1] == 1 else " ".join([fmt] * a.shape[1])) + "\n"
    with open(path, "w", encoding="latin-1") as fh:  # np.savetxt's encoding (a fast codec; cp1252 was slow)
        if header:
            fh.write("# " + header.replace("\n", "\n# ") + "\n")
        for i in range(0, len(a), 65536):
            b = a[i:i + 65536]
            fh.write((row * len(b)) % tuple(b.ravel()))


RESULT_FILES = ("_input.txt", "_mass.txt", "_peaks.csv")  # the files of a result next to the data


def _result_stem(d):
    """The name of an input's result files: its spectrum name, for a single scan with the time of
    the scan (every single scan of a file had the name <file>_<pol>_scan, so a later run on another
    scan wrote over the files of an earlier result)."""
    import re
    name = str(d.get("name") or "spectrum")
    if name.endswith("_scan"):
        m = re.search(r"at (\d+(?:\.\d+)?) min", "%s %s" % (d.get("time") or "", d.get("label") or ""))
        if m:
            name += "_%smin" % m.group(1)
    return name


def _unique_base(out, stem, tag):
    """<out>/<stem>_<tag>, or with _2, _3, ... when files of that name are there already (an earlier
    result keeps its files)."""
    base = os.path.join(out, "%s_%s" % (stem, tag))
    k = 1
    while any(os.path.exists(base + ("" if k == 1 else "_%d" % k) + s) for s in RESULT_FILES):
        k += 1
    return base if k == 1 else "%s_%d" % (base, k)


def _spec_key(e, spec):
    """Key of a spectrum of scan event e: equal for the same data (the fit of a result is drawn only on
    the spectrum it was made from)."""
    try:
        a = np.asarray(spec, float)
        if a.ndim != 2 or not len(a):
            return None
        x, y = a[:, 0], a[:, 1]
        return "%d|%d|%.17g|%.17g|%.17g|%.17g" % (int(e or 0), len(a), x[0], x[-1], float(np.sum(y)),
                                                  float(np.dot(x, y)))
    except Exception:
        return None


class WorkerClient(object):
    """Runs deconvolutions in a separate Python process (started once, kept
    for the next runs), so the window never freezes and a run can be
    cancelled. Startup failures are reported without running native code
    in the window's process."""

    clients = weakref.WeakSet()

    def __init__(self):
        self.proc = None
        self.ready = False
        self.job_id = 0
        self.handlers = None  # (job_id, on_progress, on_done, on_error)
        self.lock = threading.Lock()
        self._job = None
        self.closed = False
        self.clients.add(self)

    def _python(self):
        d = os.path.dirname(sys.executable)
        for name in ("python.exe", "python"):
            cand = os.path.join(d, name)
            if os.path.isfile(cand):
                return cand
        return sys.executable

    def _log(self):
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        for d in (os.path.join(root, "logs"), os.path.join(os.path.expanduser("~"), "UniDecPortable", "logs")):
            try:
                os.makedirs(d, exist_ok=True)
                path = os.path.join(d, "deconvolution_worker.log")
                try:  # started afresh when it has grown (UniDec prints its settings on every run)
                    if os.path.getsize(path) > 2 * 1024 * 1024:
                        os.replace(path, path[:-4] + "_old.log")
                except OSError:
                    pass
                return open(path, "ab")
            except OSError:
                continue
        return subprocess.DEVNULL

    def alive(self):
        return self.proc is not None and self.proc.poll() is None

    def start(self):
        if self.closed:
            raise RuntimeError("This processing view has been closed")
        if self.alive():
            return
        here = os.path.dirname(os.path.abspath(__file__))
        flags = 0x08000000 if os.name == "nt" else 0  # CREATE_NO_WINDOW
        log = self._log()
        try:
            self.proc = subprocess.Popen([self._python(), "-u", os.path.join(here, "deconv_worker.py")],
                                         stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=log,
                                         cwd=here, creationflags=flags)
        finally:
            if hasattr(log, "close"):
                log.close()
        self.ready = False
        threading.Thread(target=self._reader, args=(self.proc,), daemon=True).start()

    def _post(self, token, callback, *args):
        def deliver():
            with self.lock:
                valid = token == self.job_id and not self.closed
            if valid:
                callback(*args)
        wx.CallAfter(deliver)

    def _reader(self, proc):
        while True:
            hdr = _read_exact(proc.stdout, 8)
            if hdr is None:
                break
            body = _read_exact(proc.stdout, struct.unpack("<Q", hdr)[0])
            if body is None:
                break
            try:
                msg = pickle.loads(body)
            except Exception:
                continue
            with self.lock:
                if proc is not self.proc:
                    return
                h = self.handlers
                if msg[0] == "ready":
                    self.ready = True
                    continue
                if msg[0] in ("done", "error"):
                    self.handlers = self._job = None
            if h is None:
                continue
            if msg[0] == "progress":
                self._post(h[0], h[1], msg[1], msg[2], msg[3])
            elif msg[0] == "done":
                self._post(h[0], h[2], msg[1])
            elif msg[0] == "error":
                self._post(h[0], h[3], msg[1], msg[2])
        # the process ended: an error, unless it was cancelled
        with self.lock:
            if proc is not self.proc:
                return
            h = self.handlers
            self.handlers = self._job = None
        if h is not None and proc is self.proc:
            self._post(h[0], h[3], "the deconvolution process stopped unexpectedly (see logs\\"
                               "deconvolution_worker.log)", "")

    def run(self, job, on_progress, on_done, on_error):
        with self.lock:
            self.job_id += 1
            token = self.job_id
            self._job = job
            self.handlers = (token, on_progress, on_done, on_error)
        try:
            self.start()
        except Exception as ex:
            with self.lock:
                self.handlers = self._job = None
            self._post(token, on_error, "Could not start the deconvolution worker: %s" % ex, "")
            return False
        data = pickle.dumps(("run", job), protocol=4)
        proc = self.proc

        def send():  # a thread: the pipe may block until the worker has loaded its libraries
            try:
                proc.stdin.write(struct.pack("<Q", len(data)))
                proc.stdin.write(data)
                proc.stdin.flush()
            except Exception as ex:
                with self.lock:
                    if proc is not self.proc or token != self.job_id:
                        return
                    h, self.handlers = self.handlers, None
                if h is not None:
                    self._post(h[0], h[3], "could not reach the deconvolution process: %s" % ex, "")
        threading.Thread(target=send, daemon=True).start()
        return True

    def cancel(self):
        with self.lock:
            self.job_id += 1
            self.handlers = self._job = None
            self.ready = False
            proc, self.proc = self.proc, None
        if proc is not None:
            try:
                if os.name == "nt" and hasattr(proc, "pid") and proc.poll() is None:
                    subprocess.run(["taskkill", "/PID", str(proc.pid), "/T", "/F"],
                                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                                   creationflags=0x08000000, timeout=5)
                proc.kill()
            except Exception:
                try:
                    proc.kill()
                except Exception:
                    pass
            for stream in (proc.stdin, proc.stdout):  # the pipe handles
                try:
                    if stream is not None:
                        stream.close()
                except Exception:
                    pass
            try:
                proc.wait(timeout=5)
            except Exception:
                pass

    def stop(self):
        self.closed = True
        self.cancel()


try:
    import atexit
    atexit.register(lambda: [w.stop() for w in list(WorkerClient.clients)])
except Exception:
    pass


def prestart_worker(worker, delay_ms=5000, force=False):
    """Starts the deconvolution process in the background (it loads its
    libraries on another processor core meanwhile), so the first
    deconvolution does not wait for them. "prestart_deconv_worker": false in
    the settings file turns this off."""
    if worker.closed or T._load().get("prestart_deconv_worker", True) is False:
        return

    def go():
        try:
            if not worker.closed and not worker.alive():
                worker.start()
        except Exception as ex:
            print("Deconvolution process not started in advance:", ex)
    wx.CallLater(delay_ms, go)


class ProgressWindow(wx.Frame):
    """Small window while a deconvolution runs: what it does, the time, and
    Cancel (the calculation stops at once)."""

    def __init__(self, parent, title, on_cancel):
        wx.Frame.__init__(self, parent, title=title, style=wx.CAPTION | wx.FRAME_FLOAT_ON_PARENT |
                          wx.FRAME_TOOL_WINDOW | wx.FRAME_NO_TASKBAR)
        self.SetBackgroundColour(wx.Colour(C["panel"]))
        self.on_cancel = on_cancel
        self.t0 = time.time()
        self.text = wx.StaticText(self, label="Starting \u2026")
        self.text.SetFont(ui_font(9, 500))
        self.gauge = wx.Gauge(self, range=100, size=wx.Size(self.FromDIP(360), self.FromDIP(8)))
        self.clock = wx.StaticText(self, label="0 s")
        self.clock.SetForegroundColour(wx.Colour(C["muted"]))
        btn = U.flat(self, "Cancel", handler=lambda e: self.cancel(), tooltip="Stop the deconvolution now")
        vs = wx.BoxSizer(wx.VERTICAL)
        pad = self.FromDIP(14)
        vs.Add(self.text, 0, wx.LEFT | wx.RIGHT | wx.TOP, pad)
        vs.Add(self.gauge, 0, wx.EXPAND | wx.ALL, pad)
        hs = wx.BoxSizer(wx.HORIZONTAL)
        hs.Add(self.clock, 1, wx.ALIGN_CENTER_VERTICAL)
        hs.Add(btn, 0)
        vs.Add(hs, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, pad)
        self.SetSizerAndFit(vs)
        self.timer = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, self._tick, self.timer)
        self.timer.Start(250)
        self.frac = None
        self._parent = parent
        self._place()
        self.Show()
        wx.CallAfter(self._place)  # again with the final size

    def _place(self):
        """Centred on the window it belongs to, inside the screen."""
        try:
            parent = self._parent
            r = parent.GetScreenRect()
            area = wx.Display(max(0, wx.Display.GetFromWindow(parent))).GetClientArea()
            w, h = self.GetSize()
            x = r.x + (r.width - w) // 2
            y = r.y + (r.height - h) // 2
            x = max(area.x, min(x, area.x + area.width - w))
            y = max(area.y, min(y, area.y + area.height - h))
            self.SetPosition(wx.Point(x, y))
        except Exception:
            try:
                self.CentreOnParent()
            except Exception:
                pass

    def _tick(self, e):
        self.clock.SetLabel("%.0f s" % (time.time() - self.t0))
        if self.frac is None:
            self.gauge.Pulse()

    def update(self, i, n, text):
        try:
            self.text.SetLabel(text[:90])
            if n > 0:
                self.frac = min(1.0, float(i) / n)
                self.gauge.SetValue(int(100 * self.frac))
            self.Layout()
        except RuntimeError:
            pass

    def cancel(self):
        if self.on_cancel:
            self.on_cancel()

    def finish(self):
        try:
            self.timer.Stop()
            self.Destroy()
        except RuntimeError:
            pass


MAX_RESULTS = 8  # result rows kept; the oldest is dropped when a ninth comes
SHOW_Z_KEY = "show_charge_states"  # setting: the charge state tile of the result rows (off by default)


def show_charge_states():
    """True when the result rows show the charge state tile next to the
    zero charge mass spectrum (right click a result > Show charge states)."""
    return bool(U.settings().get(SHOW_Z_KEY, False))


def _grouped(res):
    """Isotope resolved result with one entry per species (average mass,
    most abundant isotope, number of isotope peaks). Maximum entropy gives
    every mass an "apex" (its grid point), so that key alone does not tell:
    envelope results were shown as resolved ones (columns "Average mass" and
    "Most abundant" for envelopes)."""
    return any(q.get("n_iso") not in (None, "") for q in (res or {}).get("peaks") or [] if not q.get("added"))


def _quality(res):
    """The quality of a result, checked: {"score": 0..1 or None, "level":
    "good", "fair", "poor" or "", "text": ...}, or None when the result has
    none (maximum entropy, IsoDec, older results) or it is not readable."""
    q = (res or {}).get("quality")
    if not isinstance(q, dict):
        return None
    level = str(q.get("level") or "").strip().lower()
    if level not in ("good", "fair", "poor"):
        level = ""
    try:
        score = float(q.get("score"))
        score = min(1.0, max(0.0, score)) if np.isfinite(score) else None
    except (TypeError, ValueError):
        score = None
    from ms_brand import display
    text = display(str(q.get("text") or "").strip())
    if not level and score is None and not text:
        return None
    return {"score": score, "level": level, "text": text}


def _quality_text(qa):
    """One line for tooltips: level, score and the reason."""
    if not qa:
        return ""
    head = "Result quality: %s" % (qa["level"] or "not rated")
    if qa["score"] is not None:
        head += " (score %.2f of 1)" % qa["score"]
    return head + (". " + qa["text"] if qa["text"] else "")


def _peak_score(q):
    """Score (0 to 1) of one mass, or None (none given, added by hand, not a
    number)."""
    if q.get("added"):
        return None
    try:
        s = float(q.get("score"))
    except (TypeError, ValueError):
        return None
    return s if np.isfinite(s) else None


def _has_scores(res):
    return any(_peak_score(q) is not None for q in (res or {}).get("peaks") or [])


def _drawn_top(mass, x, lim=None):
    """Top of the peak at x as the mass spectrum is drawn (straight lines
    between its points): from the point nearest to x uphill to the highest
    point of that peak, at most lim Da from x when given. Markers sit there.
    The apex of a peak (the centroid of its top) lies between two points;
    on a coarse mass axis it is on a flank of the drawn line (1000 to 50000
    Da at 0.01 Da gives a point every 0.12 Da, below 2 points per isotope
    peak: a marker at the apex was drawn at 83 % of the peak height)."""
    xs, ys = mass[:, 0], mass[:, 1]
    n = len(xs)
    if n == 0:
        return float(x), 0.0
    i = int(min(max(np.searchsorted(xs, x), 0), n - 1))
    if i > 0 and abs(xs[i - 1] - x) < abs(xs[i] - x):
        i -= 1
    lo, hi = 0, n - 1
    if lim:
        lo = min(i, int(np.searchsorted(xs, x - lim)))
        hi = max(i, int(np.searchsorted(xs, x + lim, side="right")) - 1)
    for _ in range(20000):
        j = i
        if i > lo and ys[i - 1] > ys[j]:
            j = i - 1
        if i < hi and ys[i + 1] > ys[j]:
            j = i + 1
        if j == i:
            break
        i = j
    return float(xs[i]), float(ys[i])


def _nearest_iso(isop, x, tol=0.3):
    """Index of the labelled isotope peak nearest to x (within tol Da), or None."""
    if not isop:
        return None
    k = min(range(len(isop)), key=lambda j: abs(isop[j][0] - x))
    return k if abs(isop[k][0] - x) <= tol else None


def _charges_of(mass, res):
    """(z, m/z) of the charge states of a mass inside the m/z range of the
    data."""
    if not res or not res.get("params"):
        return []
    p = res["params"]
    a = D.carrier_mass(p)
    data = res.get("data")
    lo, hi = (float(data[0, 0]), float(data[-1, 0])) if data is not None and len(data) else (0, 1e9)
    out = []
    for z in range(p["z_range"][0], p["z_range"][1] + 1):
        mz = (mass + z * a) / z
        if lo <= mz <= hi:
            out.append((z, mz))
    return out


def _z_text_of(mass, res):
    zs = [z for z, _ in _charges_of(mass, res)]
    if not zs:
        return ""
    return "%d to %d" % (min(zs), max(zs)) if len(zs) > 2 else ", ".join(str(z) for z in zs)


def _display_apex(q):
    """A manually clicked isotope, otherwise the deconvoluted spectrum's apex.

    Raw charge-state refinement can prefer a neighbouring isotope. That is
    useful metadata, but must not silently move a marker in the mass plot.
    """
    return float(q.get("selected_apex", q.get("apex_dec", q.get("apex", q["mass"]))))


def _isotope_selection(res, x, tolerance=0.3):
    """Resolve an isotope click to its species without snapping to its mode."""
    if not _grouped(res) or x is None:
        return None
    iso = res.get("isotope_peaks") or []
    k = _nearest_iso(iso, x, tolerance)
    if k is None:
        return None
    m = float(iso[k][0])
    candidates = [(abs(float(q["mass"]) - m), i) for i, q in enumerate(res["peaks"])
                  if abs(float(q["mass"]) - m) <= max(0.3, float(q.get("n_iso") or 0) * 1.00235)]
    if not candidates:
        return None
    _, i = min(candidates)
    return i, m


def _apex_mass(res, q):
    """Mass of the selected isotope (default: deconvoluted modal isotope),
    measured at the centroid of the top half of that individual peak."""
    ap = _display_apex(q)
    k = _nearest_iso(res.get("isotope_peaks"), ap, 0.3)
    return float(res["isotope_peaks"][k][0]) if k is not None else float(ap)


def _apex2_suffix(q, fmt=" (%.4f: %.1f %%)"):
    """Kept for saved-result compatibility; only one measured peak is displayed.

    A nearly equal neighbouring isotope belongs to the same envelope. Showing
    both numbers on one marker made one peak look like two reported results.
    """
    return ""


def _view_maxima(mass, x0, x1, step):
    """Every peak of the mass spectrum inside x0..x1 (for "Label all peaks in
    view"): its apex (centroid of the top half, as the isotope labels) and
    height, the peaks at least 1 % of the tallest one in view."""
    if mass is None or len(mass) < 5:
        return []
    m = np.asarray(mass, float)
    a, b = np.searchsorted(m[:, 0], [x0, x1])
    a, b = max(0, a - 2), min(len(m), b + 2)
    if b - a < 5:
        return []
    out = [(float(xm), float(h)) for xm, h in D.isotope_peaks(m[a:b], threshold=0.01, max_n=400)]
    # satellites (a bump half way between two isotope peaks, below a tenth of a taller peak
    # less than 0.7 Da away: the deconvolution's charge state ghosts) are not labelled
    keep = [(xm, h) for xm, h in out
            if not any(abs(xm - x2) < 0.7 and h < 0.1 * h2 for x2, h2 in out if h2 > h)]
    return [(xm, h) for xm, h in keep if x0 <= xm <= x1]


class DeconvResult(object):
    """One deconvolution result and its row of tiles: the zero charge mass
    spectrum card, the charge state card (shown on request), the selected
    peak, the zoom of the mass spectrum and the rows of the mass table. A
    result hidden from the list of open files (visible False) keeps all of
    this but is not drawn and not the active one."""

    def __init__(self, res):
        self.res = res
        self.row = None
        self.card_mass = None
        self.card_z = None
        self.sel = None
        self.view = None
        self.rows = []
        self.visible = True
        self.z_weights = None  # widths of the two tiles (as dragged) while the charge states are hidden
        self.input = None  # the input of the run (this session only): Deconvolute again uses it


class DeconvPanel(object):
    """Deconvolution in the mass spectrometry view of `tab` (a
    unilcms.MSTab): the settings window, a short section of the side panel,
    one result row (zero charge mass spectrum; charge states on request) per
    deconvolution, the newest at the top, and one mass table. The table and
    the fit drawn on the spectrum belong to the active result: the newest
    one after a run, or the one whose plot was clicked last. A result can be
    hidden (kept, not drawn) and shown again, or closed (removed); the list
    of open files shows every result of its file as a branch, and is told of
    every change through `listeners` (callbacks taking this panel)."""

    instances = weakref.WeakSet()  # every panel (a setting such as the charge states applies to all)

    def __init__(self, tab, hr=False):
        self.tab = tab
        self.frame = tab.frame
        self.hr = hr
        self.inputs = []
        self.results = []  # DeconvResult, newest first (hidden ones included)
        self.active = None  # the result shown in the table and on the spectrum (never a hidden one)
        self.table = None
        self.busy = False
        self.listeners = []  # callback(panel) after results were added, activated, hidden, shown or closed
        self.worker = WorkerClient()
        self._tab_id = tab.GetId()
        tab.Bind(wx.EVT_WINDOW_DESTROY, self._destroyed)
        self._relabel_pending = []
        DeconvPanel.instances.add(self)
        try:  # the list of open files of the window shows the results of each file
            files = getattr(wx.GetTopLevelParent(tab), "files", None)
            if files is not None and hasattr(files, "watch_results"):
                files.watch_results(self)
        except Exception:
            pass
        self.key = "deconv_hr_v2" if hr else "deconv_lc_v2"  # settings of the new window (older ones dropped)
        self.cfg = self._saved_cfg() or _default_cfg(hr)
        sp = tab.side
        sp.section("Deconvolution", T.GROUP["yellow"])
        self.run_btn = U.flat(sp, "Deconvolute…", "primary", icon="deconvolve", handler=lambda e: self.deconvolute())
        sp.buttons(self.run_btn)
        btns = [U.flat(sp, "Close results", handler=lambda e: self.close())]
        if not hr:
            btns.append(U.flat(sp, "Full window", icon="deconvolve", handler=lambda e: self.open_full(),
                               tooltip="Open in the separate Deconvolute window"))
        sp.buttons(*btns)
        import deconv_shifts  # mass shift finder on the zero charge mass spectra (3.35)
        self.shifts = deconv_shifts.ShiftPanel(self, sp)
        self.note = sp.note("")
        prestart_worker(self.worker)

    def _saved_cfg(self):
        """The settings saved by the last deconvolution of this kind (LCMS or
        HRMS, any open file or window) over the defaults, or None."""
        saved = T._load().get(self.key)
        if not isinstance(saved, dict):
            return None
        cfg = _default_cfg(self.hr)
        for k, v in saved.items():
            if k not in cfg:
                continue
            d = cfg[k]
            try:  # the saved value must have the type of the default (a hand edited file may differ)
                if isinstance(d, bool):
                    v = bool(v)
                elif isinstance(d, int):
                    v = int(float(v))
                elif isinstance(d, str):
                    v = "" if v is None else str(v)
            except (TypeError, ValueError):
                continue
            cfg[k] = v
        return cfg

    # ------------------------------------------------- the active result
    # Older code (reports, menus) reads result, sel, card_mass, card_z and
    # row of the panel: they belong to the active result.
    @property
    def result(self):
        return self.active.res if self.active is not None else None

    @property
    def sel(self):
        return self.active.sel if self.active is not None else None

    @sel.setter
    def sel(self, i):
        if self.active is not None:
            self.active.sel = i

    @property
    def card_mass(self):
        return self.active.card_mass if self.active is not None else None

    @property
    def card_z(self):
        return self.active.card_z if self.active is not None else None

    @property
    def row(self):
        return self.active.row if self.active is not None else None

    # ----------------------------------------------------------- helpers
    def _destroyed(self, event):
        if event.GetId() == self._tab_id:
            self.busy = False
            self.worker.stop()
            self._end_progress()
        event.Skip()

    def status(self, text):
        self.tab.status(text)

    def _set_note(self, text):
        sp = self.tab.side
        self.note.SetLabel(text)
        self.note.Wrap(sp.FromDIP(270))
        sp.Layout()
        sp.FitInside()

    def _notify(self):
        """Tell the listeners (the list of open files) that the results or
        the busy state changed."""
        for cb in list(self.listeners):
            try:
                cb(self)
            except RuntimeError:  # its window is gone
                try:
                    self.listeners.remove(cb)
                except ValueError:
                    pass
            except Exception as ex:
                print("[deconvolution] listener failed: %s" % ex)

    def _refresh_files(self):
        """The list of open files says "deconvoluting" while a run is busy
        and lists the results: drawn again when that changes."""
        if self.listeners:
            self._notify()
            return
        try:
            files = getattr(wx.GetTopLevelParent(self.tab), "files", None)
            if files is not None:
                files.Refresh()
        except RuntimeError:
            pass

    def visible_results(self):
        """The results drawn (not hidden), newest first."""
        return [r for r in self.results if r.visible]

    def _newest_visible(self):
        return next((r for r in self.results if r.visible), None)

    def reset(self):
        if self.busy:
            self.cancel()  # a new file: the running deconvolution belongs to the old one
        self.close()
        self.inputs = []

    def _ensure_table(self):
        if self.table is not None:
            return
        tab = self.tab
        self.table = TableCard(tab.bottom, "Masses", self._cols(False), on_select=self.select_peak,
                               empty="Deconvolution results")
        self.table.Bind(wx.EVT_MOUSEWHEEL, lambda e: U.forward_wheel(self.table, e))
        self.table.Hide()

    def _make_row(self, ent):
        """The row of tiles of one result (hidden until it is laid out). The
        charge state card always exists (it is drawn, so image export and
        older code reading card_z work) but is in the row only while the
        charge states are shown."""
        tab = self.tab
        ent.row = U.SplitBox(tab.rows, wx.HORIZONTAL)
        ent.row.tile_scale = getattr(tab, "RESULT_TILE_SCALE", None)  # HRMS: full size
        ent.card_mass = U.PlotCard(ent.row, "Deconvoluted mass spectrum", "", mode="zoom")
        ent.card_z = U.PlotCard(ent.row, "Charge states", "", mode="zoom")
        self.shifts.make_cards(ent)
        self._apply_z(ent)
        for card in (ent.card_mass, ent.card_z):
            card.tool_getter = lambda: tab.tool
            card.tools_supported = set()
            card.on_menu = lambda x, y, c=card, r=ent: self._card_menu(r, c, x)
        ent.card_mass.on_pick = lambda x, y, r=ent: self._card_pick(r, x, y)
        ent.card_z.on_pick = lambda x, y, r=ent: self.activate(r)
        ent.card_mass.on_view = lambda r=ent: self._mass_view_changed(r)
        ent.card_mass.relabel_fn = lambda r=ent: self._relabel_now(r)  # copied and saved images
        ent.card_mass.readout_fmt = lambda x, y, r=ent: "%s Da   %s" % (_mass_fmt(x, self._fmt_step(r.res)),
                                                                        U._g(y))
        ent.card_mass.tip_fmt = lambda x, y, r=ent: self._mass_tip(r, x)
        ent.card_mass.image_name = lambda r=ent: self._file_base(r.res) + "_mass"
        ent.card_z.image_name = lambda r=ent: self._file_base(r.res) + "_charges"
        ent.row.on_close_tile = lambda pane, r=ent: self._close_tile(r, pane)
        ent.row.Hide()

    def _close_tile(self, ent, pane):
        """Close this tile on a result: the mass spectrum hides the result
        (the list of open files shows it again), the charge states are
        hidden in every row (Show charge states brings them back). A row
        closed as a tile was not drawn but still the result in the table."""
        if pane is ent.card_mass:
            self.set_visible(ent, False)
            return True
        if pane is ent.card_z:
            self.set_charge_states(False)
            return True
        return self.shifts.close_tile(ent, pane)

    def _apply_z(self, ent, show=None):
        """The charge state tile of a result row shown or not (the setting by
        default). Hidden: the mass spectrum takes the whole width."""
        row = ent.row
        if row is None:
            return
        show = show_charge_states() if show is None else bool(show)
        extra, ew = self.shifts.panes(ent)  # the tables of the mass shifts, when they are on
        if show:
            row.closed.discard(ent.card_z)  # also when closed by hand (Close this tile)
            zw = ent.z_weights or [3.4, 1]
            row.set_panes([ent.card_mass] + extra + [ent.card_z], [zw[0]] + ew + [zw[1]])
        else:
            if ent.card_z in row.panes:  # widths as the user left them
                ent.z_weights = [row.weight_of(ent.card_mass, 3.4), row.weight_of(ent.card_z, 1.0)]
            row.set_panes([ent.card_mass] + extra, [3.4 if extra else 1.0] + ew)
            try:
                ent.card_z.Hide()
            except RuntimeError:
                pass

    def set_charge_states(self, flag):
        """Show or hide the charge state tile of every result row (of every
        open file); remembered for the next rows and sessions."""
        flag = bool(flag)
        T._save({SHOW_Z_KEY: flag})
        for pn in list(DeconvPanel.instances):
            try:
                pn._charge_states_changed(flag)
            except RuntimeError:
                pass  # its view is gone
        self.status("Charge states shown" if flag else "Charge states hidden")

    def _charge_states_changed(self, flag):
        if not self.results:
            return
        for ent in self.results:
            self._apply_z(ent, flag)
        self._layout()  # the tile heights follow the panes of the rows
        if flag:
            for ent in self.visible_results():
                if ent.card_z is not None:
                    ent.card_z.draw()

    @staticmethod
    def _drop_row(ent):
        row, ent.row, ent.card_mass, ent.card_z = ent.row, None, None, None
        if row is not None:
            try:
                row.Destroy()
            except RuntimeError:
                pass

    def _layout(self):
        """The rows of the results shown in the tile column, newest first,
        and the table (hidden results: their rows are hidden, kept)."""
        shown = [r.row for r in self.results if r.visible and r.row is not None]
        details = getattr(self.tab, "hrms_details", None)
        if details is not None:
            details.refresh_result(self.active if shown else None)
        for r in self.results:
            if not r.visible and r.row is not None:
                try:
                    r.row.Hide()
                except RuntimeError:
                    pass
        if shown:
            self._ensure_table()
            self.tab.set_extra("dec", shown, self.table, 3.4)
        else:
            self.tab.set_extra("dec")

    def show(self, flag):
        if flag:
            self._layout()
        else:
            self.tab.set_extra("dec")

    def _card_menu(self, ent, card, x):
        """Right click on a result plot: the result becomes the active one."""
        self.activate(ent)
        res = ent.res
        shown = show_charge_states()
        own = getattr(ent, "input", None)
        cur = next((d for d in self.tab.deconv_inputs() if d.get("e") == res.get("e")), None)
        same = cur is not None and res.get("spec_key") is not None and \
            _spec_key(cur.get("e", 0), cur.get("spec")) == res.get("spec_key")
        if own is not None or same or res.get("spec_key") is None:
            items = [("Deconvolute again with other settings…", lambda: self.deconvolute(res["e"], own))]
        else:  # the result's spectrum is not kept (an older session) and another one is shown now
            items = [("Deconvolute the spectrum shown now…", lambda: self.deconvolute(res["e"]))]
        items.append((None, None))
        if card is ent.card_mass:
            lab = bool(res.get("label_all"))
            items += self._mass_menu(x) + self.shifts.menu_items(ent, x) + [
                ("Automatic peak labels" if lab else "Label all peaks in view", lambda: self.set_label_all(ent, not lab)),
                ("Zoom to the masses found", self.zoom_peaks),
                ("Full view", card.reset_view)]
        else:
            items.append(("Full view", card.reset_view))
        items += [("Hide charge states" if shown else "Show charge states", lambda: self.set_charge_states(not shown)),
                  (None, None),
                  ("Hide this result", lambda: self.set_visible(ent, False)),
                  ("Close this result", lambda: self.close_result(ent))]
        if len(self.results) > 1:
            items.append(("Close all results", self.close))
        return items + self.tab.view_menu()

    def set_label_all(self, ent, on):
        """Label every peak of the current view of this result's mass spectrum
        (kept while zooming), or go back to the automatic labels."""
        ent.res["label_all"] = bool(on)
        self.activate(ent)
        self.plot_results(ent)

    def _card_pick(self, ent, x, y):
        self.activate(ent)
        self.pick_mass(x, y)

    def activate(self, ent):
        """Make ent the result shown in the table and on the spectrum (a
        hidden result is shown again first)."""
        if ent is None or ent not in self.results:
            return
        if not ent.visible:
            ent.visible = True
            self._layout()
            if ent.card_mass is not None:
                ent.card_mass.draw()
        if ent is self.active:
            return
        old, self.active = self.active, ent
        self._ensure_table()
        self._fill_table(ent)
        self._mark_rows()
        events = {ent.res["e"]}
        if old is not None:
            events.add(old.res["e"])
        self._replot_spec(events)
        self._notify()

    def focus(self, ent):
        """A click on a result in the list of open files: shown (if hidden),
        made the active one and its row scrolled into view."""
        if ent not in self.results:
            return
        self.activate(ent)
        try:
            self.frame.show_page(self.tab)  # LCMS Postrun: the MS view, if the PDA view is shown
        except Exception:
            pass
        if ent.row is not None:
            self.tab.reveal(ent.row)
        self._notify()

    def set_visible(self, ent, flag=True):
        """Hide a result (its row is not drawn and it is not the active one,
        but it is kept) or show it again (the list of open files, or Hide
        this result in the menu of its row)."""
        flag = bool(flag)
        if ent not in self.results or ent.visible == flag:
            return
        ent.visible = flag
        old = self.active
        if not flag and ent is old:
            self.active = self._newest_visible()  # the newest one shown takes over the table
        elif flag and old is None:
            self.active = ent
        self._layout()
        if self.active is not old:
            if self.active is not None:
                self._fill_table(self.active)
            self._replot_spec({r.res["e"] for r in (old, self.active) if r is not None})
        self._mark_rows()
        if flag and ent.row is not None:
            ent.card_mass.draw()
            self.tab.reveal(ent.row)
        n = len(self.results) - len(self.visible_results())
        self.status("Result shown again: %s" % self._subtitle(ent, False) if flag else
                    "Result hidden (%d hidden)" % n)
        self._notify()

    def set_all_visible(self, flag=True):
        """Show or hide every result (the list of open files)."""
        flag = bool(flag)
        todo = [r for r in self.results if r.visible != flag]
        if not todo:
            return
        old = self.active
        for r in todo:
            r.visible = flag
        if not flag:
            self.active = None
        elif old is None:
            self.active = self._newest_visible()
        self._layout()
        if self.active is not old:
            if self.active is not None:
                self._fill_table(self.active)
            self._replot_spec({r.res["e"] for r in (old, self.active) if r is not None})
        self._mark_rows()
        for r in todo:
            if flag and r.row is not None:
                r.card_mass.draw()
        self.status("Every result shown" if flag else "Every result hidden")
        self._notify()

    def close_result(self, ent):
        """Remove one result (its row and data); the newest remaining one
        shown becomes active if it was the active one."""
        if ent not in self.results:
            return
        self.results.remove(ent)
        was_active = ent is self.active
        events = {ent.res["e"]}
        if was_active:
            self.active = self._newest_visible()
            if self.active is not None:
                events.add(self.active.res["e"])
        self._layout()
        self._drop_row(ent)
        if was_active:
            if self.active is not None:
                self._fill_table(self.active)
            self._replot_spec(events)
        self._mark_rows()
        self._notify()

    hide_result = close_result  # older name (the row menu said "Hide this result" for removing it)

    def close(self):
        """Close every result (the rows, the table and the fit on the spectrum)."""
        ents, self.results, self.active = self.results, [], None
        self._layout()
        for ent in ents:
            self._drop_row(ent)
        if ents:
            for col in self.tab.cols:
                self.tab.plot_spec(col, keep_view=True)
        self._notify()

    def _fmt_step(self, res=None):
        res = self.result if res is None else res
        if res and res.get("tag") == "isodec":
            return 0.0001
        if res and _grouped(res):
            return 0.1
        if res:
            return res["params"]["mass_step"]
        try:
            return float(str(self.cfg.get("step", "1")).replace(",", "."))
        except ValueError:
            return 1.0

    # -------------------------------------------------------------- inputs
    def refresh_inputs(self, select=None, quiet=False):
        self.inputs = self.tab.deconv_inputs()
        for d in self.inputs:
            d.setdefault("key", _spec_key(d.get("e", 0), d.get("spec")))
        if not self.inputs and not quiet:
            self.status("No spectrum yet: average a time range on a chromatogram")

    def deconvolute(self, e=None, own=None):
        """Settings window for the spectrum of scan event e; OK runs it. own: the input of a result
        (Deconvolute again): offered first when the spectrum shown now is another one."""
        if self.busy:
            self.status("The deconvolution is still running")
            return
        self.refresh_inputs(quiet=True)
        inputs = self.inputs
        if own is not None and own.get("spec") is not None:
            if not any(d.get("key") == own.get("key") and own.get("key") is not None for d in inputs):
                lab = own.get("label0", own["label"])
                own = dict(own, label="%s (of this result)" % lab, label0=lab)
                inputs = [own] + list(inputs)
            e = own.get("e", e)
        if not inputs:
            self.status("No spectrum: average a time range on a chromatogram first")
            wx.MessageBox("No spectrum to deconvolute: average a time range on a chromatogram first.",
                          "Deconvolute", wx.ICON_INFORMATION)
            return
        # the last settings used, also when they were used on another open file (each file has
        # its own panel, which read them when the file was opened)
        self.cfg = self._saved_cfg() or self.cfg
        dlg = DeconvDialog(self, inputs, e)
        try:
            if dlg.ShowModal() != wx.ID_OK:
                return
            cfg, p, d = dlg.result_cfg, dlg.result_params, dlg.result_input
        finally:
            dlg.Destroy()
        self.cfg = cfg
        T._save({self.key: cfg})
        self.run(d, p, cfg["method"])

    # ----------------------------------------------------------------- run
    def _file_base(self, res=None):
        res = self.result if res is None else res
        if res and res.get("saved"):
            return os.path.basename(res["saved"])  # as its files (unique: _2, _3, ... after the first)
        tag = res.get("tag", "deconv") if res else METHOD_TAGS[self.cfg.get("method", 0)]
        name = (res or {}).get("name") or (self.inputs[0]["name"] if self.inputs else "spectrum")
        return "%s_%s" % (name, tag)

    def run(self, d, p, k):
        if self.busy:
            return
        out = self.frame.out_dir()
        name, tag = _result_stem(d), METHOD_TAGS[k]
        spec = np.array(d["spec"], float)
        if d.get("key") is None:
            d = dict(d, key=_spec_key(d.get("e", 0), spec))
        self.busy = True
        self.run_btn.Enable(False)
        self._refresh_files()
        self._job = {"d": d, "p": p, "k": k, "t0": time.time(), "out": out, "name": name, "tag": tag, "spec": spec}
        self.status("Deconvoluting %s (%s) \u2026" % (d.get("label0", d["label"]), METHODS[k]))
        self.progress = ProgressWindow(wx.GetTopLevelParent(self.tab), "Deconvolute: %s" % METHODS[k].split(" (")[0],
                                       self.cancel)
        self.worker.run({"method": tag, "spec": spec, "params": p, "folder": out, "name": name},
                   self._on_progress, self._on_done, self._on_error)

    def _on_progress(self, i, n, text):
        from ms_brand import display
        text = display(text)
        if self.busy and getattr(self, "progress", None) is not None:
            self.progress.update(i, n, text)
            self.status("Deconvolute: " + text)

    def _on_done(self, res):
        if not self.busy:
            return  # cancelled meanwhile
        from ms_brand import display
        for field in ('method','notes','engine'):
            if field in res:res[field]=display(res[field])
        j = self._job
        d = j["d"]
        res.update({"tag": j["tag"], "params": j["p"], "input": d.get("label0", d["label"]), "name": j["name"],
                    "e": d.get("e", 0), "seconds": time.time() - j["t0"], "spec_key": d.get("key")})
        self._save(res, np.asarray(d["spec"], float), j["out"], j["name"], j["tag"])
        # the input of this result (Deconvolute again uses it when another spectrum is shown by then)
        own = {k: v for k, v in d.items() if k != "spec"}
        own["spec"] = j.get("spec", d["spec"])
        self._done(res, None, own)

    def _on_error(self, msg, tb=""):
        if tb:
            print(tb)
        if not self.busy:
            return
        from ms_brand import display
        self._done(None, display(msg))

    def cancel(self):
        if not self.busy:
            return
        self.worker.cancel()
        prestart_worker(self.worker, 3000, force=True)
        self.busy = False
        self._end_progress()
        self._refresh_files()
        try:
            self.run_btn.Enable(True)
        except RuntimeError:
            pass
        self.status("Deconvolution cancelled")

    def _end_progress(self):
        pw = getattr(self, "progress", None)
        self.progress = None
        if pw is not None:
            pw.finish()

    def _save(self, res, spec, out, name, tag):
        """The spectrum, the engine input, the mass spectrum and the mass table next to the data, under
        a name of this result's own (an earlier result's files are never written over); the id written
        into _peaks.csv tells later edits whether the files are still this result's."""
        import uuid
        try:
            base = _unique_base(out, name, tag)
            res["saved_id"] = uuid.uuid4().hex
            # the spectrum as text (6 decimals: 0.001 ppm at m/z 1000, a quarter of
            # the size and time of the full precision export); the engine input
            # keeps every digit
            _savetxt(os.path.join(out, name + ".txt"), spec, "%.6f\t%.4f")
            inp = res.get("input_data", spec)
            _savetxt(base + "_input.txt", inp, "%.17g" if inp is not spec else "%.6f\t%.4f")
            if res.get("mass") is not None and len(res["mass"]):
                _savetxt(base + "_mass.txt", res["mass"], "%.17g", header="mass (Da)\tintensity (data scale)")
            self._write_peaks(res, base + "_peaks.csv")
            res["saved"] = base
        except Exception as ex:  # a failed save must not leave the panel busy
            res["saved"] = None
            res["save_error"] = str(ex)

    @staticmethod
    def _write_peaks(res, path):
        """The mass table as CSV: for isotope resolved species the most
        abundant isotope (as labelled in the plot), else the charges (IsoDec)
        or the charge states in range; the score of each mass when the
        method gives one."""
        grouped, scores = _grouped(res), _has_scores(res)
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("# %s\n# input: %s\n" % (res.get("notes", ""), res.get("input", "")))
            if res.get("saved_id"):
                fh.write("# result id: %s\n" % res["saved_id"])
            fh.write(("envelope centroid (Da),intensity (%),share (%),selected isotope peak (Da)" if grouped else
                      "mass (Da),intensity (%),share (%),charges") + (",score" if scores else "")
                     + "\n")
            for q in res["peaks"]:
                if grouped and "apex" in q:
                    extra = "%.4f" % _apex_mass(res, q)
                else:
                    extra = q.get("z", "")
                    if extra in ("", None) and res.get("tag") != "isodec":
                        extra = _z_text_of(q["mass"], res)
                if q.get("added"):
                    extra = (str(extra) + " (added by hand)").strip()
                line = "%.5f,%.2f,%.2f,%s" % (q["mass"], q["rel"], q["frac"], str(extra).replace(",", " "))
                if scores:
                    s = _peak_score(q)
                    line += "," + ("%.3f" % s if s is not None else "")
                fh.write(line + "\n")

    @staticmethod
    def _rewrite_peaks(res):
        """The mass table next to the data written again after an edit (a mass added or removed, undo),
        only while the files there are still this result's (its id in _peaks.csv). True when written;
        False when the result was not saved, its files are gone, or another result wrote over them
        (results of older versions carry no id: their files are left alone)."""
        base, rid = res.get("saved"), res.get("saved_id")
        if not base or not rid:
            return False
        path = base + "_peaks.csv"
        try:
            with open(path, encoding="utf-8", errors="replace") as fh:
                head = [fh.readline().strip() for _ in range(4)]
        except OSError:
            return False
        if "# result id: %s" % rid not in head:
            return False
        DeconvPanel._write_peaks(res, path)
        return True

    def _done(self, res, err, own=None):
        self.busy = False
        self._end_progress()
        self._refresh_files()
        try:
            self.run_btn.Enable(True)
        except RuntimeError:
            return
        if err:
            self.status("Deconvolution failed: " + err)
            wx.MessageBox("The deconvolution failed:\n\n%s" % err, "Deconvolute", wx.ICON_ERROR)
            return
        mass = res.get("mass")
        if mass is not None and len(mass) > 1 and np.any(np.diff(mass[:, 0]) < 0):
            res["mass"] = mass[np.argsort(mass[:, 0], kind="stable")]  # markers and labels search a rising axis
        # a new row above the older results; the oldest goes when there are too many
        ent = DeconvResult(res)
        ent.input = own
        ent.sel = int(np.argmax([q["height"] for q in res["peaks"]])) if res["peaks"] else None
        self._make_row(ent)
        self.results.insert(0, ent)
        dropped = []
        while len(self.results) > MAX_RESULTS:
            dropped.append(self.results.pop())
        old = self.active
        self.active = ent
        self._layout()
        for gone in dropped:
            self._drop_row(gone)
        self.plot_results(ent)
        self._mark_rows()
        self._replot_spec({res["e"]} | ({old.res["e"]} if old is not None else set()))
        self.tab.reveal(ent.row)
        self._notify()
        msg = "%s in %.1f s" % (res["notes"], res["seconds"])
        if res.get("r2") is not None:
            msg += "; R² %.3f" % res["r2"]
        if res.get("saved"):
            msg += "; saved as %s_*" % os.path.basename(res["saved"])
        elif res.get("save_error"):
            msg += "; SAVE FAILED: " + res["save_error"]
            wx.MessageBox("The result is available, but could not be saved:\n\n" + res["save_error"],
                          "Deconvolution export", wx.OK | wx.ICON_WARNING)
        if dropped:
            msg += "; the oldest result row was removed (%d are kept)" % MAX_RESULTS
        self.status(msg)
        self._set_note("%s · %d masses · %.1f s" % (res["method"], len(res["peaks"]), res["seconds"]))

    # ---------------------------------------------------------------- plots
    def _cols(self, iso, resolved=False, scores=False, added=False):
        # wider number column with masses added by hand ("2 (added)" was cut to "2 (...")
        cols = [("#", 78 if added else 34), ("Envelope mass (Da)" if resolved else "Mass (Da)", 170 if resolved else 124),
                ("Intensity", 96), ("Share %", 66)]
        if iso:
            cols += [("Charges", 110), ("m/z", 92)]
        elif resolved:
            cols += [("Selected isotope (Da)", 180), ("Isotope peaks", 96)]
        else:
            cols += [("Charge states in range", 170)]
        if scores:
            cols.append(("Score", 64))
        return cols

    def overlay(self, col, ax):
        """Fit and charge states of the selected mass of the active result
        on the spectrum."""
        res = self.result
        if res is None or col["e"] != res["e"]:
            return
        key = res.get("spec_key")
        if key is not None and key != _spec_key(col["e"], col.get("spec")):
            # another spectrum is shown now (another range or scan): the fit and charge states of the
            # result do not belong to it
            t = ax.text(0.01, 0.98, "Fit not shown: result of another spectrum (%s)" % res.get("input", ""),
                        transform=ax.transAxes, ha="left", va="top",
                        fontsize=7.5, color="#8A8F98", zorder=4)
            t._no_scale = True
            return
        fit = res.get("fit")
        if fit is not None and len(fit):
            ln = U.plot_trace(ax, fit[:, 0], fit[:, 1], color="#E0602F", lw=U.LW_TRACE, alpha=0.85, zorder=3)
            ln._no_scale = True
        if self.sel is None or res.get("tag") == "isodec" or self.sel >= len(res["peaks"]):
            return
        q = res["peaks"][self.sel]
        colr = U.PAL[self.sel % len(U.PAL)]
        from matplotlib.transforms import blended_transform_factory
        tr = blended_transform_factory(ax.transData, ax.transAxes)
        x0, x1 = ax.get_xlim()
        gap = (x1 - x0) / 45.0
        last = None
        for z, mz in sorted(self._charges_in_range(q["mass"], res), key=lambda t: t[1]):
            if not (x0 <= mz <= x1):
                continue
            ln = ax.axvline(mz, color=colr, lw=0.8, ls=(0, (3, 2)), alpha=0.7, zorder=1)
            ln._no_scale = True
            if last is None or mz - last >= gap:
                # charge tags just above the frame, where they cannot cover the
                # m/z labels of the peaks (near the left edge, where the
                # exponent of the axis sits, inside the frame instead)
                inside = (mz - x0) / ((x1 - x0) or 1.0) < 0.1
                ax.text(mz, 0.97 if inside else 1.012, "%d+" % z if res["params"].get("sign", 1) > 0 else "%d−" % z,
                        transform=tr, ha="center", va="top" if inside else "bottom", fontsize=7.5, color=colr,
                        zorder=6, clip_on=False,
                        bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.85))
                last = mz

    def _subtitle(self, ent, mark=True, quality=True):
        """Method, spectrum and mass range of a result (the rows are told
        apart by it), (quality) why a result is poor, and (mark) whether it
        is the one shown in the table."""
        res = ent.res
        what = res.get("input", "")
        if len(self.tab.cols) <= 1 and ": " in what:
            what = what.split(": ", 1)[1]  # one scan event: its name says nothing
        parts = [res.get("method", ""), what]
        mr = res.get("params", {}).get("mass_range")
        if mr:
            parts.append("%g to %g Da" % tuple(mr))
        n = len(res.get("peaks") or [])
        txt = ", ".join(p for p in parts if p) + " · %d mass%s" % (n, "es" if n != 1 else "")
        if res.get("r2") is not None:
            txt += ", R² %.3f" % res["r2"]
        qa = _quality(res) if quality else None
        if qa and qa["level"] == "poor" and qa["text"]:
            txt += " · " + qa["text"]  # the reason, after the "Poor result" badge
        if mark and ent is self.active and len(self.visible_results()) > 1:
            txt += " · shown in the table"
        return txt

    def _badge(self, ent):
        """(text, colour, tooltip) of the warning shown with a poor result, or
        None (no quality given, or good or fair)."""
        qa = _quality(ent.res)
        if qa and qa["level"] == "poor":
            return "Poor result", T.GROUP["red"], _quality_text(qa)
        return None

    def _mark_rows(self):
        """Subtitles of every result row (the active one is marked) and the
        warning of a poor result."""
        for ent in self.results:
            if ent.card_mass is not None:
                iso = ent.res.get("tag") == "isodec"
                ent.card_mass.set_title("Monoisotopic masses" if iso else "Deconvoluted mass spectrum",
                                        self._subtitle(ent))
                ent.card_mass.set_badge(*(self._badge(ent) or (None,)))

    def tree_text(self, ent):
        """What the list of open files shows for a result: (first line,
        second line, tooltip). First line: method and masses; second line:
        spectrum (time range) and mass range."""
        res = ent.res
        name = {"Maximum entropy": "Max. entropy"}.get(res.get("method", ""), res.get("method", "") or "Result")
        n = len(res.get("peaks") or [])
        line1 = "%s, %d mass%s" % (name, n, "es" if n != 1 else "")
        what = res.get("input", "")
        if ": " in what:
            ev, what = what.split(": ", 1)
            if len(self.tab.cols) > 1:  # two scan events (e.g. positive and negative): say which
                what = "%s %s" % (ev, what)
        mr = res.get("params", {}).get("mass_range")
        line2 = " · ".join(p for p in (what, ("%g to %g Da" % tuple(mr)) if mr else "") if p)
        tip = [self._subtitle(ent, False, False)]
        if ent is self.active:
            tip.append("Shown in the table")
        elif not ent.visible:
            tip.append("Hidden: not drawn (kept)")
        qa = _quality(res)
        if qa:
            tip.append(_quality_text(qa))
        return line1, line2, "\n".join(tip)

    def _fill_table(self, ent):
        """The mass table shows ent (columns, rows, the selected peak and the
        warning of a poor result)."""
        if self.table is None:
            return
        res = ent.res
        iso = res.get("tag") == "isodec"
        self.table.set_columns(self._cols(iso, _grouped(res), _has_scores(res),
                                          any(q.get("added") for q in res.get("peaks") or [])))
        self.table.set_rows(ent.rows, "%d masses · %s" % (len(ent.rows), res.get("method", "")))
        details = getattr(self.tab, "hrms_details", None)
        if details is not None:
            details.refresh_result(ent)
        if ent.sel is not None:
            self.table.select(ent.sel)
        badge, qa = self._badge(ent), _quality(res)
        self.table.set_badge(badge[0] if badge else None, badge[1] if badge else None, _quality_text(qa))

    def plot_results(self, ent=None):
        """Draw the plots of a result (the active one by default); the
        table follows when it is the active one."""
        ent = self.active if ent is None else ent
        if ent is None or ent.row is None:
            return
        res = ent.res
        iso = res.get("tag") == "isodec"
        resolved = _grouped(res)
        # charge states (heights in the spectrum); drawn also while the tile is
        # hidden (cheap: the canvas draws when shown), so it is ready to show
        cz = ent.card_z
        cz.reset()
        U.style_axes(cz.ax, "Charge", "Intensity")
        zd = res.get("zdist")
        if zd is not None and len(zd) and np.nanmax(zd[:, 1]) > 0:
            zz, zi = zd[:, 0], zd[:, 1]
            cz.ax.bar(zz, zi, width=0.7, color="#2A62C4", alpha=0.85)
            from matplotlib.ticker import MaxNLocator
            cz.ax.xaxis.set_major_locator(MaxNLocator(nbins=6, integer=True))
            cz.extra_scale = [(zz, zi)]
            zrange = res.get("params", {}).get("z_range")
            cz.full = (float(zrange[0]) - 1, float(zrange[1]) + 1) if zrange else (
                float(zz.min()) - 1, float(zz.max()) + 1)
            cz.ax.set_xlim(*cz.full)
            cz.autoscale_y(pad=0.12)
        else:
            cz.full = None
            cz.ax.set_xticks([])
            cz.ax.set_yticks([])
            cz.ax.text(0.5, 0.5, "see the table" if iso else "no charge data", transform=cz.ax.transAxes,
                       ha="center", va="center", fontsize=8.5, color=U.MUTED)
        cz.set_title("Charge states", "%s, tallest peak at each charge state" % res.get("method", ""))
        cz.draw()
        # zero charge mass spectrum
        cm = ent.card_mass
        cm.reset()
        U.style_axes(cm.ax, "Mass (Da)", "Intensity")
        mass = res.get("mass")
        cm.full = self._mass_range(res)
        if mass is None or not len(mass):
            if cm.full is not None:
                cm.ax.set_xlim(*cm.full)
            cm.ax.text(0.5, 0.5, "No masses found", transform=cm.ax.transAxes, ha="center", va="center",
                       fontsize=9, color=U.MUTED)
            cm.set_title("Monoisotopic masses" if iso else "Deconvoluted mass spectrum", self._subtitle(ent))
            cm.set_badge(*(self._badge(ent) or (None,)))
            ent.rows = []
            if ent is self.active:
                self._fill_table(ent)
            cm.draw()
            return
        step = self._fmt_step(res)
        peaks = res["peaks"]
        line_col = U.TRACE
        if iso:
            cm.ax.vlines([q["mass"] for q in peaks], 0, [q["height"] for q in peaks], color=line_col, lw=0.9)
            cm.extra_scale = [(np.array([q["mass"] for q in peaks]), np.array([q["height"] for q in peaks]))]
        else:
            U.plot_trace(cm.ax, mass[:, 0], mass[:, 1], fill=("#2A62C4", 0.14), color=line_col, lw=U.LW_TRACE,
                         zorder=3)
        # first shown: the whole mass range that was deconvoluted (right click >
        # Zoom to the masses found narrows it); a zoom of the user is kept
        cm.ax.set_xlim(*(ent.view if ent.view else cm.full))
        # room above the tallest peak for its label (5 pt above the marker, 8 pt bold; 16 pt touched the
        # frame) and for resolved species the average mass above it
        cm.label_room_pt = 30.0 if resolved else 20.0
        if res.get("label_all") and not iso:
            cm.label_room_pt += 21.0  # Label all peaks: two stacked labels above the tallest peaks stay inside
        cm.label_room_pt += self.shifts.add_room(ent)  # room for the brackets of the mass shifts
        cm.autoscale_y(pad=0.16)
        rows = self._draw_labels(cm, res, peaks, mass, step, iso, resolved, ent.sel)
        _labels_inside(cm.ax)
        cm.set_title("Monoisotopic masses" if iso else "Deconvoluted mass spectrum", self._subtitle(ent))
        cm.set_badge(*(self._badge(ent) or (None,)))
        ent.rows = rows
        if ent is self.active:
            self._fill_table(ent)
        cm.draw()

    @staticmethod
    def _marker_xy(q, mass, iso, resolved):
        """Where the marker of a mass sits: on the top of its peak as drawn
        (isotope resolved species: of its most abundant isotope peak);
        IsoDec: on the top of its stick."""
        if iso or mass is None or not len(mass):
            return float(q["mass"]), float(q["height"])
        if resolved:
            return _drawn_top(mass, _display_apex(q), 0.3)
        # maximum entropy gives the grid point of the top as "apex": start there
        return _drawn_top(mass, float(q.get("apex", q["mass"])))

    def _draw_labels(self, cm, res, peaks, mass, step, iso, resolved, sel=None):
        """Markers and mass labels of the peaks in the current view (drawn
        again after each zoom, without drawing the whole plot again); sel is
        the index of the selected peak. Returns the rows of the mass table.

        A marker sits on the top of its peak as drawn (isotope resolved
        species: of the selected isotope). With isotope peaks, every
        isotope peak in view is labelled with its mass (4 decimals, as
        DataAnalysis shows it); the selected species (and any labelled by
        hand) gets a dot on its selected isotope, whose label is then in
        bold. The envelope centroid is reported separately in the table;
        it is not the position of an individual isotope peak. Other results label their
        tallest masses (bold: the selected one). Labels do not overlap, on
        screen nor in the (smaller) copied image; the selected and pinned
        labels are placed first, the others where there is room."""
        for a in list(cm.ax.texts) + [ln for ln in cm.ax.lines if getattr(ln, "_dec_mark", False)]:
            if getattr(a, "_dec_label", False) or getattr(a, "_dec_mark", False):
                a.remove()
        try:
            n_lab = int(res["params"].get("labels", 12))
        except (TypeError, ValueError):
            n_lab = 12
        isop = [] if iso else list(res.get("isotope_peaks") or [])
        x0v, x1v = cm.ax.get_xlim()
        y0v, y1v = cm.ax.get_ylim()
        wv, hv = (x1v - x0v) or 1.0, (y1v - y0v) or 1.0
        # the labels are placed for the plot on screen: a journal sized copied image is
        # much narrower, and placing them for it left a zoomed view with a single label
        sw, sh = _label_room(cm.ax, screen=True)
        # isotope peaks get their own labels only when the view is narrow enough for
        # them to stand side by side (1.003 Da at least three quarters of a 10 character
        # label); in a wide view every isotope label would sit in one column above
        # the species, so the species are labelled instead (most abundant isotope)
        label_all = bool(res.get("label_all"))
        # Label all peaks on a coarse mass axis (envelopes, 0.1 Da steps and more): a selected
        # or pinned mass keeps its own label (the mass of the table, as its menu said) and
        # takes the place of the label of its drawn peak, whose top lies up to a couple of
        # steps away from it (the peak was labelled twice, e.g. 12244.3 and 12244.8)
        coarse_all = label_all and not iso and not resolved and \
            float(res.get("params", {}).get("mass_step", 1.0) or 1.0) >= 0.1
        try:  # peaks closer than the peak window are one peak of the result
            win = float(res.get("params", {}).get("peak_window", 0) or 0)
        except (TypeError, ValueError):
            win = 0.0
        if label_all and not iso:
            # right click > Label all peaks in view: every peak of the view, from the drawn curve
            isop = _view_maxima(mass, x0v, x1v, step)
            n_lab = max(n_lab, len(isop) + len(peaks))
        elif label_all:
            n_lab = max(n_lab, len(peaks))  # IsoDec: every mass in view (the menu item did nothing)
        elif isop and sw * 1.00335 / wv < 0.45 * 10 * 0.6 * 8.0:
            # too narrow for every other isotope label to fit: label the species
            isop = []
        hfrac = 1.3 * 8.0 / sh  # height of a label line in axes fractions

        def half(text, size=8.0):  # half the width of a label in axes fractions
            return 0.5 * len(text) * 0.6 * size / sw

        def clash(fx, fy, hw, nl, boxes):
            return any(abs(fx - bx) <= hw + bh and fy < by + bn * hfrac and by < fy + nl * hfrac
                       for bx, by, bh, bn in boxes)

        top = [self._marker_xy(q, mass, iso, resolved) for q in peaks]
        boxes = []  # (x, y, half width, lines) of the labels placed, axes fractions
        labels = []  # (peak index or None, x, y, text, bold, second line, lift in lines)
        special = ([sel] if sel is not None and 0 <= sel < len(peaks) else []) + \
            [n for n, q in enumerate(peaks) if q.get("pinned") and n != sel]
        used_iso = set()
        for n in special:
            q = peaks[n]
            xm, ym = top[n]
            k = _nearest_iso(isop, _display_apex(q) if resolved else float(q["mass"]),
                             max(0.3, 2.0 * step, 0.5 * win) if coarse_all else 0.3)
            second = None
            if resolved:
                text = "%.4f" % (isop[k][0] if k is not None else _display_apex(q))
                second = None
            elif k is not None and not coarse_all:
                text = "%.4f" % isop[k][0]
            else:
                text = _mass_fmt(q["mass"], step)
            if k is not None:
                used_iso.add(k)
            fx, fy = (xm - x0v) / wv, (ym - y0v) / hv
            nl = 2 if second else 1
            hw = max(half(text), half(second, 7.0) if second else 0.0)
            lift = 0
            while lift < 6 and clash(fx, fy + lift * hfrac, hw, nl, boxes):
                lift += 1
            boxes.append((fx, fy + lift * hfrac, hw, nl))
            labels.append((n, xm, ym, text, n == sel, second, lift))
        count = 0
        if isop:
            for k in sorted(range(len(isop)), key=lambda j: -isop[j][1]):
                if count >= n_lab:
                    break
                m, h = isop[k]
                if k in used_iso or not x0v <= m <= x1v:
                    continue
                fine = float(res.get("params", {}).get("mass_step", 1.0) or 1.0) < 0.1
                text = ("%.4f" % m) if (fine or not label_all) else _mass_fmt(m, step)
                fx, fy = (m - x0v) / wv, (h - y0v) / hv
                lift = 0
                def below_top(lf):  # a stacked label stays inside the frame (it was drawn across its top)
                    return lf == 0 or fy + (13.4 + 10.4 * lf) / sh <= 1.0
                while label_all and lift < 4 and clash(fx, fy + lift * hfrac, half(text), 1, boxes) and \
                        below_top(lift + 1):
                    lift += 1  # all peaks: crowded labels are stacked instead of left out
                if not clash(fx, fy + lift * hfrac, half(text), 1, boxes):
                    boxes.append((fx, fy + lift * hfrac, half(text), 1))
                    labels.append((None, m, h, text, False, None, lift))
                    count += 1
        else:
            for n in sorted(range(len(peaks)), key=lambda j: -peaks[j]["height"]):
                if count >= n_lab:
                    break
                if n in special:
                    continue
                xm, ym = top[n]
                fx, fy = (xm - x0v) / wv, (ym - y0v) / hv
                if not 0.0 <= fx <= 1.0:
                    continue
                q = peaks[n]
                text = ("%.4f" % _display_apex(q)) if resolved else _mass_fmt(q["mass"], step)
                if not clash(fx, fy, half(text), 1, boxes):
                    boxes.append((fx, fy, half(text), 1))
                    labels.append((n, xm, ym, text, False, None, 0))
                    count += 1
        cm._dec_marks = []  # (peak index, x, y) of the markers, for the tooltip
        for n, q in enumerate(peaks):
            is_sel = n == sel
            if is_sel or not isop or q.get("pinned"):
                xm, ym = top[n]
                mk = cm.ax.plot([xm], [ym], ls="none", marker="o", ms=6 if is_sel else 4.5,
                                mfc=U.PAL[n % len(U.PAL)], mec="white" if not is_sel else U.INK, mew=0.8,
                                zorder=5)[0]
                mk._no_scale = mk._dec_mark = True
                cm._dec_marks.append((n, xm, ym))
        for n, xm, ym, text, bold, second, lift in labels:
            dy = (3 if n is None else 5) + 10.4 * lift  # isotope labels sit closer to their peak
            cm.ax.annotate(text, (xm, ym), xytext=(0, dy), textcoords="offset points", ha="center", va="bottom",
                           fontsize=8, color=U.INK, fontweight="bold" if bold else "normal",
                           zorder=6)._dec_label = True
            if second:
                cm.ax.annotate(second, (xm, ym), xytext=(0, dy + 10.4), textcoords="offset points", ha="center",
                               va="bottom", fontsize=7, color="#444444", zorder=6)._dec_label = True
        ent_ = next((e for e in self.results if e.card_mass is cm), None)
        if ent_ is not None:
            self.shifts.draw(ent_)  # brackets of the mass shifts above the labels
        rows = []
        scores = _has_scores(res)
        for n, q in enumerate(peaks):
            num = "%d (added)" % (n + 1) if q.get("added") else n + 1
            if iso:
                row = (num, _mass_fmt(q["mass"], step), U._g(q["height"]), "%.1f" % q["frac"], q.get("z", ""),
                       "%.4f" % q.get("mz", 0))
            elif resolved:
                row = (num, "%.3f" % q["mass"], U._g(q["height"]), "%.1f" % q["frac"],
                       "%.4f" % _apex_mass(res, q) + _apex2_suffix(q), str(q.get("n_iso", "")))
            else:
                row = (num, _mass_fmt(q["mass"], step), U._g(q["height"]), "%.1f" % q["frac"],
                       self._z_text(q["mass"], res))
            if scores:
                s = _peak_score(q)
                row += ("%.2f" % s if s is not None else "",)
            rows.append(row)
        return rows

    def _relabel_now(self, ent):
        """Markers and labels of a result placed again for the present size of its plot (an image of
        another size: report, copied or saved image), at once."""
        res, cm = ent.res, ent.card_mass
        mass, peaks = res.get("mass"), res.get("peaks")
        if cm is None or mass is None or not len(mass) or not peaks:
            return
        self._draw_labels(cm, res, peaks, mass, self._fmt_step(res), res.get("tag") == "isodec", _grouped(res),
                          ent.sel)
        _labels_inside(cm.ax)

    def report_results(self, e=None):
        """The results drawn on screen (visible rows, newest first, as in the tile column), of scan
        event e (None: every event)."""
        out = []
        for ent in self.results:
            if not ent.visible or ent.card_mass is None:
                continue
            if e is not None and ent.res.get("e") != e:
                continue
            out.append(ent)
        return out

    def shift_state(self):
        """The mass shift finder as a plain dict (settings; per result on or off and the reference)."""
        return self.shifts.get_state()

    def set_shift_state(self, d):
        """Back to a dict of shift_state (every result drawn again)."""
        self.shifts.set_state(d)

    def charge_tile_shown(self, ent):
        try:
            return ent.card_z is not None and ent.row is not None and ent.card_z in ent.row.panes and \
                ent.card_z not in getattr(ent.row, "closed", set())
        except RuntimeError:
            return False

    def report_images(self, ent, folder, tag, w_cm=17.8, h_cm=6.0, text_scale=0.75, with_z=True, h_range=None,
                      min_scale=0.6):
        """Images of a result for a report, exactly as its row shows it: the zero charge mass spectrum
        (same mass range, curve, fill, markers and labels; the labels placed again for the size in
        the report) and the charge states when that tile is shown. Returns {"mass": path,
        "z": path or None, "w_mass", "w_z" (cm), "title", "subtitle", "badge" (text or None)}."""
        cm, res = ent.card_mass, ent.res
        show_z = with_z and self.charge_tile_shown(ent)
        wz = 0.0
        if show_z:
            try:
                wm_, wz_ = ent.row.weight_of(ent.card_mass, 3.4), ent.row.weight_of(ent.card_z, 1.0)
                wz = max(3.6, min(7.0, (w_cm - 0.2) * wz_ / (wm_ + wz_)))
            except Exception:
                wz = 4.4
        wm = w_cm - (wz + 0.2 if show_z else 0.0)

        def relabel():
            self._relabel_now(ent)
        # the height follows the proportions of the tile on screen (h_range): the labels then stack as there
        try:
            ws, hs = cm.fig.get_size_inches()
            if h_range is not None and ws > 0.4 and hs > 0.2:
                h_cm = max(h_range[0], min(h_range[1], wm * hs / ws))
        except Exception:
            pass
        out = {"mass": cm.snapshot(os.path.join(folder, "%s_mass.png" % tag), wm, h_cm, text_scale, relabel=relabel,
                                   min_scale=min_scale),
               "z": None, "w_mass": wm, "w_z": wz}
        if show_z:
            out["z"] = ent.card_z.snapshot(os.path.join(folder, "%s_z.png" % tag), wz, h_cm, text_scale,
                                           min_scale=min_scale)
        iso = res.get("tag") == "isodec"
        out["title"] = "Monoisotopic masses" if iso else "Deconvoluted mass spectrum"
        out["subtitle"] = self._subtitle(ent, mark=False)
        b = self._badge(ent)
        out["badge"] = b[0] if b else None
        return out

    def _mass_tip(self, ent, x):
        """Tooltip over the zero charge mass spectrum near a marker."""
        cm, res = ent.card_mass, ent.res
        marks = getattr(cm, "_dec_marks", None) if cm is not None else None
        peaks = res.get("peaks") or []
        if not marks or x is None:
            return ""
        x0, x1 = cm.ax.get_xlim()
        n, xm, ym = min(marks, key=lambda t: abs(t[1] - x))
        if abs(xm - x) > abs(x1 - x0) / 80.0 or n >= len(peaks):
            return ""
        q = peaks[n]
        if _grouped(res) and "apex" in q:
            lines = ["Selected peak %.4f Da (the dot)" % _apex_mass(res, q),
                     "Envelope centroid %.3f Da" % q["mass"] +
                     (", %s isotope peaks" % q["n_iso"] if q.get("n_iso") else "")]
        else:
            lines = ["%s Da" % _mass_fmt(q["mass"], self._fmt_step(res))]
            z = q.get("z") if res.get("tag") == "isodec" else self._z_text(q["mass"], res)
            if z not in ("", None):
                lines.append("Charge%s %s" % ("" if res.get("tag") == "isodec" else " states in range", z))
        lines.append("%.1f %% of the largest, %.1f %% of the total" % (q.get("rel", 0.0), q.get("frac", 0.0)))
        s = _peak_score(q)
        if s is not None:
            lines.append("Score %.2f" % s)
        if q.get("added"):
            lines.append("Added by hand")
        return "\n".join(lines)

    @staticmethod
    def _mass_range(res):
        """The requested plotting domain, even if the engine returns less data."""
        requested = res.get("params", {}).get("mass_range")
        if requested is not None:
            lo, hi = map(float, requested)
            if np.isfinite([lo, hi]).all() and hi > lo:
                return lo, hi
        mass = res.get("mass")
        if mass is None or not len(mass):
            return None
        lo, hi = float(np.min(mass[:, 0])), float(np.max(mass[:, 0]))
        if hi <= lo:
            pad = max(float(res.get("params", {}).get("mass_step", 1)), 1.0)
            return lo - pad, hi + pad
        return lo, hi

    def _peak_view(self, res):
        """Mass range around the masses found."""
        peaks = sorted(res["peaks"], key=lambda q: q["mass"])
        full = self._mass_range(res)
        resolved = _grouped(res)
        if not peaks:
            return full
        span = peaks[-1]["mass"] - peaks[0]["mass"]
        pad = max(span * 0.15, (3 if resolved else 10) * res["params"].get("peak_window", 10),
                  5 * res["params"]["mass_step"], 3.0)
        lo = max(full[0], peaks[0]["mass"] - pad)
        hi = min(full[1], peaks[-1]["mass"] + pad)
        return (lo, hi) if hi > lo else full

    def zoom_peaks(self, ent=None):
        ent = self.active if ent is None else ent
        if ent is not None and ent.row is not None:
            ent.view = self._peak_view(ent.res)
            ent.card_mass.ylock = None
            self.plot_results(ent)

    def _mass_view_changed(self, ent=None):
        """Labels follow the zoom of a mass spectrum: only the labels and
        markers of that row are drawn again, once the zooming has paused
        briefly."""
        ent = self.active if ent is None else ent
        if ent is None or ent.row is None:
            return
        ent.view = tuple(ent.card_mass.ax.get_xlim())
        if ent not in self._relabel_pending:
            self._relabel_pending.append(ent)
        t = getattr(self, "_relabel_timer", None)
        if t is not None and t.IsRunning():
            t.Restart(120)
        else:
            self._relabel_timer = wx.CallLater(120, self._relabel)

    def _relabel(self):
        pending, self._relabel_pending = self._relabel_pending, []
        for ent in pending:
            res = ent.res
            if ent not in self.results or ent.row is None or not res.get("peaks") or res.get("mass") is None:
                continue
            try:
                cm = ent.card_mass
                iso = res.get("tag") == "isodec"
                self._draw_labels(cm, res, res["peaks"], res["mass"], self._fmt_step(res), iso, _grouped(res), ent.sel)
                _labels_inside(cm.ax)
                cm.canvas.draw_idle()
            except RuntimeError:
                pass  # window closed

    def _charges_in_range(self, mass, res=None):
        return _charges_of(mass, self.result if res is None else res)

    def _z_text(self, mass, res=None):
        return _z_text_of(mass, self.result if res is None else res)

    def _replot_spec(self, events=None):
        """The spectra of the scan events given (default: the one of the
        active result) drawn again, with the overlay of the active result."""
        if events is None:
            events = {self.result["e"]} if self.result else set()
        for col in self.tab.cols:
            if col["e"] in events:
                self.tab.plot_spec(col, keep_view=True)

    def select_peak(self, i):
        if not self.result or i == self.sel:
            return
        self.sel = i
        self.plot_results()
        self._replot_spec()

    # ------------------------------------------------- masses added by hand
    # (on the active result)
    def _peak_near(self, x):
        res = self.result
        if not res or not res["peaks"] or x is None:
            return None
        hit = _isotope_selection(res, x)
        if hit is not None:
            return hit[0]
        x0, x1 = self.card_mass.ax.get_xlim()
        d = [abs(self._pick_x(res, q) - x) for q in res["peaks"]]
        i = int(np.argmin(d))
        return i if d[i] <= (x1 - x0) / 40.0 else None

    @staticmethod
    def _pick_x(res, q):
        """Where a mass is clicked: resolved species at their most abundant
        isotope (where the dot is), other masses at their mass."""
        return _display_apex(q) if _grouped(res) else float(q["mass"])

    def _local_max(self, x):
        """Mass and height of the highest point of the mass spectrum near x."""
        res = self.result
        mass = res.get("mass") if res else None
        if mass is None or not len(mass) or x is None:
            return None
        x0, x1 = self.card_mass.ax.get_xlim()
        w = (x1 - x0) / 40.0
        a, b = np.searchsorted(mass[:, 0], [x - w, x + w])
        if b <= a:
            return None
        k = a + int(np.argmax(mass[a:b, 1]))
        if mass[k, 1] <= 0:
            return None
        if 0 < k < len(mass) - 1 and mass[k - 1, 1] > 0 and mass[k + 1, 1] > 0:
            # apex between the points (parabola through the top three)
            y0, y1, y2 = mass[k - 1, 1], mass[k, 1], mass[k + 1, 1]
            den = y0 - 2 * y1 + y2
            off = 0.5 * (y0 - y2) / den if den < 0 else 0.0
            off = max(-0.5, min(0.5, off))
            step = mass[k + 1, 0] - mass[k, 0] if off >= 0 else mass[k, 0] - mass[k - 1, 0]
            return float(mass[k, 0] + off * step), float(y1)
        return float(mass[k, 0]), float(mass[k, 1])

    def _mass_menu(self, x):
        """Right click on the mass spectrum: label the mass there, or add it
        to the masses (when no mass was found there), or remove it again."""
        res = self.result
        if not res or res.get("mass") is None or x is None:
            return []
        step = self._fmt_step()
        i = self._peak_near(x)
        if i is not None:
            q = res["peaks"][i]
            txt = ("%.4f" % _apex_mass(res, q)) if _grouped(res) and "apex" in q else _mass_fmt(q["mass"], step)
            if q.get("added"):
                return [("Remove the added mass %s Da" % txt, lambda: self.remove_mass(i)), (None, None)]
            if q.get("pinned"):
                return [("Remove the label of %s Da" % txt, lambda: self.pin_mass(i, False)), (None, None)]
            return [("Label %s Da (keep the label)" % txt, lambda: self.pin_mass(i, True)), (None, None)]
        if res.get("tag") == "isodec":
            return []
        m = self._local_max(x)
        if m is None:
            return []
        return [("Add the mass here: %s Da (label and table)" % _mass_fmt(m[0], step),
                 lambda: self.add_mass(m[0], m[1])), (None, None)]

    def _masses_changed(self, msg):
        res = self.result
        D._rel_frac(res["peaks"])
        self.plot_results()
        self._replot_spec()
        if res.get("saved"):
            try:
                if not self._rewrite_peaks(res):
                    msg += "; %s_peaks.csv not updated (the file is not this result's)" % os.path.basename(
                        res["saved"])
            except OSError:
                pass
        self.status(msg)
        self._notify()  # the number of masses in the list of open files

    def pin_mass(self, i, flag):
        q = self.result["peaks"][i]
        q["pinned"] = bool(flag)
        self.sel = i
        self.plot_results()
        self._replot_spec()

    def add_mass(self, m, h):
        res = self.result
        peaks = res["peaks"]
        q = {"mass": float(m), "height": float(h), "area": 0.0, "score": 0.0, "added": True, "pinned": True}
        if _grouped(res):
            q.update(apex=float(m), n_iso="")
        ms = [p["mass"] for p in peaks]
        k = int(np.searchsorted(ms, m)) if ms == sorted(ms) else len(peaks)
        peaks.insert(k, q)
        self.sel = k
        self._masses_changed("Added %s Da (height %s) to the masses" % (_mass_fmt(m, self._fmt_step()), U._g(h)))

    def remove_mass(self, i):
        res = self.result
        q = res["peaks"].pop(i)
        if self.sel is not None:
            self.sel = None if self.sel == i else (self.sel - 1 if self.sel > i else self.sel)
        self._masses_changed("Removed %s Da" % _mass_fmt(q["mass"], self._fmt_step()))

    def pick_mass(self, x, y):
        res = self.result
        if not res or not res["peaks"]:
            return
        x0, x1 = self.card_mass.ax.get_xlim()
        d = [abs(self._pick_x(res, q) - x) for q in res["peaks"]]
        i = int(np.argmin(d))
        hit = _isotope_selection(res, x, min(0.3, max(0.02, (x1 - x0) / 40.0)))
        if hit is not None:
            i, m = hit
            res["peaks"][i]["selected_apex"] = m
        if hit is not None or d[i] <= (x1 - x0) / 40.0:
            self.sel = i
            self.plot_results()
            self._replot_spec()
            q = res["peaks"][i]
            if _grouped(res) and "apex" in q:
                self.status("Selected isotope %.4f Da; centroid %.3f Da%s: %.1f %% of the largest species" % (
                    _apex_mass(res, q), q["mass"], ", %s isotope peaks" % q["n_iso"] if q.get("n_iso") else "",
                    q["rel"]))
            else:
                self.status("%s Da: %.1f %% of the largest mass; charge states %s" % (
                    _mass_fmt(q["mass"], self._fmt_step()), q["rel"], q.get("z") or self._z_text(q["mass"])))

    # ------------------------------------------------------------- output
    def open_full(self, d=None):
        """The spectrum in the full Deconvolute window (own process)."""
        if d is None:
            self.refresh_inputs(quiet=True)
            e = self.result["e"] if self.result else None
            d = next((x for x in self.inputs if x.get("e") == e), self.inputs[0] if self.inputs else None)
        if d is None:
            self.status("No spectrum: average a time range on a chromatogram first")
            return
        if hasattr(self.tab, "open_full") and d.get("e") is not None:
            # as right click > Open this spectrum in the Deconvolute window: negative ions as
            # JCAMP-DX, which carries the polarity (a text file opened them as positive ions)
            self.tab.open_full(d["e"])
            return
        folder, base = D._unidec_place(self.frame.out_dir(), d["name"])  # Windows path limit
        os.makedirs(folder, exist_ok=True)
        path = os.path.join(folder, base + ".txt")
        np.savetxt(path, np.asarray(d["spec"], float), fmt="%.6f\t%.4f")
        U.open_in_unidec(path)
        self.status("Opening in the Deconvolute window: " + path)

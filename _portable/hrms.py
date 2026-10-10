"""
HRMS Analysis (MSpektra): high resolution LC-MS / direct infusion data from
Bruker QTOF instruments (maXis, impact, compact, micrOTOF; .d folders with
analysis.baf, or analysis.tsf from otofControl 6 and later), timsTOF .d
folders (analysis.tdf, ion mobility summed), Agilent MassHunter .D folders,
Waters .raw folders, Thermo .raw files and mzML files.

One view holds everything for the mass spectra of a run:
  chromatograms (TIC, BPC, mass chromatograms with their own windows),
  spectra on the native profile m/z axis, peak integration,
  calibration on the calibrant segment of the same run (sodium formate by
  default): right click the spectrum > Calibrate on this spectrum; linear,
  quadratic, cubic, TOF or HPC model; mass errors before and after,
  deconvolution: Bayesian deconvolution, maximum entropy or IsoDec,
  exact mass: ion m/z and isotope pattern of a formula, formula finder.

Reading: hrms_data.py (Bruker Baf2Sql library for analysis.baf, Bruker TDF SDK
for analysis.tsf and analysis.tdf, or mzML) and hrms_vendor.py (Agilent, Waters,
Thermo).
"""
import os
import json

import numpy as np
import wx

import unidec_theme as T
from unidec_theme import C
import unilcms as U
import hrms_data
import hrms_vendor
import hrms_calib as K
import ms_formula as F
from deconv_tab import TableCard

TITLE = "HRMS Analysis"
U.TOOL_ICONS.setdefault("polymer", '<path d="M3 12l4-6 5 6 5-6 4 6M3 12l4 6 5-6 5 6 4-6"/>')
T.ICONS.setdefault("polymer", U.TOOL_ICONS["polymer"])
U.TOOL_ICONS.setdefault("formula", '<path d="M9 3h6M10 3v6l-5.5 9.5A1.7 1.7 0 0 0 6 21h12a1.7 1.7 0 0 0 '
                                   '1.5-2.5L14 9V3"/><path d="M7.5 15h9"/>')
T.ICONS.setdefault("formula", U.TOOL_ICONS["formula"])
ISO_COL = "#D2456F"
CAL_COL = "#0BA064"


def _fnum(v, default):
    try:
        return float(str(v).replace(",", "."))
    except (TypeError, ValueError):
        return default


def _peak_pos(spec, lo, hi, profile=True):
    """Accurate position and height of the tallest peak between lo and hi
    (profile: the apex of the peak shape; centroid data: the tallest stick)."""
    if spec is None or not len(spec):
        return None, 0.0
    return K._peak_position(spec[:, 0], spec[:, 1], lo, hi, profile=profile)


def _pol_word(data, event):
    return "negative" if data.events[event]["polarity"] == "-" else "positive"


def _cals_of(data):
    """Every calibration applied to the file as (calibration, polarities
    "+" and/or "-"), oldest first (one per polarity with polarity switching)."""
    cals = getattr(data, "calibs", None)
    if not cals:
        if data.calib is None:
            return []
        ev = getattr(data, "calib_events", None)
        cals = [(data.calib, "", ev)]
    out = []
    for f, _, ev in cals:
        evs = range(len(data.events)) if ev is None else [k for k in ev if 0 <= k < len(data.events)]
        out.append((f, sorted(set(data.events[k]["polarity"] for k in evs))))
    return out


def _cal_of(data, event):
    """The calibration applied to the spectra of this scan event (a
    calibration made on the ions of one polarity applies to them only)."""
    f = getattr(data, "calibration_for", None)
    return f(event) if f is not None else getattr(data, "calib", None)


def _peak_near(spec, target, d, profile=True):
    """Position and height of the peak nearest to target within +- d (not
    the tallest one there: a neighbour a few ppm away is not the ion looked
    for)."""
    if spec is None or not len(spec):
        return None, 0.0
    return K.peak_near(spec[:, 0], spec[:, 1], target, d, profile=profile)


# ==========================================================================
# mass spectrometry view with calibration, deconvolution and exact mass
# ==========================================================================
class HRMSTab(U.MSTab):
    MZ_FMT = "%.4f"
    PIN_FMT = "%.4f"
    APEX_LABELS = True  # labels at the apex of the profile peak, as DataAnalysis does
    RESULT_TILE_SCALE = 1.0  # deconvolution result rows at full size (the others at half size)
    READ_FMT = "m/z %.4f   %s"
    PICK_MIN = 0.004
    TOL_DEFAULT = "0.02"
    FULL_WINDOW = False  # HRMS data are deconvoluted in the settings window only
    OPEN_HINT = "Open a .d or .raw folder, a .raw file or an mzML file"
    TOOL_HINTS = dict(U.MSTab.TOOL_HINTS)
    TOOL_HINTS["formula"] = "Formula: click a peak to find formulas"
    TOOL_HINTS["calibrate"] = "Calibrate: calibration window for the spectrum shown"
    # m/z of the labels at the apex of a profile peak: the fit of the formula check,
    # the clicked peaks and the calibrant search, so that every number of a peak agrees
    # (used by unilcms.MSTab.plot_spec where it supports it; the report uses it too)
    apex_fn = staticmethod(K.profile_apex)

    def __init__(self, parent, frame):
        U.MSTab.__init__(self, parent, frame)
        from ms_hrms_isotopes import DetailsCard
        self.hrms_details = DetailsCard(self)
        self._layout_rows()
        # tiles at half size by default (the view scrolls); right click > Tile
        # heights also offers: chromatogram and spectrum fill the window
        self.rows.fit_first = 2
        self.rows.mode_key, self.rows.default_mode = "hrms_tile_mode", "journal"
        self.rows.relayout()
        # bin width is not used: spectra keep the native profile m/z axis
        row = self.binw.GetContainingSizer()
        if row is not None:
            row.ShowItems(False)
        self.side.Layout()
        self.side.FitInside()

    def spectrum_tools(self):
        return U.MSTab.spectrum_tools(self)[:3] + [
            ("mode", "formula", "Formula", "formula", self.TOOL_HINTS["formula"]),
            ("sep", ""), ("action", "calibrate", "Calibrate\u2026", "calibrate", self.TOOL_HINTS["calibrate"])] + \
            U.MSTab.spectrum_tools(self)[4:] + [
            ("action", "polymer", "Polymer analysis ▾", "polymer", "Mass distribution, repeat units, chain lengths")]

    def tool_action(self, key):
        if key == "polymer":
            from ms_polymer_ui import menu_items
            button = self.tools.buttons["polymer"]
            menu = wx.Menu()
            actions = {}
            for label, action in menu_items(self):
                item = menu.Append(wx.ID_ANY, label)
                actions[item.GetId()] = action
            try:
                wx.ToolTip.Enable(False)
                selected = button.GetPopupMenuSelectionFromUser(menu, wx.Point(0, button.GetSize().height))
            finally:
                menu.Destroy()
                wx.ToolTip.Enable(True)
            if selected in actions:
                # Open after native menu dismissal, with no accumulating menu bindings.
                wx.CallAfter(actions[selected])
            return
        if key == "calibrate":
            self.cal.open_dialog()
            return
        U.MSTab.tool_action(self, key)

    def extra_sections(self, sp):
        self.cal = CalibPanel(self)
        self.panels.append(self.cal)
        U.MSTab.extra_sections(self, sp)
        sp.section("Exact mass", T.GROUP["red"])
        self.formula = sp.text("", 150, tooltip="Neutral molecule, e.g. C8H10N4O2")
        sp.row("Formula", self.formula)
        self.adduct = sp.choice(F.ADDUCT_LABELS)
        self.adduct.SetMinSize(wx.Size(sp.FromDIP(170), -1))
        sp.row("Ion", self.adduct)
        self.ppm = sp.text("5", 58, tooltip="For matching peaks and the formula finder")
        sp.row("Tolerance", self.ppm, "ppm")
        sp.buttons(U.flat(sp, "Check", "primary", icon="replot", handler=lambda e: self.check_formula(),
                          tooltip="m/z, isotope pattern and error of this ion"),
                   U.flat(sp, "Clear", handler=lambda e: self.clear_overlay()))
        self.exact_note = sp.note("")
        self.formula.Bind(wx.EVT_TEXT_ENTER, lambda e: self.check_formula())

    def set_data(self, data):
        U.MSTab.set_data(self, data)
        self.cal.data_opened(data)

    def chrom_extras(self, ax, scale_of, v):
        out = U.MSTab.chrom_extras(self, ax, scale_of, v)
        cal = getattr(self, "cal", None)
        if cal is None or getattr(v, "kind", "main") != "main" or self.data is None:
            return out
        # calibrant scans: those the calibration of this event was made on, else the
        # calibrant segment found in the file (on the chromatogram of its polarity)
        used = None
        if cal.cal is not None and cal.range and _cal_of(self.data, v.key) is not None and \
                v.key in cal.events_of(cal.event):  # (the other polarity may have a calibration of its own)
            used = cal.range
        rng = used or (cal.found if cal.found_event in (None, v.key) else None)
        if rng:
            t0, t1 = rng
            if t1 - t0 < 0.01:
                t0, t1 = t0 - 0.005, t1 + 0.005
            ax.axvspan(t0, t1, color=CAL_COL, alpha=0.18 if used else 0.10, lw=0, zorder=0)._range_mark = True
            out = list(out) + [(CAL_COL, "calibrant (used)" if used else "calibrant")]
        return out

    def select_range(self, a, b):
        """Average a time range, as if dragged on the chromatogram."""
        self.avg = (a, b) if b > a else (a, a)
        self.pick_t = None
        self.fill_range_fields()
        self.plot_chroms(keep_view=True)
        self.compute_average()

    def extra_chrom_menu(self, x, y, v):
        items = U.MSTab.extra_chrom_menu(self, x, y, v)
        cal = self.cal
        if not self.data:
            return items
        extra = [(None, None)]
        if cal.found:
            a, b = cal.found
            extra += [("Calibrate on the calibrant (%.2f to %.2f min)\u2026" % (a, b),
                       lambda: (self.select_range(a, b), self.cal.open_dialog(rng=(a, b), event=cal.found_event))),
                      ("Show the calibrant spectrum (%.2f to %.2f min)" % (a, b), lambda: self.select_range(a, b))]
        if self.avg and self.avg[1] > self.avg[0]:
            r = self.avg
            # the scans of the chromatogram clicked (with polarity switching: its polarity)
            extra.append(("Calibrate on the selected range (%.2f to %.2f min)\u2026" % r,
                          lambda: self.cal.open_dialog(rng=r, event=v.key)))
        if self.data.calib is not None:
            extra += [("Calibration details\u2026", self.cal.open_details),
                      ("Remove the calibration", self.cal.remove)]
        return items + extra

    def _build_cols(self, n):
        U.MSTab._build_cols(self, n)
        for col in self.cols:
            col["spec_card"].tools_supported = {"xic", "measure", "label", "formula"}

    def _nearest_peak(self, col, x):
        spec = col["spec"]
        if spec is None or not len(spec):
            return None
        x0, x1 = col["spec_card"].ax.get_xlim()
        w = max((x1 - x0) / 150.0, self.PICK_MIN, x * 5e-6)
        m, h = _peak_pos(spec, x - w, x + w, self._profile())
        return m

    def _profile(self):
        return bool(getattr(self.data, "has_profile", True))

    def _ppm(self):
        return U._num(self.ppm, 5.0) or 5.0

    # ---------------------------------------------- reading errors, exports
    def _data_error(self, what, ex):
        """A spectrum that cannot be read (file moved or changed while open,
        damaged scan): said in the window instead of only in the log."""
        import traceback
        traceback.print_exc()
        msg = "%s: %s" % (what, ex)
        self.status(msg)
        wx.MessageBox(msg, TITLE, wx.ICON_WARNING)

    def compute_average(self):
        try:
            U.MSTab.compute_average(self)
        except Exception as ex:
            self._data_error("The spectrum could not be averaged", ex)

    def chrom_pick(self, x, y, view):
        try:
            U.MSTab.chrom_pick(self, x, y, view)
        except Exception as ex:
            self._data_error("The spectrum could not be read", ex)

    def copy_spec(self, e):
        """Spectrum to the clipboard: m/z to 6 decimals (4 decimals were up to
        0.5 ppm off at m/z 100)."""
        spec = self.cols[e]["spec"]
        if spec is None:
            return
        U._clip("\n".join("%.6f\t%.10g" % (a, b) for a, b in spec))
        self.status("Spectrum copied (%d points)" % len(spec))

    def on_export_spectrum(self, e=None):
        """Text or JCAMP-DX file of the spectrum shown (calibrated m/z, as in
        the window), with the m/z to 6 decimals."""
        if not self.data:
            self.status("Open a data file first")
            return
        e = self._chosen(e)
        col = self.cols[e]
        if col["spec"] is None:
            self.status("No spectrum to export")
            return
        path = self.frame.ask_save("Export mass spectrum %s" % self.ev_short(e), self._spec_name(e) + ".txt",
                                   "Text, m/z and intensity (*.txt)|*.txt|JCAMP-DX (*.jdx)|*.jdx")
        if not path:
            return
        try:
            if path.lower().endswith((".jdx", ".dx")):
                t0, t1 = col["range"] if col["range"] else (None, None)
                hrms_data.write_spectrum_jdx(path, col["spec"], self._spec_name(e), self.pol(e) or "+", t0, t1)
            else:
                hrms_data.write_spectrum_txt(path, col["spec"])
        except OSError as ex:
            wx.MessageBox("The spectrum could not be saved:\n%s" % ex, TITLE, wx.ICON_ERROR)
            return
        self.status("Saved " + path)

    def open_full(self, e=None):
        """The spectrum in the full Deconvolute window, with the m/z to 6
        decimals (negative ions as JCAMP-DX, which carries the polarity)."""
        if not self.data:
            return
        e = self._chosen(e)
        col = self.cols[e]
        if col["spec"] is None:
            return
        import ms_deconv
        folder, base = ms_deconv._unidec_place(self.frame.out_dir(), self._spec_name(e))
        os.makedirs(folder, exist_ok=True)
        if self.data.adduct_sign(e) < 0:
            path = os.path.join(folder, base + ".jdx")
            rng = col.get("range") or (None, None)
            hrms_data.write_spectrum_jdx(path, col["spec"], self._spec_name(e), "-", rng[0], rng[1])
        else:
            path = os.path.join(folder, base + ".txt")
            hrms_data.write_spectrum_txt(path, col["spec"])
        U.open_in_unidec(path)
        self.status("Opening in the Deconvolute window: " + path)

    # ------------------------------------------------------------ overlay
    def spectrum_replaced(self, col):
        had_overlay = col.get("overlay") is not None
        U.MSTab.spectrum_replaced(self, col)
        if had_overlay:
            self.exact_note.SetLabel("")
            self._rewrap()

    def plot_spec(self, col, keep_view=False):
        U.MSTab.plot_spec(self, col, keep_view)
        ov = col.get("overlay")
        if not ov or col["spec"] is None:
            return
        ax = col["spec_card"].ax
        pat, scale, text = ov["pattern"], ov["scale"], ov["text"]
        mz = pat[:, 0]
        y = pat[:, 1] / 100.0 * scale
        ax.vlines(mz, 0, y, color=ISO_COL, lw=1.6, alpha=0.55, zorder=4)
        pts, = ax.plot(mz, y, ls="none", marker="_", ms=9, mew=1.4, color=ISO_COL, zorder=5)
        pts._no_scale = True
        ax.annotate(text, xy=(0.99, 0.97), xycoords="axes fraction", ha="right", va="top", fontsize=7.5,
                    color=ISO_COL, zorder=7, bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=ISO_COL, lw=0.6,
                                                       alpha=0.95))
        col["spec_card"].draw()

    def clear_overlay(self):
        for col in self.cols:
            if col.pop("overlay", None) is not None:
                self.plot_spec(col, keep_view=True)
        self.exact_note.SetLabel("")
        self._rewrap()

    def _rewrap(self):
        self.exact_note.Wrap(self.FromDIP(270))
        self.side.Layout()
        self.side.FitInside()

    def check_formula(self, formula=None, adduct=None, col=None):
        if not self.data:
            self.status("Open a data file first")
            return
        try:
            f = F.parse_formula(formula or self.formula.GetValue())
        except ValueError as ex:
            self.status(str(ex))
            wx.MessageBox(str(ex), "Exact mass", wx.ICON_WARNING)
            return
        k = self.adduct.GetSelection() if adduct is None else adduct
        try:
            mz = F.ion_mz(f, k)
            pat = F.isotope_pattern(f, k)
        except ValueError as ex:
            wx.MessageBox(str(ex), "Exact mass", wx.ICON_WARNING)
            return
        z = F.ADDUCTS[k][4]
        name = "%s %s" % (F.format_formula(f), F.ADDUCTS[k][0].split(" ")[0])
        if col is not None and z != 0 and (z > 0) != (self.data.adduct_sign(col["e"]) > 0):
            txt = ("%s is a %s ion, but this spectrum has %s ions: choose an ion of that polarity (e.g. %s)."
                   % (F.ADDUCTS[k][0].split(" ")[0], "positive" if z > 0 else "negative",
                      "positive" if z < 0 else "negative", "[M+H]+" if z < 0 else "[M-H]-"))
            wx.MessageBox(txt, "Exact mass", wx.ICON_INFORMATION)
            self.status(txt)
            return
        cols = [col] if col is not None else [c for c in self.cols if c["spec"] is not None and (
            z == 0 or (z > 0) == (self.data.adduct_sign(c["e"]) > 0))]
        if not cols:
            txt = "%s: m/z %.5f (no %s spectrum to compare)" % (name, mz, "positive" if z > 0 else "negative")
            self.exact_note.SetLabel(txt)
            self._rewrap()
            self.status(txt)
            return
        tol = self._ppm()
        lines = []
        for c in cols:
            spec = c["spec"]
            d = mz * max(tol, 2.0) * 3e-6
            found, h = _peak_near(spec, mz, d, self._profile())
            if found is not None and h > 0:
                err = (found - mz) / mz * 1e6
                meas = F.measured_pattern(spec, found, z or 1, n=len(pat) + 1, ppm=max(tol, 5))
                match = F.isotope_match(pat, meas)
                ok = abs(err) <= tol
                text = "%s\nm/z %.5f (calc.)\nfound %.5f, %+.2f ppm%s\nisotope match %.0f\u00a0%%" % (
                    name, mz, found, err, "" if ok else " (outside tolerance)", match)
                scale = h / (pat[0, 1] / 100.0) if pat[0, 1] > 0 else h
                # scale to the tallest calculated isotope peak
                imax = int(np.argmax(pat[:, 1]))
                mi, hi_ = _peak_near(spec, pat[imax, 0], d, self._profile())
                if hi_ > 0:
                    scale = hi_
                lines.append("%s: found %.5f (%+.2f ppm), isotope match %.0f\u00a0%%" % (
                    self.ev_short(c["e"]), found, err, match))
            else:
                text = "%s\nm/z %.5f (calc.)\nnot found (±%.0f ppm)" % (name, mz, max(tol, 2.0) * 3)
                x0, x1 = mz - 1, mz + 5
                w = (spec[:, 0] >= x0) & (spec[:, 0] <= x1)
                scale = float(spec[w, 1].max()) if np.any(w) else 1.0
                lines.append("%s: not found" % self.ev_short(c["e"]))
            ok_found = found is not None and h > 0
            info = {"name": name, "formula": F.format_formula(f), "ion": F.ADDUCTS[k][0].split(" ")[0],
                    "z": z, "mz": mz, "found": found if ok_found else None,
                    "err": err if ok_found else None, "match": match if ok_found else None,
                    "ion_formula": F.format_formula(F.ion_formula(f, k)[0])}
            if ok_found:  # a peak found outside the tolerance (searched up to 3x) is marked so in the reports
                info["ok"] = bool(abs(err) <= tol)
            c["overlay"] = {"pattern": pat, "scale": scale, "text": text, "info": info}
            card = c["spec_card"]
            span = max(3.5 / max(abs(z), 1), (pat[-1, 0] - pat[0, 0]) + 1.0)
            card.ax.set_xlim(pat[0, 0] - 0.3 * span, pat[0, 0] + span)
            card.ylock = None
            self.plot_spec(c, keep_view=True)
            card.autoscale_y(pad=0.3)
            card.draw()
        txt = "%s  m/z %.5f\n" % (name, mz) + "\n".join(lines)
        self.exact_note.SetLabel(txt)
        self._rewrap()
        self.status(txt.replace("\n", "; "))

    # --------------------------------------------------------- formula tool
    def spec_tool(self, col, tool, x, y):
        if tool == "formula":
            pk = self._nearest_peak(col, x)
            if pk is None:
                self.status("Click closer to a peak")
                return
            self.find_formulas(col, pk)
            return
        U.MSTab.spec_tool(self, col, tool, x, y)

    def spectrum_actions(self, col):
        items = [("Calibrate on this spectrum\u2026", lambda: self.cal.open_dialog(col)),
                 ("Deconvolute this spectrum\u2026", lambda: self.on_deconvolve(col["e"]))]
        if self.data is not None and self.data.calib is not None:
            items.append(("Calibration details\u2026", self.cal.open_details))
        return items

    def spec_menu(self, x, y, col):
        items = U.MSTab.spec_menu(self, x, y, col)
        pk = self._nearest_peak(col, x)
        extra = []
        if pk is not None and self.data:
            extra.append(("Find formulas for m/z %.4f\u2026" % pk, lambda m=pk: self.find_formulas(col, m)))
        if col["spec"] is not None and self.data:
            extra.append(("Check a formula\u2026", lambda: self.ask_formula(col)))
        if col.get("overlay"):
            extra.append(("Remove the isotope pattern", self.clear_overlay))
        if extra:
            k = next((i for i, it in enumerate(items) if it[0] is None), len(items))
            items = items[:k] + [(None, None)] + extra + items[k:]
        return items

    def ask_formula(self, col):
        dlg = FormulaCheckDialog(self, col)
        try:
            if dlg.ShowModal() != wx.ID_OK:
                return
            f, k, tol = dlg.values()
        finally:
            dlg.Destroy()
        self.formula.SetValue(f)
        self.adduct.SetSelection(k)
        self.ppm.SetValue(tol)
        self.check_formula(f, k, col=col)

    # ------------------------------------------- method presets (method_presets.py)
    finder_limits = None  # {element: (min, max)} of the formula finder (last used in this file), None: defaults
    finder_iso = True

    def exact_get(self):
        return {"ion": F.ADDUCT_LABELS[max(0, self.adduct.GetSelection())], "tolerance_ppm": U._num(self.ppm, None)}

    def exact_set(self, v):
        if "ion" in v:
            self.adduct.SetSelection(F.ADDUCT_LABELS.index(v["ion"]))
        if "tolerance_ppm" in v:
            self.ppm.ChangeValue("%g" % v["tolerance_ppm"])

    def finder_get(self):
        lims = dict((e, [lo, hi]) for e, lo, hi in FormulaDialog.ELEMENTS)
        lims.update((e, list(p)) for e, p in (self.finder_limits or {}).items())
        return {"limits": lims, "rank_by_isotopes": bool(self.finder_iso)}

    def finder_set(self, v):
        if "limits" in v:
            lims = dict(self.finder_limits or {})
            lims.update((e, tuple(p)) for e, p in v["limits"].items())
            self.finder_limits = lims
        if "rank_by_isotopes" in v:
            self.finder_iso = bool(v["rank_by_isotopes"])

    def find_formulas(self, col, mz):
        dlg = FormulaDialog(self, col, mz)
        dlg.Show()
        # the dialog is modeless: it goes with the file (closed with Ctrl+W)
        self.Bind(wx.EVT_WINDOW_DESTROY, lambda e, d=dlg: (e.Skip(), _destroy_quietly(d)))


class FormulaCheckDialog(wx.Dialog):
    """Formula, ion and tolerance for the exact mass check."""

    def __init__(self, tab, col=None):
        wx.Dialog.__init__(self, wx.GetTopLevelParent(tab), title="Check a formula")
        self.SetBackgroundColour(wx.Colour(C["panel"]))
        fm = U.Form(self, label_w=90, note_w=300)
        fm.section("Exact mass", T.GROUP["red"])
        self.f = fm.text(tab.formula.GetValue(), 200, tooltip="Neutral molecule, e.g. C8H10N4O2")
        fm.row("Formula", self.f)
        self.ad = fm.choice(F.ADDUCT_LABELS, 200)
        k = tab.adduct.GetSelection()
        if col is not None and tab.data is not None:
            # the ion must have the polarity of the spectrum: [M-H]- for negative ions
            neg = tab.data.adduct_sign(col["e"]) < 0
            if (F.ADDUCTS[k][4] < 0) != neg:
                k = F.ADDUCT_LABELS.index("[M-H]-" if neg else "[M+H]+")
        self.ad.SetSelection(k)
        fm.row("Ion", self.ad)
        self.tol = fm.text(tab.ppm.GetValue(), 60)
        fm.row("Tolerance", self.tol, "ppm")
        vs = wx.BoxSizer(wx.VERTICAL)
        vs.Add(fm, 1, wx.EXPAND)
        vs.Add(U.dialog_buttons(self, U.flat(self, "Cancel", handler=lambda e: self.EndModal(wx.ID_CANCEL)),
                                U.flat(self, "Check", "primary", handler=lambda e: self.EndModal(wx.ID_OK))),
               0, wx.EXPAND)
        self.SetSizerAndFit(vs)
        self.f.Bind(wx.EVT_TEXT_ENTER, lambda e: self.EndModal(wx.ID_OK))
        self.CentreOnParent()
        self.f.SetFocus()

    def values(self):
        return self.f.GetValue().strip(), self.ad.GetSelection(), self.tol.GetValue().strip() or "5"


def _destroy_quietly(win):
    try:
        if win:
            win.Destroy()
    except RuntimeError:
        pass


class FormulaDialog(wx.Dialog):
    """Formula finder for one peak: element limits, ion type, results."""

    ELEMENTS = [("C", 0, 80), ("H", 0, 160), ("N", 0, 10), ("O", 0, 20), ("S", 0, 3), ("P", 0, 2), ("F", 0, 0),
                ("Cl", 0, 0), ("Br", 0, 0), ("Na", 0, 0), ("Si", 0, 0), ("B", 0, 0)]

    def __init__(self, tab, col, mz):
        wx.Dialog.__init__(self, wx.GetTopLevelParent(tab), title="Formula finder",
                           style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        self.tab, self.col, self.mz = tab, col, mz
        self.SetBackgroundColour(wx.Colour(C["panel"]))
        self.results = []
        pos = tab.data.adduct_sign(col["e"]) > 0
        top = wx.FlexGridSizer(0, 4, self.FromDIP(6), self.FromDIP(8))
        self.mz_t = wx.TextCtrl(self, value="%.5f" % mz, size=wx.Size(self.FromDIP(100), -1))
        self.ad = wx.Choice(self, choices=F.ADDUCT_LABELS)
        default = "[M+H]+" if pos else "[M-H]-"
        k = tab.adduct.GetSelection()  # the ion of Exact mass (or of a method) when it has the polarity of the peak
        default = F.ADDUCT_LABELS[k] if k >= 0 and F.ADDUCTS[k][4] and (F.ADDUCTS[k][4] > 0) == pos else default
        self.ad.SetSelection(F.ADDUCT_LABELS.index(default))
        self.ppm_t = wx.TextCtrl(self, value="%g" % tab._ppm(), size=wx.Size(self.FromDIP(60), -1))
        for lab, ctrl in (("m/z", self.mz_t), ("Ion", self.ad), ("Tolerance (ppm)", self.ppm_t)):
            top.Add(wx.StaticText(self, label=lab), 0, wx.ALIGN_CENTER_VERTICAL)
            top.Add(ctrl, 0, wx.ALIGN_CENTER_VERTICAL)
        top.AddSpacer(1)
        top.AddSpacer(1)
        el = wx.FlexGridSizer(0, 8, self.FromDIP(4), self.FromDIP(6))
        self.lims = {}
        for e, lo, hi in self.ELEMENTS:
            lo, hi = (tab.finder_limits or {}).get(e, (lo, hi))  # the limits used last in this file (or a method)
            el.Add(wx.StaticText(self, label=e), 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT, self.FromDIP(6))
            a = wx.SpinCtrl(self, min=0, max=500, initial=lo, size=wx.Size(self.FromDIP(58), -1))
            b = wx.SpinCtrl(self, min=0, max=500, initial=hi, size=wx.Size(self.FromDIP(58), -1))
            h = wx.BoxSizer(wx.HORIZONTAL)
            h.Add(a)
            h.Add(wx.StaticText(self, label="–"), 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT | wx.RIGHT, 3)
            h.Add(b)
            el.Add(h)
            self.lims[e] = (a, b)
            el.AddSpacer(1)
            el.AddSpacer(1)
        self.iso = wx.CheckBox(self, label="Rank by the isotope pattern")
        self.iso.SetValue(bool(tab.finder_iso))
        self.list = wx.ListCtrl(self, style=wx.LC_REPORT | wx.LC_SINGLE_SEL, size=wx.Size(-1, self.FromDIP(260)))
        for n, (name, w) in enumerate([("Formula (neutral)", 150), ("Ion m/z", 100), ("Error (ppm)", 80),
                                       ("RDB", 50), ("Isotope match %", 104)]):
            self.list.InsertColumn(n, name, wx.LIST_FORMAT_LEFT if n == 0 else wx.LIST_FORMAT_RIGHT,
                                   width=self.FromDIP(w))
        self.info = wx.StaticText(self, label="")
        self.info.SetForegroundColour(wx.Colour(C["muted"]))
        btns = wx.BoxSizer(wx.HORIZONTAL)
        find = U.flat(self, "Find", "primary", handler=lambda e: self.run())
        show = U.flat(self, "Show pattern", handler=lambda e: self.show())
        copy = U.flat(self, "Copy list", handler=lambda e: self.copy())
        close = U.flat(self, "Close", handler=lambda e: self.Destroy())
        for b in (find, show, copy):
            btns.Add(b, 0, wx.RIGHT, self.FromDIP(6))
        btns.AddStretchSpacer(1)
        btns.Add(close)
        vs = wx.BoxSizer(wx.VERTICAL)
        pad = self.FromDIP(12)
        vs.Add(top, 0, wx.ALL, pad)
        vs.Add(wx.StaticText(self, label="Elements (minimum – maximum):"), 0, wx.LEFT | wx.RIGHT, pad)
        vs.Add(el, 0, wx.ALL, pad)
        vs.Add(self.iso, 0, wx.LEFT | wx.RIGHT | wx.BOTTOM, pad)
        vs.Add(self.list, 1, wx.EXPAND | wx.LEFT | wx.RIGHT, pad)
        vs.Add(self.info, 0, wx.ALL, pad)
        vs.Add(btns, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, pad)
        self.SetSizerAndFit(vs)
        self.list.Bind(wx.EVT_LIST_ITEM_ACTIVATED, lambda e: self.show())
        self.CentreOnParent()
        wx.CallAfter(self.run)

    def run(self):
        try:
            mz = float(self.mz_t.GetValue().replace(",", "."))
            ppm = float(self.ppm_t.GetValue().replace(",", "."))
        except ValueError:
            self.info.SetLabel("Check m/z and tolerance")
            return
        limits = {e: (a.GetValue(), max(a.GetValue(), b.GetValue())) for e, (a, b) in self.lims.items()}
        self.tab.finder_limits, self.tab.finder_iso = dict(limits), self.iso.GetValue()  # kept for the file
        k = self.ad.GetSelection()
        z = abs(F.ADDUCTS[k][4]) or 1
        obs = None
        if self.iso.GetValue() and self.col["spec"] is not None:
            obs = F.measured_pattern(self.col["spec"], mz, z, n=4, ppm=max(ppm, 5))
        n_comb = 1
        for e, (lo, hi) in limits.items():
            if e not in ("C", "H"):
                n_comb *= (hi - lo + 1)
        if n_comb > 3e5:
            self.info.SetLabel("Too many element combinations (%d); lower some maxima" % n_comb)
            return
        wx.BeginBusyCursor()
        try:
            self.results = F.find_formulas(mz, k, ppm, limits, observed=obs, max_results=50)
        finally:
            wx.EndBusyCursor()
        self.list.DeleteAllItems()
        for n, c in enumerate(self.results):
            i = self.list.InsertItem(n, c["formula"])
            self.list.SetItem(i, 1, "%.5f" % c["mz"])
            self.list.SetItem(i, 2, "%+.2f" % c["ppm"])
            self.list.SetItem(i, 3, "%g" % c["rdb"])
            self.list.SetItem(i, 4, "%.0f" % c["iso"] if "iso" in c else "")
        self.info.SetLabel("%d formula%s within %g ppm of m/z %.5f (%s)" % (
            len(self.results), "s" if len(self.results) != 1 else "", ppm, mz, F.ADDUCTS[k][0]) if self.results
            else "No formula within %g ppm with these elements" % ppm)
        if self.results:
            self.list.Select(0)
        self.Layout()

    def show(self):
        i = self.list.GetFirstSelected()
        if i < 0 or i >= len(self.results):
            return
        c = self.results[i]
        self.tab.formula.SetValue(c["formula"])
        self.tab.adduct.SetSelection(self.ad.GetSelection())
        self.tab.check_formula(c["formula"], self.ad.GetSelection(), col=self.col)

    def copy(self):
        rows = ["formula\tion m/z\terror ppm\tRDB\tisotope match %"]
        for c in self.results:
            rows.append("%s\t%.5f\t%.2f\t%g\t%s" % (c["formula"], c["mz"], c["ppm"], c["rdb"],
                                                   "%.0f" % c["iso"] if "iso" in c else ""))
        U._clip("\n".join(rows))
        self.info.SetLabel("List copied")


# ==========================================================================
# calibration: right click a spectrum > Calibrate on this spectrum
# ==========================================================================
CAL_DEFAULTS = {"calibrant": "Sodium formate, positive", "model": 5, "order": 0, "tol": "30", "minrel": "0.5",
                "custom": "", "auto": False}
ORDER_CHOICES = ["Automatic (cross-validated)"] + ["%d" % k for k in range(1, K.HPC_MAX_ORDER + 1)]


class CalibPanel(object):
    """Calibration state of the view `tab`: the settings, the calibrant
    segment found in the file, the calibration applied, and a short section
    of the side panel (secondary: everything is in the right click menus and
    in the calibration window).

    The calibrant spectrum is the spectrum the user chose (a range dragged on
    the chromatogram, usually the sodium formate injected at the start of the
    run). The reference ions are found in it, the model is fitted between
    measured and reference m/z, and Calibrate applies the correction to every
    spectrum and mass chromatogram of the file."""

    def __init__(self, tab):
        self.tab = tab
        self.frame = tab.frame
        self.data = None
        self.cfg = self._saved_cfg()
        self.found = None  # (t0, t1) calibrant segment found in the file
        self.found_event = None  # scan event (polarity) it was found in
        self.range = None  # (t0, t1) of the spectrum the calibration was made on
        self.cal = None  # calibration applied
        self.rows = []  # calibrants of the applied calibration
        self.spec = None
        self.event = 0
        sp = self.sp = tab.side
        sp.section("Calibration", T.GROUP["green"])
        sp.buttons(U.flat(sp, "Calibrate…", "primary", icon="calibrate", handler=lambda e: self.open_dialog()))
        sp.buttons(U.flat(sp, "Remove", handler=lambda e: self.remove(), tooltip="Back to the m/z values as recorded"),
                   U.flat(sp, "Load…", icon="open", handler=lambda e: self.load_cal(),
                          tooltip="Apply a calibration saved from another file"))
        self.auto = sp.check("Calibrate automatically on opening", bool(self.cfg["auto"]))
        self.auto.SetToolTip("With the last settings, when the calibrant is found")
        self.note = sp.note("")
        self.auto.Bind(wx.EVT_CHECKBOX, lambda e: self.save_cfg(auto=self.auto.GetValue()))

    # ------------------------------------------------------------ helpers
    @staticmethod
    def _saved_cfg():
        """The calibration settings saved last (by any open file) over the defaults."""
        cfg = dict(CAL_DEFAULTS)
        st = T._load()
        saved = st.get("hrms_cal")
        if isinstance(saved, dict):
            cfg.update({k: v for k, v in saved.items() if k in cfg})
        elif "hrms_autocal" in st:
            cfg["auto"] = bool(st.get("hrms_autocal"))
        return cfg

    def status(self, text):
        self.tab.status(text)

    def save_cfg(self, **kw):
        # over the settings saved last: each open file has its own panel, and the whole
        # settings of a file opened earlier replaced those used since on another file
        self.cfg = self._saved_cfg()
        self.cfg.update(kw)
        T._save({"hrms_cal": dict(self.cfg)})

    def method_get(self):
        """Calibration settings in a method preset (method_presets.py; not the calibration of the file)."""
        c = self._saved_cfg()
        return {"calibrant": c["calibrant"], "model": K.MODELS[int(c["model"])] if
                str(c["model"]).isdigit() and int(c["model"]) < len(K.MODELS) else None,
                "order": int(c.get("order") or 0), "tolerance_ppm": _fnum(c["tol"], None),
                "min_intensity": _fnum(c["minrel"], None), "custom_list": U._numlist(str(c.get("custom") or "")),
                "auto": bool(c["auto"])}

    def method_set(self, v):
        """The values given (the others stay): used by the next calibration (and on opening when Calibrate
        automatically is ticked), as if chosen in the calibration window."""
        kw = {k: v[k] for k in ("calibrant", "order", "auto") if k in v}
        if "model" in v:
            kw["model"] = K.MODELS.index(v["model"])
        for key, ck in (("tolerance_ppm", "tol"), ("min_intensity", "minrel")):
            if key in v:
                kw[ck] = "%g" % v[key]
        if "custom_list" in v:
            kw["custom"] = ", ".join("%.10g" % x for x in v["custom_list"])
        self.save_cfg(**kw)
        self.auto.SetValue(bool(self.cfg["auto"]))

    def _set_note(self, text):
        self.note.SetLabel(text)
        self.note.Wrap(self.sp.FromDIP(270))
        self.sp.Layout()
        self.sp.FitInside()

    def refs(self, cfg=None):
        cfg = cfg or self.cfg
        if cfg["calibrant"] == "Custom list":
            return sorted(set(U._numlist(cfg.get("custom", ""))))
        return K.calibrant_list(cfg["calibrant"])

    def event_for(self, name, data=None):
        """Scan event of the calibrant polarity."""
        data = data or self.data
        if not data:
            return 0
        want = "-" if "negative" in name else "+"
        for e, ev in enumerate(data.events):
            if ev["polarity"] == want:
                return e
        return 0

    def calibrant_for_event(self, e, data=None):
        """Calibrant of the settings, switched to the polarity of event e."""
        data = data or self.data
        pol = data.events[e]["polarity"] if data else "+"
        name = self.cfg["calibrant"]
        names = K.calibrant_names(pol)
        if name in names:
            return name
        base = name.split(",")[0]
        for n in names:
            if n.split(",")[0] == base:
                return n
        return names[0]

    def events_of(self, e):
        """The scan events with the polarity of event e: a calibration is made
        on the ions of one polarity and applies to them only (with polarity
        switching the other polarity has a calibration of its own)."""
        if not self.data:
            return [e]
        pol = self.data.events[e]["polarity"]
        return [k for k, ev in enumerate(self.data.events) if ev["polarity"] == pol]

    def _segment_args(self, data):
        name = self.calibrant_for_event(0, data)
        e = self.event_for(name, data)
        return name, e, tuple(self.refs(dict(self.cfg, calibrant=name))), max(_fnum(self.cfg["tol"], 30.0), 30.0)

    def precompute(self, data):
        """Run in the thread that reads the file (no window calls): looks for
        the calibrant segment, so that the window does not wait for the mass
        chromatograms of every reference ion when the file is shown."""
        name, e, refs, tol = self._segment_args(data)
        data._cal_segment = ((name, e, refs, tol), K.find_calibrant_segment(data, refs, e, tol_ppm=tol))

    def reset(self):
        pass  # results live in the calibration window; nothing drawn in the view

    def overlay(self, col, ax):
        pass

    def raw_spectrum(self, e, t0, t1):
        """Spectrum as recorded (no calibration) of event e, averaged over
        t0 to t1 min, or the single scan nearest to t0 when t1 == t0."""
        cal = self.data.calib
        self.data.calib = None
        try:
            if t1 is None or abs(t1 - t0) < 1e-9:
                s, k = self.data.scan_spectrum(e, t0)
                return s, 1, (float(self.data.rt[k]), float(self.data.rt[k]))
            s, n = self.data.average(e, min(t0, t1), max(t0, t1))
            return s, n, (min(t0, t1), max(t0, t1))
        finally:
            self.data.calib = cal

    def measure(self, spec, cfg):
        return K.find_calibrants(spec, self.refs(cfg), tol_ppm=_fnum(cfg["tol"], 30.0),
                                 min_rel=_fnum(cfg["minrel"], 0.5) / 100.0, profile=self.data.has_profile)

    @staticmethod
    def fit(rows, cfg):
        used = [r for r in rows if r["use"] and r["found"] is not None]
        order = int(cfg.get("order") or 0) or None
        return K.Calibration(cfg["model"], [r["found"] for r in used], [r["ref"] for r in used],
                             hpc_order=order if cfg["model"] == 5 else None)

    # ------------------------------------------------------- file opened
    def data_opened(self, data):
        self.data = data
        self.found, self.range, self.cal, self.rows, self.spec = None, None, None, [], None
        self.found_event = None
        self._set_note("")
        if data is None:
            return
        args = self._segment_args(data)
        name, e, refs, tol = args
        pre = getattr(data, "_cal_segment", None)
        if pre is not None and pre[0] == args:  # found while the file was read
            seg = pre[1]
        else:
            seg = K.find_calibrant_segment(data, refs, e, tol_ppm=tol)
        if seg is None:
            self._set_note("No %s segment found" % name.split(",")[0].lower())
            return
        self.found = (seg[0], seg[1])
        self.found_event = e
        self._set_note("%s · %.2f–%.2f min" %
                       (name.split(",")[0], seg[0], seg[1]))
        if self.auto.GetValue():
            wx.CallAfter(self.auto_calibrate, e, name)

    def auto_calibrate(self, e, name):
        if not self.data or not self.found:
            return
        cfg = dict(self.cfg, calibrant=name)
        spec, n, rng = self.raw_spectrum(e, *self.found)
        rows = self.measure(spec, cfg) if len(spec) else []
        try:
            cal = self.fit(rows, cfg)
        except (ValueError, np.linalg.LinAlgError) as ex:
            self._set_note("Automatic calibration failed: %s" % ex)
            return
        if cal.rms > cal.rms_before:
            self._set_note("Automatic calibration not applied: the %s model makes the errors larger (%.2f ppm, as "
                           "recorded %.2f ppm)." % (cal.name(), cal.rms, cal.rms_before))
            return
        self.apply(cal, rows, spec, rng, e, cfg)
        self.frame.calibration_changed()

    # ------------------------------------------------------------ actions
    def open_details(self):
        """Calibration window of the calibration applied (its spectrum and
        scan event), else as the Calibrate button."""
        if self.range:
            self.open_dialog(rng=self.range, event=self.event)
        else:
            self.open_dialog()

    def open_dialog(self, col=None, rng=None, event=None):
        """Calibration window for the spectrum of col (or the time range rng
        of the scan event given, by default that of the calibrant polarity,
        or the spectrum shown, or the last calibrant spectrum)."""
        tab = self.tab
        if not self.data:
            self.status("Open a data file first")
            return
        e, t0, t1 = None, None, None
        if rng is not None:
            if event is not None and 0 <= int(event) < len(self.data.events):
                e = int(event)
            else:
                e = self.event_for(self.calibrant_for_event(0))
            t0, t1 = rng
        else:
            cols = [col] if col is not None else [c for c in tab.cols if c["spec"] is not None]
            if col is None and cols:
                want = self.event_for(self.calibrant_for_event(0))
                cols.sort(key=lambda c: c["e"] != want)
            if cols and cols[0]["spec"] is not None:
                c = cols[0]
                e = c["e"]
                if c["range"]:
                    t0, t1 = c["range"]
                elif tab.pick_t is not None:
                    t0, t1 = tab.pick_t, tab.pick_t
            if t0 is None and self.range is not None:
                e, (t0, t1) = self.event, self.range
            if t0 is None and self.found is not None:
                e = self.found_event if self.found_event is not None else self.event_for(self.calibrant_for_event(0))
                t0, t1 = self.found
        if t0 is None:
            wx.MessageBox("Choose the calibrant spectrum first: drag across the calibrant peak on the "
                          "chromatogram.", "Calibration", wx.ICON_INFORMATION)
            return
        self.cfg = self._saved_cfg()  # the last settings, also when they were used on another open file
        try:
            self.auto.SetValue(bool(self.cfg["auto"]))
        except RuntimeError:
            pass
        dlg = CalibDialog(self, e, t0, t1)
        try:
            dlg.ShowModal()
            changed = dlg.changed
        finally:
            dlg.Destroy()
        if changed:
            self.frame.calibration_changed()

    def apply(self, cal, rows, spec, rng, e, cfg):
        self.cal, self.rows, self.spec, self.range, self.event = cal, [dict(r) for r in rows], spec, rng, e
        self.data.set_calibration(cal, cal.describe(), events=self.events_of(e))
        self.save_cfg(**{k: v for k, v in cfg.items() if k in CAL_DEFAULTS and k != "auto"})
        self.save_report(cfg["calibrant"])
        self._set_note("%s · %.2f ppm error · %d reference peaks" %
                       (cal.name(), cal.rms, len(cal.measured)))
        self.status("Calibrated: %s" % cal.describe())

    def remove(self):
        if not self.data:
            return
        self.data.set_calibration(None)
        self.cal, self.range = None, None
        self.frame.calibration_changed()
        self._set_note("Calibration removed: m/z values as recorded.")
        self.status("Calibration removed: m/z values as recorded")

    def save_report(self, calibrant):
        try:
            out = self.frame.out_dir()
            base = os.path.join(out, self.frame.base_name() + "_calibration")
            with open(base + ".txt", "w", encoding="utf-8") as fh:
                fh.write(self.cal.report(self.frame.base_name(), calibrant, self.rows))
            with open(base + ".json", "w", encoding="utf-8") as fh:
                json.dump({"model": self.cal.model, "hpc_order": getattr(self.cal, "hpc_order", None),
                           "coef": list(map(float, self.cal.coef)), "lo": self.cal.lo, "hi": self.cal.hi,
                           "rms": self.cal.rms, "cv_rms": self.cal.cv_rms,
                           "measured": list(map(float, self.cal.measured)),
                           "reference": list(map(float, self.cal.reference)),
                           "calibrant": calibrant, "range": self.range, "data": self.frame.base_name(),
                           "polarity": self.data.events[self.event]["polarity"] if self.data else ""}, fh, indent=1)
            return base + ".txt"
        except OSError as ex:
            self.status("The calibration report could not be saved: %s" % ex)
            return None

    def load_cal(self, parent=None):
        if not self.data:
            self.status("Open a data file first")
            return False
        dlg = wx.FileDialog(parent or self.tab, "Load a calibration", defaultDir=self.frame.out_dir(),
                            wildcard="MSpektra calibration (*_calibration.json)|*.json",
                            style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
        try:
            if dlg.ShowModal() != wx.ID_OK:
                return False
            path = dlg.GetPath()
        finally:
            dlg.Destroy()
        try:
            with open(path, "r", encoding="utf-8") as fh:
                d = json.load(fh)
            cal = K.Calibration(d["model"], d["measured"], d["reference"], hpc_order=d.get("hpc_order"))
        except Exception as ex:
            wx.MessageBox("This calibration could not be read:\n%s" % ex, "Calibration", wx.ICON_ERROR)
            return False
        # the ions of the polarity it was made on (files saved before: from the calibrant name)
        pol = d.get("polarity") or ("-" if "negative" in str(d.get("calibrant", "")).lower() else "+")
        events = [k for k, ev in enumerate(self.data.events) if ev["polarity"] == pol]
        where = ""
        if not events:  # the file has only the other polarity: applied to it, and said so
            events = list(range(len(self.data.events)))
            where = " (made on %s ions, applied to the %s ions of this file)" % (
                "negative" if pol == "-" else "positive", "positive" if pol == "-" else "negative")
        self.cal, self.rows, self.spec, self.range = cal, [], None, None
        self.event = events[0]
        self.data.set_calibration(cal, cal.describe() + " (from %s)" % d.get("data", os.path.basename(path)),
                                  events=events)
        if parent is None:
            self.frame.calibration_changed()
        self._set_note("Applied the calibration of %s: %s%s" % (d.get("data", "?"), cal.describe(), where))
        self.status("Calibration loaded from " + path + where)
        return True


class CalibDialog(wx.Dialog):
    """Calibration on one spectrum: settings on the left; the calibrant
    spectrum, the mass errors, the models compared and the calibrant table
    on the right. The fit is shown as soon as the window opens; Calibrate
    applies it to the whole file. Close returns to the main window, which
    then shows the calibrated spectra."""

    def __init__(self, cp, e, t0, t1):
        frame = wx.GetTopLevelParent(cp.tab)
        wx.Dialog.__init__(self, frame, title="Calibration", style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER |
                           wx.MAXIMIZE_BOX)
        self.cp, self.data, self.e = cp, cp.data, e
        self.changed = False  # calibration of the file changed (applied or removed)
        self.applied_sig = None
        self.rows, self.spec, self.cal, self.rng, self.n = [], None, None, (t0, t1), 0
        self.SetBackgroundColour(wx.Colour(C["bg"]))
        cfg = dict(cp.cfg, calibrant=cp.calibrant_for_event(e))

        # ---------------------------------------------------------- settings
        lw = self.leftwin = wx.ScrolledWindow(self, style=wx.VSCROLL | wx.BORDER_NONE)
        lw.SetBackgroundColour(wx.Colour(C["panel"]))
        lw.SetScrollRate(0, self.FromDIP(12))
        fm = self.form = U.Form(lw, label_w=118, note_w=292)
        fm.section("Calibrant spectrum", T.GROUP["blue"])
        self.t0, self.t1 = fm.text("%.3f" % t0, 64), fm.text("%.3f" % t1, 64)
        # the range as given (a drag on the chromatogram, the calibrant segment found): the
        # fields show it to 3 decimals, which would leave out the scans at its ends
        self._given = ((self.t0.GetValue(), self.t1.GetValue()), (t0, t1))
        fm.row("Time", self.t0, "to", self.t1, "min")
        self.spec_note = fm.note("")
        fm.section("Calibrant", T.GROUP["green"])
        pol = self.data.events[e]["polarity"] or "+"
        self.calib = fm.choice(K.calibrant_names(pol), 292)
        self.calib.SetStringSelection(cfg["calibrant"])
        fm.full(self.calib)
        self.custom = fm.text(cfg.get("custom", ""), 292, tooltip="Reference m/z values separated by commas")
        fm.full(self.custom)
        self.tol = fm.text(str(cfg["tol"]), 60)
        fm.row("Search window \u00b1", self.tol, "ppm")
        self.minrel = fm.text(str(cfg["minrel"]), 60)
        fm.row("Minimum intensity", self.minrel, "% of the base peak")
        fm.buttons(U.flat(fm, "Reference list from a file\u2026", icon="open", height=28,
                          handler=lambda e_: self.load_list(),
                          tooltip="One m/z per line"))
        fm.section("Model", T.GROUP["yellow"])
        self.model = fm.choice(K.MODELS, 292)
        self.model.SetSelection(int(cfg["model"]))
        fm.full(self.model)
        self.order = fm.choice(ORDER_CHOICES, 174, tooltip="Automatic: chosen by leave one out prediction")
        self.order.SetSelection(int(cfg.get("order") or 0))
        self.row_order = fm.row("HPC order", self.order)
        self.go = U.flat(fm, "Calibrate", "primary", icon="calibrate", height=36, handler=lambda e_: self.calibrate(),
                         tooltip="Apply to the whole file")
        fm.buttons(self.go)
        self.result = fm.note("")
        self.result.SetFont(T.ui_font(9, 600))
        self.auto = fm.check("Calibrate automatically when a file is opened", bool(cp.auto.GetValue()))
        ls = wx.BoxSizer(wx.VERTICAL)
        ls.Add(fm, 1, wx.EXPAND)
        lw.SetSizer(ls)
        sb = wx.SystemSettings.GetMetric(wx.SYS_VSCROLL_X, self)
        lw.SetMinSize(wx.Size(fm.GetBestSize().width + max(sb, self.FromDIP(16)), -1))

        # ----------------------------------------------------------- results
        right = wx.Panel(self, style=wx.BORDER_NONE)
        right.SetBackgroundColour(wx.Colour(C["bg"]))
        self.scroller = U.RowScroller(right)
        self.split = U.StackBox(self.scroller)
        self.split.shape = U.tile_shape
        # spectrum and mass error fill the window; the calibrant table below (scroll)
        self.split.fit_first = 2
        self.split.mode_key, self.split.default_mode = "calib_tile_mode", "two"
        self.card_spec = U.PlotCard(self.split, "Calibrant spectrum", "", mode="zoom")
        mid = U.SplitBox(self.split, wx.HORIZONTAL, min_size=200)
        self.card_err = U.PlotCard(mid, "Mass error", "", mode="zoom")
        self.models = TableCard(mid, "Models", [("Model", 116), ("RMS (ppm)", 76), ("CV (ppm)", 70)],
                                on_activate=self.pick_model, empty="double click a model to use it")
        mid.set_panes([self.card_err, self.models], [2.4, 1])
        self.table = TableCard(self.split, "Calibrants", [
            ("Use", 44), ("Reference m/z", 110), ("Found m/z", 110), ("Intensity", 90), ("As recorded (ppm)", 118),
            ("Calibrated (ppm)", 110), ("Left out (ppm)", 100)], on_activate=self.toggle,
            empty="no calibrants", multi=True)
        self.table.extra_menu = self.table_menu
        self.table.on_delete = lambda rows: self.set_use(rows, False)
        # the calibrants selected (in a plot or the table) are ringed in both plots
        self.table.on_select = lambda i: self._later_highlight()
        self.table.list.Bind(wx.EVT_LIST_ITEM_DESELECTED, lambda e: self._later_highlight())
        self.split.set_panes([self.card_spec, mid, self.table], [1.25, 1.35, 1.0])
        self.scroller.set_box(self.split)
        self.card_spec.on_view = lambda: self.plot_spec(keep=True)
        self.card_spec.readout_fmt = lambda x, y: "m/z %.4f" % x
        self.card_err.readout_fmt = lambda x, y: "m/z %.1f   %.2f ppm" % (x, y)
        self.card_spec.image_name = lambda: frame.base_name() + "_calibrant_spectrum"
        self.card_err.image_name = lambda: frame.base_name() + "_calibration_error"
        self.card_spec.on_menu = lambda x, y: self.point_menu(x, None, self.card_spec)
        self.card_err.on_menu = lambda x, y: self.point_menu(x, y, self.card_err)
        self.card_spec.on_pick = lambda x, y: self.pick_point(x, None)
        self.card_err.on_pick = lambda x, y: self.pick_point(x, y)
        rs = wx.BoxSizer(wx.VERTICAL)
        rs.Add(self.scroller, 1, wx.EXPAND | wx.ALL, self.FromDIP(12))
        right.SetSizer(rs)

        # ----------------------------------------------------------- buttons
        bs = wx.BoxSizer(wx.HORIZONTAL)
        bs.Add(U.flat(self, "Remove calibration", handler=lambda e_: self.remove(),
                      tooltip="Back to the m/z values as recorded"), 0, wx.RIGHT, self.FromDIP(8))
        bs.Add(U.flat(self, "Load calibration…", icon="open", handler=lambda e_: self.load(),
                      tooltip="Apply a calibration saved from another file"), 0, wx.RIGHT, self.FromDIP(8))
        self.status_line = wx.StaticText(self, label="")
        self.status_line.SetFont(T.ui_font(9, 400))
        self.status_line.SetForegroundColour(wx.Colour(C["muted"]))
        bs.Add(self.status_line, 1, wx.ALIGN_CENTER_VERTICAL | wx.LEFT, self.FromDIP(8))
        bs.Add(U.flat(self, "Close", "primary", height=34, min_width=110, handler=lambda e_: self.close()), 0)
        body = wx.BoxSizer(wx.HORIZONTAL)
        body.Add(lw, 0, wx.EXPAND)
        body.Add(right, 1, wx.EXPAND)
        vs = wx.BoxSizer(wx.VERTICAL)
        vs.Add(body, 1, wx.EXPAND)
        line = wx.Panel(self, size=wx.Size(-1, 1))
        line.SetBackgroundColour(wx.Colour(C["line"]))
        vs.Add(line, 0, wx.EXPAND)
        vs.Add(bs, 0, wx.EXPAND | wx.ALL, self.FromDIP(12))
        self.SetSizer(vs)
        area = T.display_area(self)  # that of the main window
        self.SetSize(wx.Size(int(area.width * 0.96), int(area.height * 0.94)))  # large: the plots get the room
        self.SetMinSize(wx.Size(min(self.FromDIP(900), area.width), min(self.FromDIP(600), area.height)))
        self.CentreOnParent()  # (inside the screen when it is shown: T.fit_to_screen)

        self.calib.Bind(wx.EVT_CHOICE, lambda e_: self._calib_changed(remeasure=True))
        self.model.Bind(wx.EVT_CHOICE, lambda e_: self.refit())
        self.order.Bind(wx.EVT_CHOICE, lambda e_: self.refit())
        for c in (self.t0, self.t1):
            c.Bind(wx.EVT_TEXT_ENTER, lambda e_: self.load_spectrum())
        for c in (self.tol, self.minrel, self.custom):
            c.Bind(wx.EVT_TEXT_ENTER, lambda e_: self.remeasure())
        self.auto.Bind(wx.EVT_CHECKBOX, lambda e_: self._auto())
        self.Bind(wx.EVT_CLOSE, lambda e_: self.close())
        self.Bind(wx.EVT_CHAR_HOOK, self._keys)
        self._calib_changed(remeasure=False)
        wx.CallAfter(self.load_spectrum)

    # ------------------------------------------------------------ helpers
    def ask_save(self, *a):
        return self.cp.frame.ask_save(*a)

    def base_name(self):
        return self.cp.frame.base_name()

    def SetStatusText(self, text, i=0):
        self.status_line.SetLabel(text)
        self.Layout()

    def _keys(self, e):
        k = e.GetKeyCode()
        if k == wx.WXK_ESCAPE:
            self.close()
        elif k in (wx.WXK_DELETE, wx.WXK_BACK) and not isinstance(wx.Window.FindFocus(), wx.TextCtrl):
            sel = self.table.selected()
            if sel:  # a calibrant picked in a plot or the table: Delete removes it from the fit
                self.set_use(sel, False)
            else:
                e.Skip()
        else:
            e.Skip()

    def _wrap(self, st, text):
        self.form.set_note(st, text)
        self._relayout()

    def _relayout(self):
        self.form.Layout()
        self.leftwin.FitInside()
        self.leftwin.Layout()

    def cfg(self):
        return {"calibrant": self.calib.GetStringSelection(), "model": self.model.GetSelection(),
                "order": self.order.GetSelection(), "tol": self.tol.GetValue().strip() or "30",
                "minrel": self.minrel.GetValue().strip() or "0.5", "custom": self.custom.GetValue().strip()}

    def _sig(self):
        c = self.cfg()
        return (tuple(sorted(c.items())), self.rng, tuple((r["ref"], r["use"]) for r in self.rows))

    def _auto(self):
        self.cp.auto.SetValue(self.auto.GetValue())
        self.cp.save_cfg(auto=self.auto.GetValue())

    def _calib_changed(self, remeasure=True):
        custom = self.calib.GetStringSelection() == "Custom list"
        self.custom.Show(custom)
        if not custom:  # 6 decimals: the list shown becomes the custom list when it is edited
            self.custom.SetValue(", ".join("%.6f" % v for v in K.calibrant_list(self.calib.GetStringSelection())))
        self._relayout()
        if remeasure:
            self.remeasure()

    def load_list(self):
        dlg = wx.FileDialog(self, "Reference m/z list", defaultDir=self.cp.frame.folder(),
                            wildcard="Text or csv (*.txt;*.csv;*.ref)|*.txt;*.csv;*.ref|All files (*.*)|*.*",
                            style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
        try:
            if dlg.ShowModal() != wx.ID_OK:
                return
            vals = K.read_list(dlg.GetPath())
        finally:
            dlg.Destroy()
        if not vals:
            self.SetStatusText("No m/z values found in the file")
            return
        self.calib.SetStringSelection("Custom list")
        self._calib_changed(remeasure=False)
        self.custom.SetValue(", ".join("%.6f" % v for v in vals))
        self.remeasure()

    # --------------------------------------------------------------- work
    def load_spectrum(self):
        t0, t1 = U._num(self.t0), U._num(self.t1)
        shown, exact = self._given
        if (self.t0.GetValue(), self.t1.GetValue()) == shown:  # not edited: the range itself
            t0, t1 = exact
        if t0 is None:
            self.SetStatusText("Enter the time range of the calibrant")
            return
        if t1 is None:
            t1 = t0
        wx.BeginBusyCursor()
        try:
            self.spec, self.n, self.rng = self.cp.raw_spectrum(self.e, t0, t1)
        finally:
            wx.EndBusyCursor()
        if not len(self.spec):
            self._wrap(self.spec_note, "No scans in %.3f to %.3f min." % (t0, t1))
            self.rows = []
            self.refit()
            return
        what = ("scan at %.3f min" % self.rng[0]) if self.n == 1 and self.rng[0] == self.rng[1] else \
            "%.3f to %.3f min, %d scan%s averaged" % (self.rng[0], self.rng[1], self.n, "" if self.n == 1 else "s")
        self._wrap(self.spec_note, "%s, %s, m/z as recorded." % (self.data.event_label(self.e), what))
        self.card_spec.ylock = None
        self.remeasure(keep_view=False)

    def remeasure(self, keep_view=True):
        if self.spec is None or not len(self.spec):
            return
        self.rows = self.cp.measure(self.spec, self.cfg())
        self.refit(keep_view=keep_view)

    def refit(self, keep_view=True):
        self.row_order.ShowItems(self.model.GetSelection() == 5)
        self._relayout()
        self.cal, err = None, ""
        if self.rows:
            try:
                self.cal = self.cp.fit(self.rows, self.cfg())
            except (ValueError, np.linalg.LinAlgError) as ex:
                err = str(ex)
        self._result(err)
        self.plot_spec(keep=keep_view)
        self.plot_errors()
        self.fill_tables()

    def _result(self, err=""):
        c = self.cal
        applied = self.applied_sig is not None and self.applied_sig == self._sig()
        if c is None:
            found = sum(1 for r in self.rows if r["found"] is not None)
            txt = "Not calibrated: %s" % (err or ("%d of %d calibrant peaks found" % (found, len(self.rows))
                                                  if self.rows else "no calibrant peaks in this spectrum"))
            self.go.Enable(False)
        else:
            cv = "%.2f ppm" % c.cv_rms if c.cv_rms is not None else "n/a"
            txt = "%s: %.2f ppm RMS on %d calibrants (cross-validated %s); as recorded %.2f ppm." % (
                c.name(), c.rms, len(c.measured), cv, c.rms_before)
            if c.rms > c.rms_before:
                txt += " This model makes the errors larger."
            txt += " Applied to the file." if applied else " Not applied yet."
            self.go.Enable(True)
        self.result.SetForegroundColour(wx.Colour("#0B7A4B" if applied else C["text"]))
        self._wrap(self.result, txt)

    def calibrate(self):
        if self.cal is None:
            return
        cfg = self.cfg()
        self.cp.apply(self.cal, self.rows, self.spec, self.rng, self.e, cfg)
        self.changed = True
        self.applied_sig = self._sig()
        self._result()
        self.SetStatusText("Calibrated: %s" % self.cal.describe())

    def remove(self):
        if _cal_of(self.data, self.e) is None:
            self.SetStatusText("The file is not calibrated (m/z as recorded)" if self.data.calib is None else
                               "These ions are not calibrated (m/z as recorded)")
            return
        # the calibration of this polarity: one of the other polarity stays
        evs = self.cp.events_of(self.e)
        self.data.set_calibration(None, events=evs)
        if self.data.calib is None or self.cp.event in evs:
            self.cp.cal, self.cp.range = None, None
        left = self.data.calib is not None
        word = _pol_word(self.data, self.e)
        self.cp._set_note("Calibration of the %s ions removed: their m/z values as recorded." % word if left else
                          "Calibration removed: m/z values as recorded.")
        self.changed = True
        self.applied_sig = None
        self._result()
        self.SetStatusText("Calibration of the %s ions removed" % word if left else
                           "Calibration removed: the file has the m/z values as recorded")

    def load(self):
        if self.cp.load_cal(parent=self):
            self.changed = True
            self.applied_sig = None
            self._result()
            self.SetStatusText(self.cp.cal.describe() + " (loaded)")

    def toggle(self, i):
        if 0 <= i < len(self.rows) and self.rows[i]["found"] is not None:
            self.set_use([i], not self.rows[i]["use"])

    def set_use(self, idx, flag):
        """Remove calibrants from the fit (flag False) or use them again."""
        idx = [i for i in idx if 0 <= i < len(self.rows) and self.rows[i]["found"] is not None]
        if not idx:
            return
        for i in idx:
            self.rows[i]["use"] = bool(flag)
        self.refit()
        for n, i in enumerate(idx):
            self.table.select(i, only=(n == 0))
        self.SetStatusText("%d calibrant%s %s" % (len(idx), "s" if len(idx) > 1 else "",
                                                  "used again" if flag else "removed from the fit"))

    def _nearest(self, x, y=None):
        """Row of the calibrant closest to a click (m/z, and error in ppm on
        the mass error plot)."""
        best, dist = None, None
        for i, r in enumerate(self.rows):
            if r["found"] is None:
                continue
            d = abs(r["found"] - x) / max(x, 1.0) * 1e3  # per mil of m/z
            if y is not None:
                errs = [r["err"]]
                if self.cal is not None:
                    errs.append((float(self.cal(np.array([r["found"]]))[0]) - r["ref"]) / r["ref"] * 1e6)
                d += min(abs(e - y) for e in errs) * 0.2
            if dist is None or d < dist:
                best, dist = i, d
        return best

    def _hit(self, card, x, y):
        """Calibrant whose marker is under the pointer (within about 14
        pixels), on the mass error plot or the calibrant spectrum."""
        import math
        ax = card.ax
        try:
            px, py = ax.transData.transform((x, y if y is not None else 0.0))
        except Exception:
            return None
        best, bd = None, None
        for i, r in enumerate(self.rows):
            if r["found"] is None:
                continue
            if y is None:  # spectrum: the triangle at the found m/z
                pts = [(r["found"], r["intensity"])]
            else:
                pts = [(r["found"], r["err"])]
                if self.cal is not None:
                    pts.append((r["found"], (float(self.cal(np.array([r["found"]]))[0]) - r["ref"]) / r["ref"] * 1e6))
            for a, b in pts:
                try:
                    qx, qy = ax.transData.transform((a, b))
                except Exception:
                    continue
                d = abs(qx - px) if y is None else math.hypot(qx - px, qy - py)
                if bd is None or d < bd:
                    best, bd = i, d
        return best if bd is not None and bd <= card.FromDIP(14) else None

    # ------------------------------------------------------- highlighting
    def _later_highlight(self):
        if not getattr(self, "_hl_pending", False):
            self._hl_pending = True
            wx.CallAfter(self._draw_highlight)

    def _hl_rows(self):
        try:
            sel = self.table.selected()
        except Exception:
            return []
        return [i for i in sel if 0 <= i < len(self.rows) and self.rows[i]["found"] is not None]

    def _hl_card(self, card):
        """Orange rings around the selected calibrants on one plot."""
        for a in getattr(card, "_hl_arts", []):
            try:
                a.remove()
            except Exception:
                pass
        card._hl_arts = []
        rows = self._hl_rows()
        if not rows or card.full is None:
            return
        xs, ys = [], []
        for i in rows:
            r = self.rows[i]
            if card is self.card_spec:
                xs.append(r["found"])
                ys.append(r["intensity"])
            else:
                xs.append(r["found"])
                ys.append(r["err"])
                if self.cal is not None:
                    xs.append(r["found"])
                    ys.append((float(self.cal(np.array([r["found"]]))[0]) - r["ref"]) / r["ref"] * 1e6)
        ring, = card.ax.plot(xs, ys, ls="none", marker="o", ms=13, mfc="none", mec="#E0602F", mew=1.8, zorder=8)
        ring._no_scale = True
        card._hl_arts.append(ring)
        if len(rows) == 1:
            r = self.rows[rows[0]]
            ln = card.ax.axvline(r["found"], color="#E0602F", lw=0.8, ls=(0, (3, 2)), alpha=0.8, zorder=1)
            ln._no_scale = True
            card._hl_arts.append(ln)

    def _draw_highlight(self):
        self._hl_pending = False
        try:
            for card in (self.card_spec, self.card_err):
                self._hl_card(card)
                card.canvas.draw_idle()
        except RuntimeError:
            pass  # window closed

    def pick_point(self, x, y):
        card = self.card_err if y is not None else self.card_spec
        i = self._hit(card, x, y)
        if i is None:
            i = self._nearest(x, y)
        if i is not None:
            self.table.select(i)
            self._later_highlight()
            r = self.rows[i]
            self.SetStatusText("m/z %.4f (reference %.5f): %+.2f ppm as recorded; %s" % (
                r["found"], r["ref"], r["err"], "used in the fit" if r["use"] else "removed from the fit"))

    def point_menu(self, x, y, card):
        """Right click: on a point, remove that calibrant from the fit (or use
        it again); anywhere, the calibrants near the click and the bulk
        actions."""
        items = []
        i = self._hit(card, x, y)
        if i is not None:
            r = self.rows[i]
            self.table.select(i)
            self._draw_highlight()  # ringed at once, while the menu is open
            if r["use"]:
                items.append(("Remove this point from the fit: m/z %.4f (%+.2f ppm)" % (r["ref"], r["err"]),
                              lambda: self.set_use([i], False)))
            else:
                items.append(("Use this point in the fit again: m/z %.4f" % r["ref"], lambda: self.set_use([i], True)))
        else:
            j = self._nearest(x, y)
            if j is not None:
                r = self.rows[j]
                items.append((("Remove the nearest calibrant, m/z %.4f" if r["use"] else
                               "Use the nearest calibrant again, m/z %.4f") % r["ref"],
                              lambda: self.set_use([j], not r["use"])))
        items += self.bulk_menu() + [(None, None), ("Full view", card.reset_view)]
        return items

    def bulk_menu(self):
        items = [(None, None), ("Remove calibrants with an error above\u2026", self.remove_above)]
        if any(r["found"] is not None and not r["use"] for r in self.rows):
            items.append(("Use every calibrant again", lambda: self.set_use(range(len(self.rows)), True)))
        return items

    def table_menu(self):
        sel = self.table.selected()
        items = []
        if sel:
            items += [("Remove the selected calibrants from the fit (Delete)", lambda: self.set_use(sel, False)),
                      ("Use the selected calibrants again", lambda: self.set_use(sel, True))]
        return items + self.bulk_menu() + [(None, None)]

    def remove_above(self):
        v = U.ask_form(self, "Remove calibrants", [
            dict(key="ppm", label="Error above", value="3", unit="ppm"),
            dict(key="which", label="Error", kind="choice", choices=["after calibration", "as recorded",
                                                                     "left out (cross-validated)"], value=0,
                 width=190)], note=None, ok="Remove")
        if v is None:
            return
        try:
            lim = abs(float(v["ppm"].replace(",", ".")))
        except ValueError:
            return
        used = [i for i, r in enumerate(self.rows) if r["found"] is not None and r["use"]]
        idx = []
        for n, i in enumerate(used):
            r = self.rows[i]
            if v["which"] == 1:
                e = r["err"]
            elif v["which"] == 2 and self.cal is not None and self.cal.cv_err is not None:
                e = self.cal.cv_err[n]
            elif self.cal is not None:
                e = (float(self.cal(np.array([r["found"]]))[0]) - r["ref"]) / r["ref"] * 1e6
            else:
                e = r["err"]
            if abs(e) > lim:
                idx.append(i)
        if not idx:
            self.SetStatusText("No calibrant has an error above %g ppm" % lim)
            return
        self.set_use(idx, False)

    def pick_model(self, i):
        keys = getattr(self, "_model_keys", [])
        if 0 <= i < len(keys):
            k, order = keys[i]
            # the row was fitted with the order chosen here (or automatic): keep it, so the
            # calibration gets the numbers the row showed
            self.model.SetSelection(k)
            self.refit()

    def close(self):
        if self.cal is not None and self.applied_sig != self._sig() and self.cal.rms <= self.cal.rms_before:
            msg = "Apply this calibration to the file (%s, %.2f ppm)?" % (self.cal.name(), self.cal.rms)
            if self.applied_sig is not None:
                msg = ("The settings changed after Calibrate. Apply the new calibration (%s, %.2f ppm)?"
                       % (self.cal.name(), self.cal.rms))
            a = wx.MessageBox(msg, "Calibration", wx.YES_NO | wx.CANCEL | wx.ICON_QUESTION, self)
            if a == wx.CANCEL:
                return
            if a == wx.YES:
                self.calibrate()
        self.EndModal(wx.ID_OK)

    # -------------------------------------------------------------- plots
    def plot_spec(self, keep=False):
        cs = self.card_spec
        view = cs.ax.get_xlim() if (keep and cs.full) else None
        cs.reset()
        U.style_axes(cs.ax, "m/z", "Intensity")
        if self.spec is None or not len(self.spec):
            cs.full = None
            cs.set_title("Calibrant spectrum", "")
            cs.draw()
            return
        s = self.spec
        U.plot_trace(cs.ax, s[:, 0], s[:, 1], color=U.TRACE, lw=U.LW_SPEC)
        cs.full = (float(s[0, 0]), float(s[-1, 0]))
        cs.ax.set_xlim(*(view or cs.full))
        cs.autoscale_y(pad=0.16)
        x0, x1 = cs.ax.get_xlim()
        sep = (x1 - x0) / 16.0
        last = None
        found = [r for r in self.rows if r["found"] is not None]
        for r in self.rows:
            if r["found"] is None:
                ln = cs.ax.axvline(r["ref"], color="#C4CAD3", lw=0.7, ls=(0, (2, 2)), zorder=1)
                ln._no_scale = True
                continue
            col = CAL_COL if r["use"] else "#9AA3AF"
            pt, = cs.ax.plot([r["found"]], [r["intensity"]], ls="none", marker="v", ms=6, mfc=col, mec="white",
                             mew=0.6, zorder=5)
            pt._no_scale = True
            if x0 <= r["found"] <= x1 and (last is None or r["found"] - last >= sep):
                cs.ax.annotate("%.4f" % r["ref"], (r["found"], r["intensity"]), xytext=(0, 7),
                               textcoords="offset points", ha="center", va="bottom", fontsize=7, color=col, zorder=6)
                last = r["found"]
        cs.set_title("Calibrant spectrum", "%s: %d of %d reference ions found" % (
            self.calib.GetStringSelection(), len(found), len(self.rows)))
        cs._hl_arts = []
        self._hl_card(cs)
        cs.draw()

    def plot_errors(self):
        ce = self.card_err
        ce.reset()
        U.style_axes(ce.ax, "m/z", "Error (ppm)")
        found = [r for r in self.rows if r["found"] is not None]
        if not found:
            ce.full = None
            ce.set_title("Mass error", "no calibrant peaks found")
            ce.set_legend([])
            ce.draw()
            return
        x = np.array([r["found"] for r in found])
        ref = np.array([r["ref"] for r in found])
        before = np.array([r["err"] for r in found])
        use = np.array([r["use"] for r in found])
        ce.ax.axhline(0, color=U.INK, lw=0.7, zorder=1)._no_scale = True
        ce.ax.plot(x[use], before[use], ls="none", marker="o", ms=5, mfc="white", mec="#8A94A3", mew=1.1, zorder=3)
        ce.ax.plot(x[~use], before[~use], ls="none", marker="x", ms=5, color="#B8BEC7", zorder=3)
        legend = [("#8A94A3", "as recorded")]
        vals = [before]
        if self.cal is not None:
            after = (self.cal(x) - ref) / ref * 1e6
            vals.append(after)
            ce.ax.plot(x[use], after[use], ls="none", marker="o", ms=5, mfc="#2A62C4", mec="white", mew=0.6,
                       zorder=4)
            ce.ax.plot(x[~use], after[~use], ls="none", marker="x", ms=5, color="#9DB7E6", zorder=4)
            lo, hi = float(self.spec[0, 0]), float(self.spec[-1, 0])
            xx = np.linspace(lo, hi, 400)
            yy = self.cal(xx)
            ln, = ce.ax.plot(xx, (xx - yy) / yy * 1e6, color="#E0602F", lw=0.8, alpha=0.85, zorder=2)
            ln._no_scale = True
            ce.ax.axvspan(self.cal.lo, self.cal.hi, color=CAL_COL, alpha=0.06, lw=0, zorder=0)
            legend += [("#2A62C4", "calibrated"), ("#E0602F", "model")]
            ce.set_title("Mass error", "%.2f ppm RMS (as recorded %.2f ppm)" % (self.cal.rms, self.cal.rms_before))
        else:
            ce.set_title("Mass error", "no fit")
        ce.set_legend(legend)
        ce.full = (float(self.spec[0, 0]), float(self.spec[-1, 0]))
        ce.ax.set_xlim(*ce.full)
        v = np.concatenate([np.asarray(a, float) for a in vals])
        lim = max(2.0, float(np.nanmax(np.abs(v))) * 1.25)
        ce.ax.set_ylim(-lim, lim)
        ce.ylock = (-lim, lim)
        ce._hl_arts = []
        self._hl_card(ce)
        ce.draw()

    def fill_tables(self):
        rows = []
        k = 0
        cv = {}
        if self.cal is not None and self.cal.cv_err is not None:
            used = [r for r in self.rows if r["use"] and r["found"] is not None]
            cv = {id(r): e for r, e in zip(used, self.cal.cv_err)}
        for r in self.rows:
            if r["found"] is None:
                rows.append(("", "%.5f" % r["ref"], "not found", "", "", "", ""))
                continue
            a = (float(self.cal(np.array([r["found"]]))[0]) - r["ref"]) / r["ref"] * 1e6 if self.cal is not None \
                else None
            rows.append(("yes" if r["use"] else "no", "%.5f" % r["ref"], "%.5f" % r["found"], U._g(r["intensity"]),
                         "%+.2f" % r["err"], "%+.2f" % a if a is not None else "",
                         "%+.2f" % cv[id(r)] if id(r) in cv else ""))
            k += 1
        n_used = sum(1 for r in self.rows if r["use"] and r["found"] is not None)
        self.table.set_rows(rows, "%d of %d used in the fit" % (n_used, len(self.rows)))
        # models compared on the same calibrants
        used = [r for r in self.rows if r["use"] and r["found"] is not None]
        mrows, keys = [], []
        best = None
        cfg = self.cfg()
        for m in range(len(K.MODELS)):
            try:
                c = self.cp.fit(self.rows, dict(cfg, model=m, order=cfg["order"] if m == 5 else 0))
            except (ValueError, np.linalg.LinAlgError):
                mrows.append((K.MODELS[m], "", "too few peaks"))
                keys.append((m, None))
                continue
            short = ["Single point", "Linear", "Quadratic", "Cubic", "TOF"][m] if m < 5 else \
                "HPC (order %d)" % c.hpc_order
            mrows.append((short, "%.2f" % c.rms, "%.2f" % c.cv_rms if c.cv_rms is not None else ""))
            keys.append((m, getattr(c, "hpc_order", None)))
            if c.cv_rms is not None and (best is None or c.cv_rms < best[1]):
                best = (m, c.cv_rms)
        self._model_keys = keys
        info = ("lowest CV: %s" % K.MODELS[best[0]].split(" (")[0]) if best else ""
        if used:
            info += "; double click to use"
        self.models.set_rows(mrows, info)
        self.models.select(self.model.GetSelection())


# ==========================================================================
# window
# ==========================================================================
HELP_TEXT = """HRMS Analysis reads Bruker QTOF data (.d folders with analysis.baf or
analysis.tsf), timsTOF data (analysis.tdf, ion mobility summed, MS/MS left
out), Agilent MassHunter .D folders, Waters .raw folders, Thermo .raw files
and mzML files. Sciex .wiff files: convert them to mzML with ProteoWizard
msconvert.

Chromatograms and spectra
  * Drag on a chromatogram to average a range; Shift + drag for the
    background. Right click > Add for a mass chromatogram.

Calibration
  1. Drag across the calibrant peak (sodium formate; found and marked
     green when the file opens).
  2. Right click the spectrum > Calibrate on this spectrum.
  3. Choose calibrant and model; Calibrate applies it to the whole file.
     Delete or right click removes calibrants from the fit.
  4. The top bar shows the error; click it to open the window again.

Deconvolution
  * Right click the spectrum > Deconvolute this spectrum.
  * Isotope patterns need resolved isotopes, a mass step of 0.01 to
    0.02 Da and a narrow mass range.
  * Mass shifts: right click the mass spectrum > Mass shifts.

Exact mass
  * Right click the spectrum > Check a formula.
  * Right click a peak > Find formulas.

Methods, undo, files and reports
  * Method (top bar): save, apply, export or import every setting.
  * Undo and Redo: Ctrl+Z and Ctrl+Y. F1 lists every shortcut.
  * Report... (Ctrl+R): compound report, deconvolution report or
    Supporting Information page, as PDF or Word.
  * Everything saved for a data file goes into <name>_analysis next to it.

Mouse and images
  * Ctrl + drag zooms a box, double click gives the full view.
  * Copy image (Ctrl+C) and Save image as (Ctrl+S): 8.5 x 6 cm, 300 dpi.

Bruker data are read with Bruker's Baf2Sql library (analysis.baf; freely
distributed by Bruker, see _portable\\baf2sql; it writes analysis.sqlite
into the .d folder the first time) and Bruker's TDF SDK (analysis.tsf,
analysis.tdf; _portable\\timsdata). The m/z values as recorded are those
of the last calibration saved with the file.
This software uses TDF Software Development Kit software. Copyright (c)
2019 by Bruker Daltonik GmbH. All rights reserved.
Agilent data are read with the rainbow library (LGPL), Waters data with
Waters' MassLynx library and Thermo data with Thermo's RawFileReader.

HPC follows the published High Precision Calibration (Gobom et al., Anal.
Chem. 2002, 74, 3915); results can differ slightly from DataAnalysis.

Bayesian deconvolution: M. T. Marty et al., Anal. Chem. 2015, 87, 4370."""


class HRMSFrame(U.PostrunFrame):
    TITLE = TITLE
    HELP = HELP_TEXT
    START_HINT = "Open a .d or .raw folder, a .raw file or an mzML file (Ctrl+O)"
    FOLDER_KEY = "hrms_folder"

    def __init__(self, path=None):
        U.PostrunFrame.__init__(self, path)
        chip = self.chips["cal"]
        chip.SetCursor(wx.Cursor(wx.CURSOR_HAND))
        chip.SetToolTip("Click for the calibration window")
        chip.Bind(wx.EVT_LEFT_UP, lambda e: self.ms.cal.open_details())

    def build_tabs(self, book):
        self.ms = HRMSTab(book, self)
        return [("Mass spectrometry", self.ms)]

    def chip_specs(self):
        return [("file", "File"), ("ms", "MS"), ("cal", "Calibration")]

    def accepts(self, path):
        return hrms_data.find_d_folder(path) is not None or path.lower().endswith(".mzml") or \
            hrms_vendor.detect_kind(path) in hrms_vendor.VENDOR_KINDS

    def browse_accepts(self, path):
        low = path.lower().rstrip("\\/")
        if low.endswith(".mzml"):
            return os.path.isfile(path)
        if low.endswith(".raw"):
            return hrms_vendor.detect_kind(path) in ("thermo_raw", "waters_raw")
        return low.endswith(".d") and os.path.isdir(path) and (
            hrms_data.is_bruker_d(path) or os.path.isfile(os.path.join(path, "AcqData", "MSScan.bin")))

    def write_spectrum(self, path, spec):
        hrms_data.write_spectrum_txt(path, spec)

    def on_open(self, e=None):
        menu = wx.Menu()
        a = menu.Append(wx.ID_ANY, "Bruker or Agilent .d folder…")
        c = menu.Append(wx.ID_ANY, "Waters .raw folder…")
        b = menu.Append(wx.ID_ANY, "mzML or Thermo .raw file…")
        self.Bind(wx.EVT_MENU, lambda ev: self.open_d(), a)
        self.Bind(wx.EVT_MENU, lambda ev: self.open_waters(), c)
        self.Bind(wx.EVT_MENU, lambda ev: self.open_mzml(), b)
        btn = e.GetEventObject() if e is not None else None
        if isinstance(btn, wx.Window) and btn is not self and not isinstance(btn, wx.Frame):
            U.popup_menu(btn, menu, wx.Point(0, btn.GetSize()[1]))
        else:
            self.open_d()
        menu.Destroy()

    def open_d(self):
        dlg = wx.DirDialog(self.window(), "Open a Bruker or Agilent .d folder", defaultPath=self.folder(),
                           style=wx.DD_DEFAULT_STYLE | wx.DD_DIR_MUST_EXIST)
        try:
            if U.show_open_dialog(dlg) != wx.ID_OK:
                return
            path = dlg.GetPath()
        finally:
            dlg.Destroy()
        d = hrms_data.find_d_folder(path) or hrms_vendor.agilent_folder(path)
        if d is None:
            wx.MessageBox("%s is not a Bruker or Agilent .d folder." % path, TITLE, wx.ICON_WARNING)
            return
        self.load(d)

    def open_waters(self):
        dlg = wx.DirDialog(self.window(), "Open a Waters .raw folder", defaultPath=self.folder(),
                           style=wx.DD_DEFAULT_STYLE | wx.DD_DIR_MUST_EXIST)
        try:
            if U.show_open_dialog(dlg) != wx.ID_OK:
                return
            path = dlg.GetPath()
        finally:
            dlg.Destroy()
        d = hrms_vendor.waters_folder(path)
        if d is None:
            wx.MessageBox("%s is not a Waters .raw folder." % path, TITLE, wx.ICON_WARNING)
            return
        self.load(d)

    def open_mzml(self):
        dlg = wx.FileDialog(self.window(), "Open data files (several: Ctrl or Shift + click)", defaultDir=self.folder(),
                            wildcard="mzML or Thermo .raw (*.mzML;*.raw)|*.mzML;*.mzml;*.raw;*.RAW|"
                                     "Sciex .wiff (*.wiff)|*.wiff;*.wiff2|All files (*.*)|*.*",
                            style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST | wx.FD_MULTIPLE)
        try:
            if U.show_open_dialog(dlg) == wx.ID_OK:
                for p in dlg.GetPaths():
                    self.load(p)
        finally:
            dlg.Destroy()

    def read_data(self, path, progress):
        try:
            data = hrms_data.open_hrms(path, progress)
        except Exception as ex:
            import traceback
            traceback.print_exc()
            return {"ms": None}, [str(ex)]
        try:  # still in the reading thread: the calibrant segment (mass chromatograms)
            self.ms.cal.precompute(data)
        except Exception as ex:
            print("[%s] calibrant segment not looked for: %s" % (TITLE, ex))
        return {"ms": data}, []

    def data_loaded(self, path, data, errors):
        ms = data.get("ms")
        old = self.ms_data
        self.ms_data = ms
        if old is not None and old is not ms:
            try:
                old.close()
            except Exception:
                pass
        pols = [e["polarity"] for e in ms.events]
        self.chips["ms"].SetValue("ESI " + " / ".join({"+": "+", "-": "−"}.get(p, "?") for p in pols))
        self.chips["cal"].SetValue(self._file_cal(ms))
        self.ms.set_data(ms)
        info = [ms.summary()] + (["; ".join(errors)] if errors else [])
        self.SetStatusText(" | ".join(info), 0)
        self.set_file_status(os.path.basename(path.rstrip("\\/")))

    @staticmethod
    def _file_cal(d):
        note = getattr(d, "recal_note", "")
        return ("file recal. (%s)" % note) if note else "as recorded"

    def calibration_changed(self):
        d = self.ms_data
        if d is None:
            return
        if d.calib is not None:
            # with both polarities in the file: which ions each calibration is for (each
            # polarity can have its own)
            both = len(set(e["polarity"] for e in d.events)) > 1
            parts = []
            for cal, mine in _cals_of(d):
                which = (" (%s)" % {"+": "+", "-": "\u2212"}.get(mine[0], mine[0])) if both and len(mine) == 1 else ""
                parts.append("%.2f ppm, %s%s" % (cal.rms, cal.short(), which))
            self.chips["cal"].SetValue("; ".join(parts) + (" + file recal." if getattr(d, "recal_note", "") else ""))
        else:
            self.chips["cal"].SetValue(self._file_cal(d))
        try:  # the chips left of it move when its text changes: draw the whole bar again
            self.window().chips["cal"].GetParent().Refresh()
        except Exception:
            pass
        ms = self.ms
        ms.on_update()
        if ms.avg and ms.avg[1] != ms.avg[0]:
            ms.compute_average()
        elif ms.pick_t is not None:
            ms.chrom_pick(ms.pick_t, None, None)
        for col in ms.cols:
            if col.get("overlay"):
                ms.check_formula()
                break


def open_window(path=None):
    """The HRMS Analysis window; path: a data file (.d or .raw folder, .raw or mzML file) or a list
    of them (each opens as a file of the window)."""
    import time
    t0 = time.perf_counter()
    paths = [p for p in (path if isinstance(path, (list, tuple)) else [path]) if p]
    f = HRMSFrame(paths[0] if paths else None)
    for p in paths[1:]:
        wx.CallAfter(f.load, p)
    f.Show()
    f.Raise()
    print("[HRMS Analysis] window built in %.1f s" % (time.perf_counter() - t0))
    return f


def run(path=None):
    app = wx.App(False)
    app.SetAppName(U.APP)
    open_window(path)
    app.MainLoop()


if __name__ == "__main__":
    import sys
    run(sys.argv[1] if len(sys.argv) > 1 else None)

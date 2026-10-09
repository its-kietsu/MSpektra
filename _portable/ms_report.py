"""
Reports of MSpektra (LCMS Analysis and HRMS Analysis).

Five kinds, chosen when the report is made:
  hrms    HRMS compound report     chromatograms, spectrum, isotope pattern,
                                   exact mass, deconvolution, calibration
  lcms    LC-MS purity report      PDA chromatogram with peak table, UV and
                                   MS spectra of the main peak, TICs
  deconv  Deconvolution report     spectrum with fit, zero charge spectrum,
                                   masses (with expected masses), settings
  si      Supporting Information   one figure (a to d), caption and a
                                   characterisation text to paste
  compare Comparison report        the runs of the Compare view of LCMS
                                   Analysis as shown, the runs, their peaks
Each as PDF (reportlab) or Word (python-docx). The content is taken from the
file shown in the window, as it is on screen (ranges, peaks, results); the
comparison report from the Compare view of the window (all its runs).

The report libraries are shipped as wheels in _portable\\wheels and unpacked
the first time a report is made (no installation, no admin rights).
"""
import os
import re
import sys
import html
import tempfile
import datetime

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
KINDS = [
    ("hrms", "HRMS compound report", "Chromatograms, mass spectrum, masses and calibration."),
    ("lcms", "LC-MS purity report", "Chromatogram with peak table, UV and MS spectra of the main peak, TICs."),
    ("deconv", "Deconvolution report", "Spectrum with fit, deconvoluted spectrum, masses and settings."),
    ("si", "Supporting Information page", "One figure with caption and a characterisation text."),
    ("compare", "Comparison report", "Comparison plot, runs, peaks and a methods text."),
]
KIND_NAMES = {k: n for k, n, _ in KINDS}


# ==========================================================================
# libraries and fonts
# ==========================================================================
def ensure_libs():
    """reportlab and python-docx importable (unpacked from the wheels in
    _portable\\wheels the first time, see pylibs_loader). Returns the folder used."""
    import pylibs_loader
    try:
        return pylibs_loader.ensure_pylibs("reportlab", "docx")
    except RuntimeError as ex:
        raise RuntimeError("the report libraries could not be loaded (%s)" % ex)


def _font_files():
    """(regular, bold, italic, bold italic) for a sans and a serif family:
    Arial and Times New Roman on Windows, else Liberation, else DejaVu."""
    win = os.path.join(os.environ.get("WINDIR", r"C:\Windows"), "Fonts")
    lib = "/usr/share/fonts/truetype/liberation"
    try:
        import matplotlib
        dv = os.path.join(os.path.dirname(matplotlib.__file__), "mpl-data", "fonts", "ttf")
    except Exception:
        dv = ""
    cands = {
        "sans": [(win, ("arial.ttf", "arialbd.ttf", "ariali.ttf", "arialbi.ttf")),
                 (lib, ("LiberationSans-Regular.ttf", "LiberationSans-Bold.ttf", "LiberationSans-Italic.ttf",
                        "LiberationSans-BoldItalic.ttf")),
                 (dv, ("DejaVuSans.ttf", "DejaVuSans-Bold.ttf", "DejaVuSans-Oblique.ttf",
                       "DejaVuSans-BoldOblique.ttf"))],
        "serif": [(win, ("times.ttf", "timesbd.ttf", "timesi.ttf", "timesbi.ttf")),
                  (lib, ("LiberationSerif-Regular.ttf", "LiberationSerif-Bold.ttf", "LiberationSerif-Italic.ttf",
                         "LiberationSerif-BoldItalic.ttf")),
                  (dv, ("DejaVuSerif.ttf", "DejaVuSerif-Bold.ttf", "DejaVuSerif-Italic.ttf",
                        "DejaVuSerif-BoldItalic.ttf"))],
    }
    out = {}
    for fam, lst in cands.items():
        for d, names in lst:
            paths = [os.path.join(d, n) for n in names]
            if d and all(os.path.isfile(p) for p in paths):
                out[fam] = paths
                break
    return out


_FONTS_DONE = {}


def _register_fonts():
    if _FONTS_DONE:
        return _FONTS_DONE
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
    from reportlab.lib.fonts import addMapping
    files = _font_files()
    names = {}
    for fam, tag in (("sans", "RSans"), ("serif", "RSerif")):
        if fam in files:
            for suffix, path in zip(("", "-Bold", "-Italic", "-BoldItalic"), files[fam]):
                pdfmetrics.registerFont(TTFont(tag + suffix, path))
            addMapping(tag, 0, 0, tag)
            addMapping(tag, 1, 0, tag + "-Bold")
            addMapping(tag, 0, 1, tag + "-Italic")
            addMapping(tag, 1, 1, tag + "-BoldItalic")
            names[fam] = tag
        else:  # built in (no Greek letters)
            names[fam] = "Helvetica" if fam == "sans" else "Times-Roman"
    _FONTS_DONE.update(names)
    return names


# ==========================================================================
# figures (matplotlib, without any window)
# ==========================================================================
INK, NAVY, ACCENT, ORANGE, GREEN, GREY = "#1B2330", "#1F3F7A", "#2A62C4", "#E0602F", "#1E9E6A", "#8A94A6"
AX, BLUE = "#000000", "#1D4ED8"  # figures: black axes, blue traces
CM = 1 / 2.54
W_FULL = 17.8


def _mpl():
    import matplotlib
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    return matplotlib, Figure, FigureCanvasAgg


def new_fig(w_cm, h_cm, ncols=1, nrows=1, **kw):
    matplotlib, Figure, FigureCanvasAgg = _mpl()
    with matplotlib.rc_context(_RC):
        f = Figure(figsize=(w_cm * CM, h_cm * CM), facecolor="white")
        FigureCanvasAgg(f)
        axs = f.subplots(nrows, ncols, **kw)
    return f, axs


_RC = {"font.family": ["Arial", "Liberation Sans", "DejaVu Sans"], "font.size": 7.5, "axes.linewidth": 0.7,
       "axes.edgecolor": AX, "axes.labelcolor": AX, "xtick.color": AX, "ytick.color": AX,
       "mathtext.default": "regular", "pdf.fonttype": 42}


def style(ax, xlabel="", ylabel=""):
    from matplotlib.ticker import AutoMinorLocator, MaxNLocator, ScalarFormatter
    for s in ax.spines.values():
        s.set_linewidth(0.7)
    ax.tick_params(which="both", direction="out", top=False, right=False, width=0.6)
    ax.tick_params(which="major", length=3.0, labelsize=7, pad=1.5)
    ax.tick_params(which="minor", length=1.6)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.xaxis.set_major_locator(MaxNLocator(nbins=8, steps=[1, 2, 2.5, 5, 10]))
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10], min_n_ticks=3))
    f = ScalarFormatter(useMathText=True)
    f.set_powerlimits((-3, 4))
    ax.yaxis.set_major_formatter(f)
    ax.yaxis.get_offset_text().set_fontsize(6.5)
    ax.set_xlabel(xlabel, fontsize=7.5, labelpad=1.5)
    ax.set_ylabel(ylabel, fontsize=7.5, labelpad=2)


def xbins(ax, n):
    from matplotlib.ticker import MaxNLocator
    ax.xaxis.set_major_locator(MaxNLocator(nbins=n, steps=[1, 2, 2.5, 5, 10]))


def headroom(ax, y, pad=0.2):
    y = np.asarray(y, float)
    top = float(np.nanmax(y)) if len(y) else 1.0
    ax.set_ylim(0, (top if top > 0 else 1.0) * (1 + pad))


_APEX = {"on": False}  # HRMS profile spectra: label positions at the apex of the peak


def _apex_x(x, y, i):
    """Label position: the apex of the profile peak (HRMS), found as in the
    window (formula check, clicked peaks, calibrant search), else the point."""
    if not _APEX["on"]:
        return float(x[i])
    try:
        import hrms_calib
        lo, hi = max(0, i - 25), min(len(x), i + 26)
        return hrms_calib.profile_apex(x[lo:hi], y[lo:hi], i - lo)
    except Exception:
        return float(x[i])


def local_max(x, y, n=6, sep=None, xlim=None, min_rel=0.0):
    x, y = np.asarray(x, float), np.asarray(y, float)
    if xlim:
        m = (x >= xlim[0]) & (x <= xlim[1])
        x, y = x[m], y[m]
    if len(x) < 3:
        return []
    sep = sep or (x[-1] - x[0]) / 18.0
    k = np.where((y[1:-1] >= y[:-2]) & (y[1:-1] > y[2:]))[0] + 1
    top = float(y.max()) if len(y) else 0
    out = []
    for i in k[np.argsort(y[k])[::-1]]:
        if len(out) >= n or y[i] < min_rel * top or y[i] <= 0:
            break
        if all(abs(x[i] - a) > sep for a, _ in out):
            out.append((_apex_x(x, y, i), float(y[i])))
    return out


def label_peaks(ax, x, y, n=6, fmt="%.4f", sep=None, xlim=None, min_rel=0.0, fs=6.3):
    x0, x1 = xlim or ax.get_xlim()
    for px, py in local_max(x, y, n, sep, (x0, x1), min_rel):
        f = (px - x0) / ((x1 - x0) or 1.0)
        ha = "left" if f < 0.06 else "right" if f > 0.94 else "center"
        ax.annotate(fmt % px, (px, py), xytext=(1 if ha == "left" else -1 if ha == "right" else 0, 2),
                    textcoords="offset points", ha=ha, va="bottom", fontsize=fs, color="#222222")


def save_fig(f, path, **adj):
    if adj:
        f.subplots_adjust(**adj)
    f.savefig(path, dpi=300, facecolor="white")
    return path


TEXT_SCALE = 0.75  # text in the report / text on screen (axis titles: 10 pt on screen, 7.5 pt here)


def _tile_png(card, path, w_cm, h_cm, **kw):
    """Image of a tile of the window exactly as it shows it (same view, traces, labels, markers,
    overlays, shaded ranges), for a report figure of w_cm x h_cm. None when there is no such tile
    on screen (or its image fails): the report then draws the figure from the data itself."""
    if card is None or not hasattr(card, "snapshot"):
        return None
    try:
        if hasattr(card, "shown") and not card.shown():
            return None
        return card.snapshot(path, w_cm, h_cm, TEXT_SCALE, **kw)
    except Exception as ex:
        print("report: no image of the tile %r (%s); the report draws it" % (getattr(card, "title", ""), ex))
        return None


def _is_mark(a):
    return bool(getattr(a, "_range_mark", False))


def _mark_opts(mark):
    """Tile image options: the shaded time ranges (spectrum, background, calibrant scans) kept
    (as on screen) or left out (report window: the shading option unticked)."""
    return {} if mark else {"hide": _is_mark, "legend_keep": lambda lab: "calibrant" not in lab}


def _sticks_or_line(ax, s, lw=0.6):
    """Profile spectra as a line; centroid (stick) spectra as sticks."""
    s = np.asarray(s, float)
    if len(s) > 2 and np.median(np.diff(s[:, 0])) > 0.2 and np.mean(s[:, 1] == 0) < 0.3:
        ax.vlines(s[:, 0], 0, s[:, 1], color=BLUE, lw=lw)
    else:
        ax.plot(s[:, 0], s[:, 1], color=BLUE, lw=lw)


# ==========================================================================
# content blocks (rendered to PDF or Word)
# ==========================================================================
def esc(t):
    return html.escape(str(t), quote=False)


class Report(object):
    def __init__(self, kind, meta):
        self.kind, self.meta, self.blocks = kind, meta, []

    def add(self, *b):
        self.blocks.append(b)


_TITLES = ("<b>Masses</b>", "<b>Isotope peaks</b>", "<b>Mass shifts</b>", "<b>Degree of conjugation</b>")


def _plain_report(rep, fields):
    """Keep measurements, tables and short headings (HRMS and deconvolution
    reports): no captions or text blocks, only the table titles; a heading
    with nothing under it is left out."""
    def clean(blocks):
        out = []
        for b in blocks:
            kind = b[0]
            if kind in ("caption", "text") or (kind == "small" and b[1] not in _TITLES):
                continue
            if kind == "pair":
                out.append(("pair", clean(b[1]), clean(b[2]), b[3]))
                continue
            out.append(b)
        return [b for i, b in enumerate(out) if not (b[0] == "section" and
                (i == len(out) - 1 or out[i + 1][0] == "section"))]
    rep.blocks = clean(rep.blocks)
    return rep


def g(v):
    try:
        return ("%.2e" % float(v)).replace("e+0", "e").replace("e+", "e").replace("e-0", "e-")
    except (TypeError, ValueError):
        return str(v)


def fmt_date(d):
    if isinstance(d, datetime.datetime):
        return d.strftime("%d.%m.%Y %H:%M")
    return str(d or "")


# ==========================================================================
# what is on screen
# ==========================================================================
def sample_info(frame):
    """Sample information of the file shown (LabSolutions or Bruker)."""
    doc = frame.active
    cached = doc.attrs.get("sample_info")
    if cached is not None:
        return cached
    info = {}
    path = frame.path or ""
    try:
        if path.lower().endswith(".lcd"):
            import lcms_data
            info = lcms_data.read_sample_info(path)
        else:
            props = getattr(frame.ms_data, "props", {}) or {}
            for key, names in (("sample_name", ("SampleName", "Sample Name", "SampleID")),
                               ("operator", ("OperatorName", "Operator")),
                               ("instrument", ("InstrumentName", "Instrument")),
                               ("method_file", ("MethodName", "Method", "AcquisitionMethod")),
                               ("acquired", ("AcquisitionDateTime", "AcquisitionDate", "StartTime"))):
                for n in names:
                    if props.get(n):
                        info[key] = str(props[n])
                        break
            if info.get("acquired"):  # ISO 8601 as Bruker writes it, e.g. 2026-09-15T12:14:42.300+02:00
                info["acquired"] = _iso_date(info["acquired"])
    except Exception:
        info = {}
    if "acquired" not in info and path:
        try:
            info["acquired_file"] = datetime.datetime.fromtimestamp(os.path.getmtime(path))
        except OSError:
            pass
    doc.sample_info = info
    return info


def _iso_date(text):
    """datetime of an ISO 8601 date (local time as written), or the text."""
    m = re.match(r"\s*(\d{4})-(\d{2})-(\d{2})[T ](\d{2}):(\d{2})(?::(\d{2}))?", str(text))
    if not m:
        return text
    try:
        return datetime.datetime(*[int(v) for v in m.groups(default="0")])
    except ValueError:
        return text


def _instrument(frame):
    path = frame.path or ""
    if path.lower().endswith(".lcd"):
        try:
            import olefile
            f = olefile.OleFileIO(path)
            names = ["/".join(e) for e in f.listdir(streams=False, storages=True)]
            f.close()
            for n in names:
                m = re.search(r"ShimadzuLCMS(\d+)", n)
                if m:
                    return "Shimadzu LCMS-%s" % m.group(1)
        except Exception:
            pass
        return "Shimadzu LC-MS (LabSolutions)"
    try:
        si = sample_info(frame)
    except Exception:  # (the Compare report: a path only)
        import lcms_sources
        si = lcms_sources.read_sample_info(path)
    props = getattr(getattr(frame, "ms_data", None), "props", None) or {}
    if si.get("instrument") or props.get("InstrumentName"):
        return si.get("instrument") or str(props["InstrumentName"])
    import data_formats
    return "Bruker QTOF" if data_formats.detect(path)[1] == "bruker_d" else ""


def _spec_col(tab, e=None):
    cols = [c for c in tab.cols if c.get("spec") is not None and len(c["spec"])]
    if e is not None:
        cols = [c for c in cols if c["e"] == e] or cols
    over = [c for c in cols if c.get("overlay")]
    return (over or cols or [None])[0]


def _traces(tab, e):
    out = []
    for v in tab.views:
        if v.key == e:
            for tr in v.traces:
                if tr.get("t") is not None and len(tr["t"]):
                    out.append(tr)
    return out


def _formula_info(col):
    ov = col.get("overlay") if col else None
    if not ov:
        return None
    info = dict(ov.get("info") or {})
    info["pattern"] = ov.get("pattern")
    info["scale"] = ov.get("scale")
    if "name" not in info:  # older results: read the text
        lines = (ov.get("text") or "").split("\n")
        info["name"] = lines[0] if lines else ""
        m = re.search(r"m/z ([\d.]+) \(calc", ov.get("text") or "")
        if m:
            info["mz"] = float(m.group(1))
        m = re.search(r"found ([\d.]+), ([+\-−][\d.]+) ppm", ov.get("text") or "")
        if m:
            info["found"] = float(m.group(1))
            info["err"] = float(m.group(2).replace("−", "-"))
        m = re.search(r"isotope match (\d+)", ov.get("text") or "")
        if m:
            info["match"] = float(m.group(1))
    return info


def availability(frame):
    """kind -> (available, reason if not)."""
    out = {}
    lc = hasattr(frame, "pda") and getattr(frame, "pda", None) is not None
    tab = getattr(frame, "ms", None)
    data = getattr(frame, "ms_data", None)
    col = _spec_col(tab) if (tab is not None and data is not None) else None
    pda_ok = lc and getattr(frame, "pda_data", None) is not None
    if lc:
        out["hrms"] = (False, "only in HRMS Analysis")
    else:
        out["hrms"] = (col is not None, "show a spectrum first")
    if lc:
        out["lcms"] = (pda_ok or data is not None, "open a data file first")
    else:
        out["lcms"] = (False, "only in LCMS Analysis")
    dec = getattr(tab, "dec", None) if tab is not None else None
    res = getattr(dec, "result", None)
    hidden = res is None and bool(getattr(dec, "results", None))  # every result hidden in the list of open files
    out["deconv"] = (res is not None, "show a hidden result first (eye in the file list)" if hidden else
                     "deconvolute a spectrum first")
    out["si"] = (pda_ok or col is not None, "open a data file and show a spectrum first")
    # the Compare view belongs to the window: available whichever view is shown
    drawn = _compare_drawn(frame)[1]
    if not lc:
        out["compare"] = (False, "only in LCMS Analysis")
    else:
        out["compare"] = (bool(drawn), "show two or more files in the Compare view first")
    return out


# ==========================================================================
# building the four reports
# ==========================================================================
PAL = [BLUE, "#E0602F", "#1E9E6A", "#8E44AD", "#C0392B", "#2A62C4"]


def _event_name(msd, e):
    pol = None
    try:
        pol = msd.events[e].get("polarity") if hasattr(msd, "events") else None
    except Exception:
        pol = None
    if pol == "+":
        return "ESI(+)"
    if pol == "-":
        return "ESI(−)"
    try:
        return msd.event_label(e)
    except Exception:
        return "event %d" % (e + 1)


def figure_options(frame, kind):
    """[(key, label, default)]: the figures a report of this kind can hold
    for the file shown (the report window lists them as check boxes)."""
    lc = getattr(frame, "pda", None) is not None
    msd = getattr(frame, "ms_data", None)
    pda_ok = lc and getattr(frame, "pda_data", None) is not None
    events = []
    if msd is not None:
        for e in range(int(getattr(msd, "n_events", 0) or 0)):
            events.append((e, _event_name(msd, e)))
    tab = getattr(frame, "ms", None)
    out = []
    if kind == "lcms":
        out.append(("chrom", "PDA chromatogram with the integrated peaks" if pda_ok else
                    "Chromatogram with the integrated peaks", True))
        if pda_ok:
            out.append(("uv", "UV spectrum of the main peak", True))
        for e, lab in events:
            out.append(("ms%d" % e, "Mass spectrum %s of the main peak" % lab, True))
        for e, lab in events:
            out.append(("tic%d" % e, "Total ion chromatogram %s" % lab, True))
        user = getattr(tab, "xics", None) or {}
        names = dict(events)
        made = [(e, x) for e in sorted(user) for x in user[e]]
        for e, x in made:
            out.append(("xic%d_%g" % (e, x["mz"]), "Extracted ion chromatogram %s, m/z %g" % (names.get(e, ""),
                                                                                            x["mz"]), True))
        if not made:
            for e, lab in events:
                out.append(("xicb%d" % e, "Extracted ion chromatogram %s of the main peak's base peak" % lab, True))
        if events:
            out.append(("shade", "Shade the time ranges on the chromatograms", True))
    elif kind == "hrms":
        out += [("chrom", "Chromatograms", True), ("spec", "Mass spectrum", True)]
        if tab is not None and getattr(getattr(tab, "dec", None), "result", None) is not None:
            out.append(("dec", "Deconvolution", True))
        if msd is not None and getattr(msd, "calib", None) is not None:
            out.append(("cal", "Calibration errors", True))
        out.append(("shade", "Shade the time ranges on the chromatograms", True))
    elif kind == "deconv":
        out += [("spec", "Mass spectrum with the fit", True), ("mass", "Zero charge mass spectrum", True)]
    elif kind == "si":
        if lc:
            out.append(("chrom", "Chromatogram", True))
            if pda_ok:
                out.append(("uv", "UV spectrum", True))
            for e, lab in events[:2]:
                out.append(("ms%d" % e, "Mass spectrum %s" % lab, True))
        else:
            out += [("chrom", "Chromatogram", True), ("spec", "Mass spectrum", True), ("iso", "Isotope pattern", True)]
            if tab is not None and getattr(getattr(tab, "dec", None), "result", None) is not None:
                out.append(("dec", "Zero charge mass spectrum", True))
    elif kind == "compare":
        ctab = _compare_tab(frame)
        guides = (getattr(ctab, "s", None) or {}).get("guides") if ctab is not None else None
        out += [("plot", "Comparison plot", True)]
        if _cmp_mz(ctab) is not None:  # the m/z tool shows the spectra of a peak picked on a trace
            out.append(("mzspec", "Spectra of the picked peak", True))
        out += [("files", "Table of the runs", True),
                ("peaks", "Peaks of each trace" + (" with area % at the guide lines" if guides else ""), True),
                ("methods", "Methods text", True)]
    return out


def want(fields, key, default=True):
    """Is the figure key chosen in the report window?"""
    figs = fields.get("figs")
    if not isinstance(figs, dict) or key not in figs:
        return default
    return bool(figs[key])


def _meta(frame, kind, sample, fields):
    import unidec_theme as T
    ver = getattr(T, "APP_VERSION", "")
    base = os.path.basename((frame.path or "").rstrip("\\/"))
    footer = base
    if kind in ("hrms", "deconv"):
        result = getattr(getattr(getattr(frame, "ms", None), "dec", None), "result", None) or {}
        method = str(result.get("method") or "")
        if method:
            from ms_brand import display
            footer += "  ·  deconvolution: " + display(method)
    return {"kind": KIND_NAMES[kind], "sample": sample, "version": ver,
            "date": datetime.datetime.now().strftime("%d.%m.%Y %H:%M"),
            "footer": footer}


def _info_common(frame, fields, si, extra=()):
    pairs = []
    sname = si.get("sample_name")
    if fields.get("compound"):
        pairs.append(("Compound", esc(fields["compound"])))
    if sname:
        pairs.append(("Sample name", esc(sname) + ((" (ID %s)" % esc(si["sample_id"])) if si.get("sample_id") else "")))
    pairs.append(("Data file", esc(os.path.basename((frame.path or "").rstrip("\\/")))))
    inst = _instrument(frame)
    if inst:
        pairs.append(("Instrument", esc(inst)))
    acq = si.get("acquired") or si.get("acquired_file")
    if acq:
        pairs.append(("Acquired", esc(fmt_date(acq)) + ("" if si.get("acquired") else " (file date)")))
    op = fields.get("operator") or si.get("operator")
    if op:
        pairs.append(("Operator", esc(op)))
    if si.get("method_file"):
        pairs.append(("Method", esc(si["method_file"])))
    if si.get("vial") or si.get("inj_vol"):
        v = []
        if si.get("vial"):
            v.append("vial %s" % esc(si["vial"]))
        if si.get("inj_vol"):
            v.append("%s µL injected" % esc(si["inj_vol"]))
        pairs.append(("Injection", ", ".join(v)))
    pairs += list(extra)
    if fields.get("notes"):
        pairs.append(("Notes", esc(fields["notes"]).replace("\n", "<br/>")))
    return pairs


def _chrom_axes(ax, tr, ranges=(), unit="Intensity", colour=None, lw=0.8):
    """One chromatogram: trace, shaded ranges, integrated peaks with their
    retention times, label at the top right."""
    t, y = np.asarray(tr["t"], float), np.asarray(tr["y"], float)
    ax.plot(t, y, color=colour or tr.get("color") or BLUE, lw=lw)
    for a, b, colr, lab in ranges:
        if b is None or abs(b - a) < 1e-9:
            ax.axvline(a, color=colr, lw=0.8, ls=(0, (3, 2)))
        else:
            ax.axvspan(a, b, color=colr, alpha=0.16, lw=0)
    for p in tr.get("peaks", []):
        i0, i1 = int(p["i0"]), int(p["i1"])
        ts, ys = t[i0:i1 + 1], y[i0:i1 + 1]
        if len(ts) > 1:
            base = p["b0"] + (p["b1"] - p["b0"]) * (ts - ts[0]) / max(ts[-1] - ts[0], 1e-9)
            ax.fill_between(ts, base, np.maximum(ys, base), color=ACCENT, alpha=0.22, lw=0)
            ax.plot([ts[0], ts[-1]], [p["b0"], p["b1"]], color=AX, lw=0.5)
            ax.annotate("%.2f" % p["rt"], (p["rt"], p.get("apex_y", ys.max())), xytext=(0, 2),
                        textcoords="offset points", ha="center", va="bottom", fontsize=6.3,
                        fontweight="bold" if p.get("_main") else "normal")
    headroom(ax, y, 0.18)
    style(ax, "", unit)
    ax.text(0.995, 0.93, tr.get("label") or tr.get("name", ""), transform=ax.transAxes, fontsize=6.8,
            va="top", ha="right")
    ax.set_xlim(float(t[0]), float(t[-1]))


H_PANEL = 5.2  # cm: every graph of the LC-MS report has this height and the full width


def _one_fig(path, draw, h=H_PANEL):
    """A full width graph of the common height; the same margins for every
    graph, so their axes line up down the page."""
    f, ax = new_fig(W_FULL, h)
    draw(ax)
    return save_fig(f, path, left=1.35 / W_FULL, right=0.975, top=1 - 0.42 / h, bottom=0.95 / h)


def _chrom_fig(path, traces, ranges=(), picks=(), w=W_FULL, per=2.5, xlabel="Time (min)", unit="Intensity"):
    """Stacked chromatogram panels (one per trace, at most 3)."""
    traces = traces[:3]
    n = max(1, len(traces))
    f, axs = new_fig(w, per * n + 1.0, 1, n, sharex=True, squeeze=False)
    axs = axs[:, 0]
    for i, (ax, tr) in enumerate(zip(axs, traces)):
        _chrom_axes(ax, tr, ranges, unit, tr.get("color") or PAL[i % len(PAL)])
    axs[-1].set_xlabel(xlabel, fontsize=7.5, labelpad=1.5)
    top = 1 - 0.55 / (per * n + 1.0)
    return save_fig(f, path, left=0.075, right=0.975, top=top, bottom=0.8 / (per * n + 1.0),
                    hspace=0.35 if n > 1 else 0.2)


def _spec_fig(path, spec, xlim=None, fmt="%.4f", n=6, tags=(), w=W_FULL, h=5.6, label=None, sep=None):
    f, ax = new_fig(w, h)
    _sticks_or_line(ax, spec)
    style(ax, "m/z", "Intensity")
    x0, x1 = xlim or (float(spec[0, 0]), float(spec[-1, 0]))
    ax.set_xlim(x0, x1)
    m = (spec[:, 0] >= x0) & (spec[:, 0] <= x1)
    headroom(ax, spec[m, 1] if np.any(m) else spec[:, 1], 0.22)
    label_peaks(ax, spec[:, 0], spec[:, 1], n=n, fmt=fmt, sep=sep or (x1 - x0) / 16.0)
    for mz, txt in tags:
        if x0 < mz < x1:
            ax.axvline(mz, color=ACCENT, lw=0.6, ls=(0, (3, 2)), alpha=0.8, zorder=0)
            ax.text(mz, 1.01, txt, transform=ax.get_xaxis_transform(), ha="center", va="bottom", fontsize=6.3,
                    color=ACCENT)
    if label:
        ax.text(0.99, 0.95, label, transform=ax.transAxes, fontsize=6.8, va="top", ha="right")
    if w < 10:
        xbins(ax, 4)
    return save_fig(f, path, left=0.075 if w > 12 else 0.2, right=0.975, top=0.88, bottom=0.17 if h > 5 else 0.22)


def _zoom_axes(ax, spec, lo, hi, pattern=None, scale=None, title=None):
    m = (spec[:, 0] >= lo) & (spec[:, 0] <= hi)
    ax.plot(spec[m, 0], spec[m, 1], color=BLUE, lw=0.6)
    ax.fill_between(spec[m, 0], 0, spec[m, 1], color=ACCENT, alpha=0.10, lw=0)
    style(ax, "m/z", "Intensity")
    ax.set_xlim(lo, hi)
    xbins(ax, 5)
    top = spec[m, 1].max() if np.any(m) else 1.0
    if pattern is not None and len(pattern):
        pat = np.asarray(pattern, float)
        s = (scale or top) / 100.0 * (100.0 / max(pat[:, 1].max(), 1e-9))
        ax.vlines(pat[:, 0], 0, pat[:, 1] * s, color=ORANGE, lw=1.1, alpha=0.9, label="calculated")
        ax.plot([], [], color=BLUE, lw=0.6, label="measured")
        ax.legend(fontsize=6, frameon=False, loc="upper right")
        top = max(top, (pat[:, 1] * s).max())
    ax.set_ylim(0, top * 1.22)
    label_peaks(ax, spec[:, 0], spec[:, 1], n=5, fmt="%.4f", sep=(hi - lo) / 7.0, xlim=(lo, hi))
    if title:
        ax.text(0.02, 0.95, title, transform=ax.transAxes, fontsize=6.8, va="top")


def _grouped(res):
    """Isotope resolved species (as deconv_tab._grouped): one entry per species
    with its average mass, most abundant isotope ("apex") and isotope peaks."""
    return any(q.get("n_iso") not in (None, "") for q in (res or {}).get("peaks") or [] if not q.get("added"))


def _mass_dec(res):
    """Decimals of the masses of a result as the window shows them (deconv_tab: _fmt_step of the
    panel and _mass_fmt): IsoDec 4, isotope resolved species 2, else by the mass step."""
    if res and res.get("tag") == "isodec":
        step = 0.0001
    elif _grouped(res):
        step = 0.1
    else:
        try:
            step = float((res or {}).get("params", {}).get("mass_step", 1) or 1)
        except (TypeError, ValueError):
            step = 1.0
    return 4 if step < 0.05 else 2 if step < 0.5 else 1


def _signed(v, dec=4):
    """A difference with its sign; a difference that rounds to zero without one (not "-0.00")."""
    s = "%+.*f" % (dec, v)
    return "%.*f" % (dec, 0.0) if float(s) == 0 else s


def _apex_of(res, q):
    """Selected isotope, or the deconvoluted modal isotope, as labelled in the UI."""
    ap = float(q.get("selected_apex", q.get("apex_dec", q.get("apex", q["mass"]))))
    iso = res.get("isotope_peaks") or []
    if iso:
        m = min((t[0] for t in iso), key=lambda v: abs(v - ap))
        if abs(m - ap) <= 0.3:
            return float(m)
    return ap


def _main_peak(res, panel=None, selected=None):
    """Use the active selected species; otherwise the strongest engine result."""
    peaks = (res or {}).get("peaks") or []
    if selected is None and panel is not None and getattr(panel, "result", None) is res:
        selected = getattr(panel, "sel", None)
    if isinstance(selected, (int, np.integer)) and 0 <= selected < len(peaks):
        return peaks[selected]
    return max(peaks, key=lambda q: q["height"]) if peaks else None


def _peak_basis(q):
    if q.get("added"):
        return "Added peak"
    return "Selected peak" if "selected_apex" in q else "Strongest peak"


def _report_isotopes(res):
    """Every detected isotope, including results that did not cache peak labels."""
    iso = res.get("isotope_peaks")
    if iso is None and _grouped(res):
        mass = res.get("mass")
        if mass is not None and len(mass):
            import ms_deconv as D
            iso = D.isotope_peaks(np.asarray(mass, float))
    return sorted([(float(m), float(h)) for m, h in (iso or [])], key=lambda t: t[0])


def _isotope_height(res, q, isotopes=None):
    ap = _apex_of(res, q)
    iso = _report_isotopes(res) if isotopes is None else isotopes
    if iso:
        m, h = min(iso, key=lambda t: abs(t[0] - ap))
        if abs(m - ap) <= 0.3:
            return h
    return float(q.get("height", 0))


def _isotope_table(rep, res, exp=None, selected=None):
    """Individual isotope masses of ALL components, independent of viewport."""
    iso = _report_isotopes(res)
    if not iso:
        return
    main = _main_peak(res, selected=selected)
    chosen = _apex_of(res, main) if main is not None else None
    hit = min(range(len(iso)), key=lambda i: abs(iso[i][0] - chosen)) if chosen is not None else None
    if hit is not None and abs(iso[hit][0] - chosen) > 0.3:
        hit = None
    top = max(h for _, h in iso) or 1.0
    rows = []
    for i, (m, h) in enumerate(iso):
        row = [str(i + 1), "%.4f" % m, g(h), "%.2f" % (100.0 * h / top),
               _peak_basis(main) if i == hit else ""]
        if exp:
            e = min(exp, key=lambda t: abs(t[0] - m))
            row += [_exp_text(e[0]), esc(e[1]), _signed(m - e[0])] if abs(e[0] - m) < 0.5 else ["", "", ""]
        rows.append(row)
    hdr = ["Peak", "Isotope mass (Da)", "Intensity", "Relative (%)", "Reported peak"]
    wid = [1.3, 3.2, 2.5, 2.4, 4.4]
    if exp:
        hdr += ["Expected (Da)", "Species", "Δ (Da)"]
        wid += [2.6, 3.0, 2.0]
    rep.add("table", hdr, rows, wid, hit)


def _iso_names(res, iso):
    """Names of the isotope peaks in the tables: M+0, M+1, ... only when they
    are the pattern of one species starting at its monoisotopic peak (a small
    molecule: every peak 1.003 Da after the one before); otherwise numbers.
    The isotope peaks of a result are those of the whole mass spectrum, so with
    several species (a protein and its adducts) they were all numbered as
    isotopes of one, and for a resolved protein the first peak seen is not its
    monoisotopic peak. Returns (names, True when named M+i)."""
    ms = [float(t[0]) for t in iso]
    one = bool(ms) and all(abs(b - a - 1.00235) < 0.12 for a, b in zip(ms, ms[1:]))
    if one and not _grouped(res):
        return ["M+%d" % i for i in range(len(ms))], True
    return ["%d" % (i + 1) for i in range(len(ms))], False


def _added(q):
    return " (added)" if q.get("added") else ""


def _species_table(rep, res, exp=None, dec=2):
    """One representative isotope per component, without envelope centroids."""
    peaks = sorted(res.get("peaks") or [], key=lambda q: _apex_of(res, q))
    iso = _report_isotopes(res)
    rows = []
    for q in peaks:
        mass = _apex_of(res, q)
        r = ["%.4f" % mass + _added(q), g(_isotope_height(res, q, iso)), str(q.get("n_iso", "") or ""),
             "%.1f" % q.get("frac", 0), _peak_basis(q)]
        if exp:
            e = min(exp, key=lambda t: abs(t[0] - mass))
            r += [_exp_text(e[0]), esc(e[1]), _signed(mass - e[0])]
        rows.append(r)
    hdr = ["Mass (Da)", "Isotope intensity", "Isotope peaks", "Species share (%)", "Peak choice"]
    wid = [3.0, 2.5, 2.2, 2.7, 4.2]
    if exp:
        hdr += ["Expected (Da)", "Species", "Δ (Da)"]
        wid = [2.8, 2.2, 1.8, 2.5, 3.3, 2.3, 3.2, 1.9]
    rep.add("table", hdr, rows, wid, None)


def _mass_fig(path, res, w=W_FULL, h=5.0, selected=None):
    mass = np.asarray(res["mass"], float)
    peaks = res.get("peaks") or []
    if mass.ndim != 2 or not len(mass):  # a method that found nothing (e.g. IsoDec)
        f, ax = new_fig(w, h)
        style(ax, "Mass (Da)", "Intensity")
        ax.text(0.5, 0.5, "no masses found", transform=ax.transAxes, ha="center", va="center", fontsize=7,
                color="#666666")
        return save_fig(f, path, left=0.075, right=0.975, top=0.9, bottom=0.2)
    zd = res.get("zdist")
    zd = np.asarray(zd, float) if zd is not None and len(zd) else None
    if zd is not None and zd.ndim == 2 and len(zd):
        f, (ax, az) = new_fig(w, h, 2, 1, gridspec_kw={"width_ratios": [2.8, 1]})
    else:
        f, ax = new_fig(w, h)
        az = None
    iso = _report_isotopes(res)
    if iso:
        ax.plot(mass[:, 0], mass[:, 1], color=BLUE, lw=0.6)
        ax.fill_between(mass[:, 0], 0, mass[:, 1], color=ACCENT, alpha=0.13, lw=0)
        xs = [m for m, _ in iso]
        lo, hi = min(xs) - 1.5, max(xs) + 1.5
        # the tallest isotope peaks that have room for their label (every one was labelled: with
        # the isotope peaks of a protein and its adducts the labels ran into each other)
        sep, done = (hi - lo) / 8.0, []
        main = _main_peak(res, selected=selected) if _grouped(res) else None
        ap = _apex_of(res, main) if main is not None else None
        chosen = min(iso, key=lambda t: abs(t[0] - ap)) if ap is not None else None
        if chosen is not None and abs(chosen[0] - ap) > 0.3:
            chosen = None
        ordered = ([chosen] if chosen is not None else []) + [t for t in sorted(iso, key=lambda t: -t[1]) if t != chosen]
        for m_, h_ in ordered:
            if len(done) >= 8:
                break
            if all(abs(m_ - a) > sep for a in done):
                done.append(m_)
                reported = chosen is not None and m_ == chosen[0]
                if reported:
                    ax.plot([m_], [h_], "o", ms=4, mfc=ORANGE, mec=INK, mew=0.6)
                ax.annotate("%.4f" % m_, (m_, h_), xytext=(0, 5 if reported else 2), textcoords="offset points", ha="center",
                            va="bottom", fontsize=7 if reported else 6.2, fontweight="bold" if reported else "normal")
    else:
        ax.plot(mass[:, 0], mass[:, 1], color=BLUE, lw=0.6)
        ax.fill_between(mass[:, 0], 0, mass[:, 1], color=ACCENT, alpha=0.13, lw=0)
        pm = [q["mass"] for q in peaks] or [float(mass[np.argmax(mass[:, 1]), 0])]
        span = max(max(pm) - min(pm), 50.0)
        lo, hi = min(pm) - 0.6 * span, max(pm) + 0.6 * span
        dec = _mass_dec(res)
        lab = sorted(peaks, key=lambda q: -q["height"])[:8]
        lab += [q for q in peaks if q.get("pinned") and all(q is not r for r in lab)]  # labelled by hand
        for q in lab:
            ax.annotate(("%%.%df" % dec) % q["mass"], (q["mass"], q["height"]), xytext=(0, 2),
                        textcoords="offset points", ha="center", va="bottom", fontsize=6.3)
    lo, hi = max(lo, float(mass[0, 0])), min(hi, float(mass[-1, 0]))
    style(ax, "Mass (Da)", "Intensity")
    ax.set_xlim(lo, hi)
    m = (mass[:, 0] >= lo) & (mass[:, 0] <= hi)
    headroom(ax, mass[m, 1] if np.any(m) else mass[:, 1], 0.2)
    if az is not None:
        az.bar(zd[:, 0], zd[:, 1], width=0.65, color="#4A74C9", edgecolor=AX, lw=0.4)
        style(az, "Charge", "")
        from matplotlib.ticker import MaxNLocator
        az.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=5))
        headroom(az, zd[:, 1], 0.15)
    return save_fig(f, path, left=0.075, right=0.975, top=0.9, bottom=0.2, wspace=0.2)


def _labelled_peaks(card, spec):
    """(m/z, intensity) of the peaks labelled on a spectrum tile (automatic labels and labels set by
    hand, at the positions written), tallest first; [] without such a tile."""
    if card is None or not hasattr(card, "ax"):
        return []
    from matplotlib.text import Annotation
    out = []
    try:
        x0, x1 = card.ax.get_xlim()
        for t in card.ax.texts:
            if not isinstance(t, Annotation) or not t.get_visible():
                continue
            try:
                m = float(t.get_text().strip())
            except ValueError:
                continue  # charge tags, measurements
            if not x0 <= t.xy[0] <= x1 or abs(m - t.xy[0]) > max(0.02, 1e-4 * abs(m)):
                continue
            out.append((float(t.xy[0]), float(t.xy[1])))  # a label sits on the top of its peak
    except Exception:
        return []
    out = sorted(set(out), key=lambda t: -t[1])
    return out[:12]


def _result_table(rep, res, selected=None):
    """The mass table of a result, as the window lists it: one row per species (isotope resolved),
    the isotope peaks (a small molecule), or the masses."""
    iso = res.get("isotope_peaks")
    if _grouped(res):  # one row per species (the isotope peaks of every species were listed as M+0, M+1, ...)
        rep.add("small", "<b>Masses</b>")
        _species_table(rep, res)
        rep.add("small", "<b>Isotope peaks</b>")
        _isotope_table(rep, res, selected=selected)
    elif iso:
        tp = max(h for _, h in iso)
        names, mono = _iso_names(res, iso)
        rows = [[names[i], "%.4f" % m, g(h), "%.1f" % (100 * h / tp)] for i, (m, h) in enumerate(iso)]
        rep.add("table", ["Isotope" if mono else "Isotope peak", "Mass (Da)", "Intensity", "Relative (%)"], rows,
                [2.2, 3.0, 2.6, 2.4], None)
    else:
        d = _mass_dec(res)
        rows = [[("%.*f" % (d, q["mass"])) + _added(q), g(q["height"]), "%.1f" % q.get("frac", 0)]
                for q in sorted(res.get("peaks", []), key=lambda q: q["mass"])]
        rep.add("table", ["Mass (Da)", "Intensity", "Share %"], rows, [3.0, 2.6, 2.0], None)


def _shift_blocks(rep, tab, ent, heading=None):
    """The mass shifts of a result as its row shows them (the table of the pairs and, with a tag, the
    degree of conjugation), when they are shown for it (deconv_shifts.ShiftPanel.report_tables); heading:
    a section title before them (None: a bold title in the text)."""
    sp = getattr(getattr(tab, "dec", None), "shifts", None)
    try:
        t = sp.report_tables(ent) if sp is not None and ent is not None else None
    except Exception as ex:
        print("report: the mass shifts were left out (%s)" % ex)
        t = None
    if not t:
        return False
    if heading:
        rep.add("section", heading)
    hdr, rows, wid = t["pairs"]
    if rows:
        if not heading:
            rep.add("small", "<b>Mass shifts</b>")
        rep.add("table", [esc(h) for h in hdr], [[esc(c) for c in r] for r in rows], wid, None)
    if t.get("conj"):
        hdr, rows, wid = t["conj"]
        rep.add("small", "<b>Degree of conjugation</b>")
        rep.add("table", [esc(h) for h in hdr], [[esc(c) for c in r] for r in rows], wid, None)
    return True


def _result_figures(rep, tab, ent, tmp, tag, h=6.0):
    """The figure(s) of one result as its row shows it (zero charge mass spectrum, and the charge
    states when that tile is shown); returns the image description of DeconvPanel.report_images,
    or None (no row on screen)."""
    dec = getattr(tab, "dec", None)
    try:
        im = dec.report_images(ent, tmp, tag, W_FULL, h, TEXT_SCALE)
    except Exception as ex:
        print("report: no image of the result row (%s); the report draws it" % ex)
        return None
    if im.get("z"):
        rep.add("pair", [("figure", im["mass"], im["w_mass"])], [("figure", im["z"], im["w_z"] - 0.1)],
                im["w_mass"])
    else:
        rep.add("figure", im["mass"], W_FULL)
    return im


def _result_blocks(rep, tab, ents, tmp, tag, section):
    """A section with every result given (figures as on screen and the mass table); returns the
    results drawn."""
    done = []
    for i, ent in enumerate(ents):
        if not done:
            rep.add("section", section)
        im = _result_figures(rep, tab, ent, tmp, "%s%d" % (tag, i))
        if im is None:
            if not done:
                rep.blocks.pop()  # the heading without a figure
            continue
        _result_table(rep, ent.res, selected=getattr(ent, "sel", None))
        _shift_blocks(rep, tab, ent)
        done.append(ent)
    return done


def _calib_rows(cal):
    """Every model on the calibrants of the calibration applied: the row of
    the model in use is that calibration itself (with its HPC order), the
    other HPC row uses the same order, as the Models table of the window."""
    import hrms_calib as HC
    out = []
    order = getattr(cal, "hpc_order", None) if cal.model == 5 else None
    for m in range(6):
        try:
            c = cal if m == cal.model else HC.Calibration(m, cal.measured, cal.reference,
                                                          hpc_order=order if m == 5 else None)
            nm = c.name().replace(" (quadratic + order ", " (order ").replace(" correction)", ")")
            out.append((m, nm, c.rms, c.cv_rms))
        except Exception:
            continue
    return out


def _hrms_state(frame):
    tab, data = frame.ms, frame.ms_data
    col = _spec_col(tab)
    e = col["e"]
    spec = np.asarray(col["spec"], float)
    res = getattr(tab.dec, "result", None)
    dec = res if (res is not None and res.get("e") == e) else None
    # the calibration of this spectrum (one made on the other polarity is not used for it)
    cal = data.calibration_for(e) if hasattr(data, "calibration_for") else getattr(data, "calib", None)
    ranges = []
    evs = getattr(getattr(tab, "cal", None), "events_of", None)
    if cal is not None and getattr(tab, "cal", None) is not None and getattr(tab.cal, "range", None) and \
            (evs is None or e in evs(getattr(tab.cal, "event", e))):  # the range of this polarity's calibration
        a, b = tab.cal.range
        ranges.append((a, b, GREEN, "calibrant"))
    if tab.avg:
        ranges.append((tab.avg[0], tab.avg[1], ACCENT, "averaged"))
    elif tab.pick_t is not None:
        ranges.append((tab.pick_t, None, ACCENT, "scan"))
    return {"tab": tab, "data": data, "col": col, "e": e, "spec": spec, "dec": dec, "cal": cal,
            "fi": _formula_info(col), "ranges": ranges, "traces": _traces(tab, e)}


def build_hrms(frame, fields, tmp):
    st = _hrms_state(frame)
    si = sample_info(frame)
    spec, fi, dec, cal = st["spec"], st["fi"], st["dec"], st["cal"]
    name = fields.get("compound") or si.get("sample_name") or frame.base_name()
    rep = Report("hrms", _meta(frame, "hrms", name, fields))
    sub = [esc(name)] + ([esc(fi["name"])] if fi else [])
    rep.add("title", "HRMS compound report", "  ·  ".join(sub))
    recal = getattr(st["data"], "recal_note", "")
    extra = [("Calibration", "Yes" if cal is not None or recal else "No")]
    if st["tab"].avg:
        extra.append(("Time range (min)", "%.2f to %.2f" % tuple(st["tab"].avg[:2])))
    elif st["tab"].pick_t is not None:
        extra.append(("Scan time (min)", "%.2f" % st["tab"].pick_t))
    if dec and dec.get("method"):
        extra.append(("Deconvolution", esc(dec["method"])))
    rep.add("info", _info_common(frame, fields, si, extra))
    top = local_max(spec[:, 0], spec[:, 1], n=8, sep=0.5)
    items = []
    if dec and dec.get("peaks"):
        q = _main_peak(dec, st["tab"].dec)
        if _grouped(dec):
            items.append(("%.4f Da" % _apex_of(dec, q), "mass: " + _peak_basis(q).lower()))
        else:
            items.append(("%.*f Da" % (_mass_dec(dec), q["mass"]),  # decimals as the window shows them
                          "monoisotopic mass (IsoDec)" if dec.get("tag") == "isodec" else "deconvoluted mass"))
    if fi and fi.get("found") is not None and fi.get("ok") is False:  # the nearest peak is outside the tolerance
        items += [("outside tolerance", "m/z found, %s" % fi.get("ion", "")),
                  ("%.4f" % fi["found"] + (", %+.2f ppm" % fi["err"] if fi.get("err") is not None else ""),
                   "nearest peak vs calculated %.4f" % fi["mz"]),
                  ("%.0f %%" % fi["match"] if fi.get("match") is not None else "–", "isotope pattern match")]
    elif fi and fi.get("found") is not None:
        items += [("%.4f" % fi["found"], "m/z found, %s" % fi.get("ion", "")),
                  ("%+.2f ppm" % fi["err"], "error vs calculated %.4f" % fi["mz"]),
                  ("%.0f %%" % fi["match"], "isotope pattern match")]
    else:
        items += [("%.4f" % top[0][0] if top else "–", "base peak, m/z"),
                  (g(top[0][1]) if top else "–", "base peak intensity")]
    items.append(("%.2f ppm" % cal.rms if cal is not None else ("file recal." if recal else "as recorded"),
                  "calibration"))
    rep.add("result", items)

    mark = want(fields, "shade", True)
    if want(fields, "chrom"):
        # the chromatogram tiles as HRMS Analysis shows them (zoom, traces, mass chromatograms, shading)
        imgs = [p for p in (_tile_png(getattr(v, "card", None), os.path.join(tmp, "a_chrom%d.png" % i), W_FULL, 4.4,
                                      **_mark_opts(mark))
                            for i, v in enumerate(getattr(st["tab"], "views", []) or []))
                if p]
        if imgs:
            rep.add("section", "Chromatograms")
            for p in imgs:
                rep.add("figure", p, W_FULL)
        elif st["traces"]:
            rep.add("section", "Chromatograms")
            rep.add("figure", _chrom_fig(os.path.join(tmp, "a_chrom.png"), st["traces"],
                                         st["ranges"] if mark else ()), W_FULL)
    ev = st["data"].events[st["e"]] if hasattr(st["data"], "events") else {}
    lo = ev.get("mz_low") or float(spec[0, 0])
    hi = ev.get("mz_high") or float(spec[-1, 0])
    tags = []
    if fi and fi.get("mz"):
        tags.append((fi["mz"], fi.get("ion", "")))
    if want(fields, "spec"):
        rep.add("section", "Mass spectrum")
        # the spectrum tiles as shown (zoom, labels, formula pattern, the fit and charge states of the
        # deconvolution shown in the table); the spectrum of this report first
        cols = [st["col"]] + [c for c in (getattr(st["tab"], "cols", None) or [])
                              if c is not st["col"] and c.get("spec") is not None and len(c["spec"])]
        shown = []
        for i, c in enumerate(cols):
            p = _tile_png(c.get("spec_card"), os.path.join(tmp, "a_spec%d.png" % i), W_FULL, 5.8)
            if p:
                shown.append((p, c))
        if shown:
            for p, c in shown:
                rep.add("figure", p, W_FULL)
        else:
            rep.add("figure", _spec_fig(os.path.join(tmp, "a_spec.png"), spec, (lo, hi), "%.4f", 6, tags,
                                        sep=(hi - lo) / 14.0), W_FULL)
    # The duplicate isotope-pattern zoom panel is intentionally omitted,
    # including when an older saved report selection still has iso=True.
    rep.add("section", "Exact mass")
    if fi:
        out = fi.get("found") is not None and fi.get("ok") is False  # nearest peak outside the tolerance
        rows = [[esc(fi.get("name", "")), str(fi.get("z", "")), "%.5f" % fi["mz"],
                 "outside tolerance" if out else "%.5f" % fi["found"] if fi.get("found") is not None else "not found",
                 "%+.2f" % fi["err"] if fi.get("err") is not None else "–",
                 "%.0f %%" % fi["match"] if fi.get("match") is not None else "–"]]
        rep.add("table", ["Ion", "z", "m/z calcd", "m/z found", "Δ (ppm)", "Isotope match"], rows,
                [5.4, 1.2, 2.8, 2.8, 2.2, 3.4], None)
    else:
        shown = _labelled_peaks(st["col"].get("spec_card"), spec)  # the peaks labelled on the spectrum tile
        lst = shown or top
        tt = max(t for _, t in lst) if lst else 1.0
        rows = [["%.4f" % m, g(h), "%.1f" % (100 * h / tt)] for m, h in lst]
        rep.add("table", ["m/z", "Intensity", "Relative (%)"], rows, [3.4, 3.0, 3.0], None)
    if want(fields, "dec"):
        # every result drawn on screen (visible rows, newest first), each as its row shows it
        ents = st["tab"].dec.report_results() if hasattr(getattr(st["tab"], "dec", None), "report_results") else []
        drawn = _result_blocks(rep, st["tab"], ents, tmp, "a_dec", "Deconvolution")
        if not drawn and dec:
            rep.add("section", "Deconvolution")
            rep.add("figure", _mass_fig(os.path.join(tmp, "a_dec.png"), dec, selected=getattr(st["tab"].dec, "sel", None)), W_FULL)
            _result_table(rep, dec, selected=getattr(st["tab"].dec, "sel", None))
    if cal is not None and want(fields, "cal"):
        rep.add("section", "Calibration")
        f, ax = new_fig(8.6, 4.8)
        ax.axhline(0, color=AX, lw=0.5)
        ax.plot(cal.reference, cal.before, "o", ms=2.8, mfc="white", mec=GREY, mew=0.7, label="as recorded")
        ax.plot(cal.reference, cal.after, "o", ms=2.8, color=ACCENT, label="calibrated")
        style(ax, "m/z", "Error (ppm)")
        from matplotlib.ticker import ScalarFormatter
        ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.legend(fontsize=6.3, frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2,
                  handletextpad=0.2, columnspacing=1.2)
        figp = save_fig(f, os.path.join(tmp, "a_cal.png"), left=0.16, right=0.97, top=0.86, bottom=0.22)
        rows, hl = [], None
        for i, (m, nm, rms, cv) in enumerate(_calib_rows(cal)):
            rows.append([esc(nm), "%.2f" % rms, "%.2f" % cv if cv is not None else "–"])
            if m == cal.model:
                hl = i
        rep.add("pair", [("figure", figp, 8.6)],
                [("table", ["Model", "RMS (ppm)", "CV (ppm)"], rows, [4.6, 1.9, 1.9], hl)], 8.8)
    return _plain_report(rep, fields)


def _lc_state(frame):
    pda, pdata, msd = frame.pda, frame.pda_data, frame.ms_data
    st = {"auto": False}
    tr = None
    if pdata is not None and pda.view.traces:
        tr = pda.view.traces[0]
        tab_peaks = pda.peaks
        st.update(src="pda", unit="Absorbance (mAU)", wl=pda.wl_sel, bw=pda.bw.GetValue().strip() or "4")
    elif msd is not None and frame.ms.views and frame.ms.views[0].traces:
        tr = frame.ms.views[0].traces[0]
        tab_peaks = frame.ms.peaks
        st.update(src="ms", unit="Intensity", wl=None, bw=None)
    if tr is None:
        raise RuntimeError("no chromatogram to report")
    t, y = np.asarray(tr["t"], float), np.asarray(tr["y"], float)
    keys = {tr.get("key"), tr.get("name")}
    peaks = [dict(p) for p in tab_peaks if p.get("key") in keys or p.get("trace") in keys]
    if not peaks:
        import lcms_integrate as LI
        peaks = LI.finish(LI.auto_integrate(t, y, threshold_pct=1.0, min_width_s=2.0))
        st["auto"] = True
    if not peaks:
        raise RuntimeError("no peak found in the chromatogram (integrate it first)")
    main = max(peaks, key=lambda p: p["area"])
    main["_main"] = True
    label = (tr.get("label") or tr.get("name") or "").replace(" (slider)", "").strip()
    st.update(t=t, y=y, label=label, peaks=peaks, main=main,
              total=sum(p["area"] for p in peaks))
    if pdata is not None:
        # the UV spectrum as shown in LCMS Analysis (PDA view), else at the apex
        shown_uv = getattr(pda, "uv", None)
        if shown_uv is not None and len(shown_uv[0]):
            st["uv"] = np.column_stack([np.asarray(shown_uv[0], float), np.asarray(shown_uv[1], float)])
            st["uv_desc"] = getattr(pda, "uv_desc", "") or ""
            st["uv_src"] = "window"
        else:
            wl, a = pdata.spectrum(main["rt"])
            st["uv"] = np.column_stack([wl, a])
            st["uv_desc"] = "at the apex, %.2f min" % main["rt"]
            st["uv_src"] = "peak"
    delay = 0.0
    try:
        delay = float(pda.delay.GetValue().replace(",", ".")) if st["src"] == "pda" else 0.0
    except (ValueError, AttributeError):
        pass
    st["delay"] = delay
    ms = []
    if msd is not None:
        # the mass spectra as shown in LCMS Analysis (same time range,
        # background and bin width); an event without a spectrum in the
        # window: averaged over the main peak with the window's settings
        tab = getattr(frame, "ms", None)
        shown = {}
        for c in (getattr(tab, "cols", None) or []):
            if c.get("spec") is not None and len(c["spec"]):
                shown[c["e"]] = c
        bg, binw = None, 0.05
        try:
            bg = tab.bg if (tab.use_bg.GetValue() and tab.bg) else None
            v = float(tab.binw.GetValue().strip().replace(",", "."))
            binw = v if v > 0 else 0.05
        except Exception:
            pass
        for e in range(msd.n_events):
            c = shown.get(e)
            if c is not None:
                s = np.asarray(c["spec"], float)
                desc, src = c.get("desc") or "", "window"
            else:
                s, n = msd.average(e, main["t0"] + delay, main["t1"] + delay, bg=bg, binw=binw)
                s = np.asarray(s, float)
                desc = "%.2f to %.2f min, %d scans, averaged over the main peak%s" % (
                    main["t0"] + delay, main["t1"] + delay, n, ", background %.2f to %.2f min subtracted" % bg
                    if bg else "")
                src = "peak"
            pol = (msd.events[e].get("polarity") if hasattr(msd, "events") else None) or "+"
            lab = "ESI(+)" if pol == "+" else "ESI(−)" if pol == "-" else msd.event_label(e)
            tt, yy = msd.chromatogram(e, "tic")
            ms.append({"e": e, "label": lab, "spec": s, "desc": desc, "src": src,
                       "tic": (np.asarray(tt), np.asarray(yy)),
                       "tops": local_max(s[:, 0], s[:, 1], n=5, sep=1.5) if len(s) else []})
    st["ms"] = ms
    st["xics"] = _lc_xics(frame, ms)
    return st


def _lc_xics(frame, ms):
    """Extracted ion chromatograms: those made in LCMS Analysis (Mass chrom.),
    else one of the base peak of the main peak per polarity. Each: key (for
    the figure choice), label, t, y."""
    msd, tab = frame.ms_data, getattr(frame, "ms", None)
    out = []
    if msd is None:
        return out
    import unilcms as U
    try:
        n = int(tab.smooth.GetValue()) if tab is not None else 0
    except Exception:
        n = 0
    names = {m["e"]: m["label"] for m in ms}
    user = getattr(tab, "xics", None) or {}
    for e in sorted(user):
        for x in user[e]:
            try:
                t, y = msd.chromatogram(e, "xic", x["mz"], U.xic_tol(x))
                y = U.LI.smooth(np.asarray(y, float), n) if n else np.asarray(y, float)
            except Exception:
                continue
            out.append({"key": "xic%d_%g" % (e, x["mz"]), "e": e, "t": np.asarray(t, float), "y": y,
                        "label": "%s  m/z %g %s" % (names.get(e, ""), x["mz"], U.xic_window_text(x)),
                        "mz": x["mz"], "user": True})
    if out:
        return out
    try:
        w, ppm = tab.default_window() if tab is not None else (0.5, False)
    except Exception:
        w, ppm = 0.5, False
    for m in ms:
        if not m["tops"]:
            continue
        mz = float(m["tops"][0][0])
        tol = mz * w * 1e-6 if ppm else w
        try:
            t, y = msd.chromatogram(m["e"], "xic", mz, tol)
        except Exception:
            continue
        out.append({"key": "xicb%d" % m["e"], "e": m["e"], "t": np.asarray(t, float), "y": np.asarray(y, float),
                    "label": "%s  m/z %.1f \u00b1%g%s" % (m["label"], mz, w, " ppm" if ppm else ""), "mz": mz,
                    "user": False})
    return out


def build_lcms(frame, fields, tmp):
    st = _lc_state(frame)
    si = sample_info(frame)
    name = fields.get("compound") or si.get("sample_name") or frame.base_name()
    main, peaks, total = st["main"], st["peaks"], st["total"] or 1.0
    rep = Report("lcms", _meta(frame, "lcms", name, fields))
    rep.add("title", "LC-MS purity report", esc(name) + ("  ·  " + esc(si["sample_id"]) if si.get("sample_id")
                                                         else ""))
    det = ("%s nm, bandwidth %s nm" % (("%.0f" % st["wl"]) if st["wl"] else "?", st["bw"])) if st["src"] == "pda" \
        else "MS: " + esc(st["label"])
    extra = [("Detection", det), ("Integration", "automatic (threshold 1 %, width 2 s)"
              if st["auto"] else "LCMS Analysis, %d peaks" % len(peaks))]
    rep.add("info", _info_common(frame, fields, si, extra))
    items = [("%.2f min" % main["rt"], "retention time, main peak"),
             ("%.1f %%" % (100 * main["area"] / total), "area %s (%d peaks)" % (
                 "at %.0f nm" % st["wl"] if st["wl"] else "", len(peaks)))]
    for m in st["ms"][:2]:
        if m["tops"]:
            items.append(("%.1f" % m["tops"][0][0], "base peak %s, m/z" % m["label"]))
    rep.add("result", items)
    rep.add("section", ("PDA chromatogram, %.0f nm" % st["wl"]) if st["src"] == "pda" else "Chromatogram")
    mark = want(fields, "shade", True)
    chrom_tile = None
    if want(fields, "chrom") and not st["auto"]:
        # the chromatogram tile as shown (zoom, integrated peaks, labels); with peaks integrated
        # only for this report the report draws them
        card = getattr(frame.pda, "chrom", None) if st["src"] == "pda" else \
            getattr(frame.ms.views[0], "card", None) if frame.ms.views else None
        chrom_tile = _tile_png(card, os.path.join(tmp, "b_chrom.png"), W_FULL, H_PANEL, **_mark_opts(mark))
        if chrom_tile:
            rep.add("figure", chrom_tile, W_FULL)
    if want(fields, "chrom") and chrom_tile is None:
        tr = {"t": st["t"], "y": st["y"], "label": st["label"], "peaks": peaks, "color": BLUE}

        def draw_chrom(ax):
            _chrom_axes(ax, tr, (), st["unit"], BLUE)
            ax.set_xlabel("Time (min)", fontsize=7.5, labelpad=1.5)
        rep.add("figure", _one_fig(os.path.join(tmp, "b_chrom.png"), draw_chrom), W_FULL)
    rows = [[str(i + 1), "%.3f" % p["rt"], "%.3f" % p["t0"], "%.3f" % p["t1"], "%.1f" % p["height"],
             "%.1f" % p["area"], "%.2f" % (100 * p["area"] / total)] for i, p in enumerate(peaks)]
    unit = "mAU" if st["src"] == "pda" else "counts"
    rep.add("table", ["#", "RT (min)", "Start", "End", "Height (%s)" % unit, "Area (%s·s)" % unit, "Area %"],
            rows, [1.0, 2.4, 2.2, 2.2, 3.0, 3.2, 2.0], peaks.index(main))
    rng = [(main["t0"] + st["delay"], main["t1"] + st["delay"], ACCENT, "")] if mark else []
    ms_tab = getattr(frame, "ms", None)
    main_views = {v.key: v for v in (getattr(ms_tab, "views", None) or []) if getattr(v, "kind", "main") == "main"}

    def chrom_drawer(tr_):
        def draw(ax):
            _chrom_axes(ax, tr_, rng, "Intensity", BLUE)
            ax.set_xlabel("Time (min)", fontsize=7.5, labelpad=1.5)
        return draw
    # total ion chromatograms, each a graph of its own (same size as the others)
    tics = [m for m in st["ms"] if want(fields, "tic%d" % m["e"])]
    if tics:
        rep.add("section", "Total ion chromatograms" if len(tics) > 1 else "Total ion chromatogram")
        for m in tics:
            v = main_views.get(m["e"])
            p = _tile_png(getattr(v, "card", None), os.path.join(tmp, "b_tic%d.png" % m["e"]), W_FULL, H_PANEL,
                          **_mark_opts(mark))
            if not p:
                tr_ = {"t": m["tic"][0], "y": m["tic"][1], "label": "TIC " + m["label"]}
                p = _one_fig(os.path.join(tmp, "b_tic%d.png" % m["e"]), chrom_drawer(tr_))
            rep.add("figure", p, W_FULL)
    # the main peak: UV and mass spectra, each a graph of its own
    rep.add("section", "Main peak, %.2f min" % main["rt"])
    shown = []
    if st.get("uv") is not None and want(fields, "uv"):
        p = _tile_png(getattr(frame.pda, "spec_card", None), os.path.join(tmp, "b_uv.png"), W_FULL, H_PANEL) \
            if st.get("uv_src") == "window" else None
        rep.add("figure", p or _one_fig(os.path.join(tmp, "b_uv.png"),
                                        lambda ax: _panel(ax, "uv", st["uv"], "UV", nb=8)), W_FULL)
        shown.append("uv")
    spec_cards = {c["e"]: c.get("spec_card") for c in (getattr(ms_tab, "cols", None) or [])}
    for m in st["ms"]:
        if len(m["spec"]) and want(fields, "ms%d" % m["e"]):
            p = _tile_png(spec_cards.get(m["e"]), os.path.join(tmp, "b_ms%d.png" % m["e"]), W_FULL, H_PANEL) \
                if m["src"] == "window" else None
            rep.add("figure", p or _one_fig(os.path.join(tmp, "b_ms%d.png" % m["e"]),
                                            lambda ax, m=m: _panel(ax, "ms", m["spec"], m["label"], nb=8)), W_FULL)
            shown.append(m)
    if shown:
        what = []
        if "uv" in shown:
            what.append("UV spectrum %s" % esc(st.get("uv_desc", "")))
        for m in shown:
            if m != "uv":
                what.append("%s mass spectrum %s" % (m["label"], esc(m["desc"])))
        txt = "; ".join(what) + "."
        rep.add("caption", txt[:1].upper() + txt[1:])
    if st["ms"]:
        cols = st["ms"][:2]
        n = max(len(m["tops"]) for m in cols)
        rows = []
        for i in range(n):
            r = []
            for m in cols:
                r += (["%.2f" % m["tops"][i][0], g(m["tops"][i][1])] if i < len(m["tops"]) else ["", ""])
            rows.append(r)
        hdr = []
        for m in cols:
            hdr += ["m/z %s" % m["label"], "Intensity"]
        rep.add("table", hdr, rows, [3.0, 2.6] * len(cols), None)
        if fields.get("assignment"):
            rep.add("text", "<b>Assignment:</b> " + esc(fields["assignment"]))
    # extracted ion chromatograms
    xics = [x for x in st.get("xics", []) if want(fields, x["key"])]
    if xics:
        rep.add("section", "Extracted ion chromatograms" if len(xics) > 1 else "Extracted ion chromatogram")
        xic_cards = {}
        for v in (getattr(ms_tab, "xic_views", None) or []):  # mass chromatograms in tiles of their own
            if getattr(v, "xic", None):
                xic_cards["xic%d_%g" % (v.key, v.xic["mz"])] = v.card
        for k, x in enumerate(xics):
            p = _tile_png(xic_cards.get(x["key"]), os.path.join(tmp, "b_xic%d.png" % k), W_FULL, H_PANEL,
                          **_mark_opts(mark))
            rep.add("figure", p or _one_fig(os.path.join(tmp, "b_xic%d.png" % k), chrom_drawer(x)), W_FULL)
    return rep


def _panel(ax, kind, s, lab, nb=4):
    s = np.asarray(s, float)
    if s.ndim != 2 or not len(s):
        style(ax, "Mass (Da)" if kind == "mass" else "", "")
        ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center", va="center", fontsize=7, color="#666666")
        return
    s = np.asarray(s, float)
    if kind == "uv":
        ax.plot(s[:, 0], s[:, 1], color=BLUE, lw=0.6)
        style(ax, "Wavelength (nm)", "Absorbance (mAU)")
        lo = max(float(s[0, 0]), 200.0)
        ax.set_xlim(lo, float(s[-1, 0]))
        m = s[:, 0] >= lo
        headroom(ax, s[m, 1], 0.2)
        label_peaks(ax, s[m, 0], s[m, 1], n=3, fmt="%.0f", sep=25, xlim=(lo, float(s[-1, 0])), min_rel=0.1)
    elif kind == "chrom":
        ax.plot(s[:, 0], s[:, 1], color=BLUE, lw=0.6)
        style(ax, "Time (min)", "Intensity")
        ax.set_xlim(float(s[0, 0]), float(s[-1, 0]))
        headroom(ax, s[:, 1], 0.15)
    elif kind == "mass":
        ax.plot(s[:, 0], s[:, 1], color=BLUE, lw=0.6)
        ax.fill_between(s[:, 0], 0, s[:, 1], color=ACCENT, alpha=0.13, lw=0)
        style(ax, "Mass (Da)", "Intensity")
        sig = np.where(s[:, 1] > 0.02 * s[:, 1].max())[0]
        if len(sig):
            a, b = float(s[sig[0], 0]), float(s[sig[-1], 0])
            pad = max(0.25 * (b - a), 1.0)
            ax.set_xlim(max(float(s[0, 0]), a - pad), min(float(s[-1, 0]), b + pad))
        headroom(ax, s[:, 1], 0.2)
    else:
        _sticks_or_line(ax, s)
        style(ax, "m/z", "Intensity")
        lo, hi = float(s[0, 0]), float(s[-1, 0])
        ax.set_xlim(lo, hi)
        headroom(ax, s[:, 1], 0.22)
        label_peaks(ax, s[:, 0], s[:, 1], n=7 if nb >= 8 else 4, fmt="%.1f" if (hi - lo) > 300 else "%.4f",
                    sep=(hi - lo) / (30.0 if nb >= 8 else 20.0), min_rel=0.05 if nb >= 8 else 0.08)
    xbins(ax, nb)
    ax.text(0.98, 0.95, lab, transform=ax.transAxes, fontsize=6.8, va="top", ha="right")


_EXP_NUM = r"([0-9]+(?:,[0-9]+)+\.[0-9]+|[0-9]+(?:[.,][0-9]+)?)"


def _exp_float(txt):
    """A number as typed: commas before a decimal point group the thousands (148,056.3), a comma
    without a point after it is the decimal mark (1234,5)."""
    return float(txt.replace(",", "") if "." in txt else txt.replace(",", "."))


def _parse_expected(text):
    out = []
    for line in (text or "").splitlines():
        m = re.match(r"\s*%s\s*$" % _EXP_NUM, line)  # a number alone first ("1234,5" is not 1234 and "5")
        if m:
            out.append((_exp_float(m.group(1)), ""))
        else:
            m = re.match(r"\s*%s\s*[,;:\t ]\s*(.*)$" % _EXP_NUM, line)
            if m:
                out.append((_exp_float(m.group(1)), m.group(2).strip()))
    return out


def _exp_text(v):
    """An expected mass with all the decimals it was typed with (not rounded to 6 digits as by %g)."""
    s = repr(float(v))
    return s[:-2] if s.endswith(".0") else s


def build_deconv(frame, fields, tmp):
    res = frame.ms.dec.result
    si = sample_info(frame)
    name = fields.get("compound") or si.get("sample_name") or frame.base_name()
    p = res.get("params", {})
    peaks = sorted(res.get("peaks", []), key=lambda q: -q["height"])
    rep = Report("deconv", _meta(frame, "deconv", name, fields))
    rep.add("title", "Deconvolution report", esc(name) + "  ·  " + esc(res.get("method", "")))
    extra = [("Method", esc(res.get("method", "")))]
    rep.add("info", _info_common(frame, fields, si, extra))
    step = float(p.get("mass_step", 1) or 1)
    dec = _mass_dec(res)  # as the window shows the masses
    mf = "%%.%df Da" % dec
    items = []
    iso = res.get("isotope_peaks")
    grouped = _grouped(res)
    names, mono = _iso_names(res, iso) if iso else ([], False)
    if grouped:
        main = _main_peak(res, getattr(frame.ms, "dec", None))
        items.append(("%.4f Da" % _apex_of(res, main), "mass: " + _peak_basis(main).lower()))
        items.append((g(_isotope_height(res, main)), "peak intensity"))
    elif iso and mono:
        items.append(("%.4f Da" % iso[0][0], "first isotope peak (M+0)"))
        top_iso = max(iso, key=lambda t: t[1])
        items.append(("%.4f Da" % top_iso[0], "most abundant isotope"))
    else:
        for q in peaks[:2]:
            items.append((mf % q["mass"], "mass, %.0f %% of the total" % q.get("frac", 0)))
    z = p.get("z_range")
    if z:
        items.append(("%d to %d" % tuple(z), "charge states"))
    if res.get("r2") is not None:
        items.append(("%.3f" % res["r2"], "fit R²"))
    rep.add("result", items)
    data = res.get("input_data")
    if data is None:
        data = res.get("data")
    data = np.asarray(data, float) if data is not None else None
    tab = frame.ms
    dpanel = getattr(tab, "dec", None)
    ent = getattr(dpanel, "active", None)
    spec_tile = None
    if want(fields, "spec"):
        # the spectrum tile as shown: its zoom and labels, the fit and the charge states of the mass
        # selected in the table
        col = next((c for c in (getattr(tab, "cols", None) or []) if c.get("e") == res.get("e")), None)
        spec_tile = _tile_png(col.get("spec_card") if col else None, os.path.join(tmp, "c_spec.png"), W_FULL, 5.8)
        if spec_tile:
            rep.add("section", "Mass spectrum and fit")
            rep.add("figure", spec_tile, W_FULL)
    if spec_tile is None and data is not None and len(data) and want(fields, "spec"):
        rep.add("section", "Mass spectrum and fit")
        f, ax = new_fig(W_FULL, 5.4)
        _sticks_or_line(ax, data, 0.7)
        fit = res.get("fit")
        if fit is not None and len(fit):
            fit = np.asarray(fit, float)
            ax.plot(fit[:, 0], fit[:, 1], color=ORANGE, lw=0.7, alpha=0.9, label="fit")
            ax.plot([], [], color=BLUE, lw=0.7, label="spectrum")
            ax.legend(fontsize=6.3, frameon=False, loc="upper right", bbox_to_anchor=(1, 0.97))
        style(ax, "m/z", "Intensity")
        mzr = p.get("mz_range") or (float(data[0, 0]), float(data[-1, 0]))
        ax.set_xlim(*mzr)
        m = (data[:, 0] >= mzr[0]) & (data[:, 0] <= mzr[1])
        headroom(ax, data[m, 1] if np.any(m) else data[:, 1], 0.14)
        if peaks and z:
            M = peaks[0]["mass"]
            sign = p.get("sign", 1)
            import ms_deconv as _D
            adm = _D.carrier_mass(p)
            last = None  # high charges crowd at low m/z: a tag only where it has room (they overlapped)
            for zz in range(int(z[1]), int(z[0]) - 1, -1):
                mz = (M + zz * adm) / zz
                if mzr[0] < mz < mzr[1] and (last is None or mz - last >= (mzr[1] - mzr[0]) / 30.0):
                    inside = (mz - mzr[0]) / ((mzr[1] - mzr[0]) or 1.0) < 0.1  # where the axis exponent sits
                    ax.text(mz, 0.97 if inside else 1.01, "%d%s" % (zz, "+" if sign > 0 else "−"),
                            transform=ax.get_xaxis_transform(), ha="center", va="top" if inside else "bottom",
                            fontsize=6, color=ACCENT)
                    last = mz
        rep.add("figure", save_fig(f, os.path.join(tmp, "c_spec.png"), left=0.075, right=0.975, top=0.88,
                                   bottom=0.17), W_FULL)
    if want(fields, "mass"):
        rep.add("section", "Deconvoluted mass spectrum")
        im = _result_figures(rep, tab, ent, tmp, "c_mass") if (ent is not None and ent.res is res) else None
        if im is None:
            rep.add("figure", _mass_fig(os.path.join(tmp, "c_mass.png"), res, selected=getattr(dpanel, "sel", None)), W_FULL)
    rep.add("section", "Masses")
    exp = _parse_expected(fields.get("expected", ""))
    peaks_for_table = peaks
    if grouped:
        rep.add("small", "<b>Masses</b>")
        _species_table(rep, res, exp, max(dec, 2))
        rep.add("small", "<b>Isotope peaks</b>")
        _isotope_table(rep, res, exp, selected=getattr(dpanel, "sel", None))
        peaks_for_table = []
        iso = None
    if iso:
        tp = max(h for _, h in iso)
        irows = []
        for i, (m_, h_) in enumerate(iso):
            r = [names[i], "%.4f" % m_, g(h_), "%.1f" % (100 * h_ / tp)]
            if exp:
                e = min(exp, key=lambda t: abs(t[0] - m_))
                r += [_exp_text(e[0]), esc(e[1]), _signed(m_ - e[0])] if abs(e[0] - m_) < 0.5 else ["", "", ""]
            irows.append(r)
        hdr = ["Isotope" if mono else "Isotope peak", "Mass (Da)", "Intensity", "Relative (%)"]
        wid = [2.0, 2.8, 2.4, 2.2]
        if exp:
            hdr += ["Expected (Da)", "Species", "Δ (Da)"]
            wid += [2.6, 4.0, 2.0]
        rep.add("table", hdr, irows, wid, None)
        exp = []
        peaks_for_table = []
    rows = []
    for q in sorted(peaks_for_table, key=lambda q: q["mass"]):
        r = [mf.replace(" Da", "") % q["mass"] + _added(q), g(q["height"]), "%.1f" % q.get("frac", 0)]
        if exp:
            e = min(exp, key=lambda t: abs(t[0] - q["mass"]))
            r += [_exp_text(e[0]), esc(e[1]), _signed(q["mass"] - e[0], dec)]
        rows.append(r)
    if rows and exp:
        rep.add("table", ["Mass found (Da)", "Intensity", "Share %", "Expected (Da)", "Species", "Δ (Da)"], rows,
                [3.0, 2.4, 1.8, 2.6, 5.4, 2.2], None)
    elif rows:
        rep.add("table", ["Mass found (Da)", "Intensity", "Share %"], rows, [3.0, 2.6, 2.0], None)
    if ent is not None and ent.res is res:
        _shift_blocks(rep, tab, ent, "Mass shifts")
    rep.add("section", "Settings")
    pairs = []
    if z:
        pairs.append(("Charge range", "%d to %d" % tuple(z)))
    if p.get("mass_range"):
        pairs.append(("Mass range", "%g to %g Da" % tuple(p["mass_range"])))
    pairs.append(("Mass step", "%g Da" % step))
    di = res.get("data_info") or {}
    if di.get("mz"):
        pairs.append(("m/z range", "%.2f to %.2f" % tuple(di["mz"])))
    elif p.get("mz_range"):
        pairs.append(("m/z range", "%.2f to %.2f" % tuple(p["mz_range"])))
    if di:
        import ms_deconv as _D
        pairs.append(("Minimum intensity", _D._fmt_int(di["threshold"]) if di.get("threshold", 0) > 0 else "None"))
    if p.get("resolution"):
        pairs.append(("Resolving power", "%.0f" % p["resolution"]))
    if p.get("peak_width"):
        pairs.append(("Peak width", "%g m/z" % p["peak_width"]))
    rep.add("info", pairs)
    return _plain_report(rep, fields)


def _si_draw(ax, kind, s, lab):
    """One panel of the Supporting Information page drawn from the data (when the window shows no
    tile of it, e.g. the isotope pattern)."""
    if kind == "chromp":
        ax.plot(s["t"], s["y"], color=BLUE, lw=0.6)
        for p in s["peaks"]:
            i0, i1 = int(p["i0"]), int(p["i1"])
            ts, ys = s["t"][i0:i1 + 1], s["y"][i0:i1 + 1]
            if len(ts) > 1:
                base = p["b0"] + (p["b1"] - p["b0"]) * (ts - ts[0]) / max(ts[-1] - ts[0], 1e-9)
                ax.fill_between(ts, base, np.maximum(ys, base), color=ACCENT, alpha=0.22, lw=0)
        ax.annotate("%.2f" % s["main"]["rt"], (s["main"]["rt"], s["main"].get("apex_y", s["y"].max())),
                    xytext=(0, 2), textcoords="offset points", ha="center", va="bottom", fontsize=6.3,
                    fontweight="bold")
        style(ax, "Time (min)", s["unit"])
        ax.set_xlim(float(s["t"][0]), float(s["t"][-1]))
        headroom(ax, s["y"], 0.15)
    elif kind == "zoom":
        spec, fi = s
        if fi and fi.get("pattern") is not None:
            pat = np.asarray(fi["pattern"], float)
            z = abs(int(fi.get("z") or 1)) or 1
            _zoom_axes(ax, spec, pat[0, 0] - 0.8 / z, pat[-1, 0] + 0.8 / z, pat, fi.get("scale"))
        else:
            pk = local_max(spec[:, 0], spec[:, 1], n=1)
            c = pk[0][0] if pk else float(spec[np.argmax(spec[:, 1]), 0])
            _zoom_axes(ax, spec, c - 1.6, c + 2.8)
    else:
        _panel(ax, kind, s, lab, nb=5)


def _si_panel_png(path, kind, s, lab, w, h):
    f, ax = new_fig(w, h)
    _si_draw(ax, kind, s, lab)
    return save_fig(f, path, left=1.3 / w, right=1 - 0.25 / w, top=1 - 0.45 / h, bottom=0.95 / h)


def _compose(path, pngs, ncol, pw, ph, gap=0.6):
    """Panels (images of pw x ph cm at 300 dpi) on one sheet, ncol per row, lettered a, b, c, d."""
    matplotlib, Figure, FigureCanvasAgg = _mpl()
    from PIL import Image as PILImage
    nrow = (len(pngs) + ncol - 1) // ncol
    W, H = W_FULL, nrow * ph + (nrow - 1) * 0.3
    with matplotlib.rc_context(_RC):
        f = Figure(figsize=(W * CM, H * CM), dpi=300, facecolor="white")
        FigureCanvasAgg(f)
        px = 300.0 * CM  # pixels per cm
        for i, p in enumerate(pngs):
            r, c = divmod(i, ncol)
            im = PILImage.open(p).convert("RGB")
            tw, th = int(round(pw * px)), int(round(ph * px))
            if im.size != (tw, th):
                im = im.resize((tw, th), PILImage.LANCZOS)
            x0 = c * (pw + gap) * px
            y0 = (H - (r + 1) * ph - r * 0.3) * px
            f.figimage(np.asarray(im), xo=int(round(x0)), yo=int(round(y0)), origin="upper", zorder=1)
            f.text((c * (pw + gap) + 0.05) / W, 1 - (r * (ph + 0.3) + 0.05) / H, "abcd"[i], fontsize=9,
                   fontweight="bold", va="top", ha="left", zorder=2)
        f.savefig(path, dpi=300, facecolor="white")
    return path


def build_si(frame, fields, tmp):
    si = sample_info(frame)
    name = fields.get("compound") or si.get("sample_name") or frame.base_name()
    rep = Report("si", _meta(frame, "si", name, fields))
    rep.add("title", "Supporting Information page", esc(name))
    panels, cap, text = [], [], ""
    if getattr(frame, "pda", None) is not None and (frame.pda_data is not None or frame.ms_data is not None):
        st = _lc_state(frame)
        main = st["main"]
        ms_tab = getattr(frame, "ms", None)
        spec_cards = {c["e"]: c.get("spec_card") for c in (getattr(ms_tab, "cols", None) or [])}

        def tile(card, **kw):  # the tile as shown (None: the report draws the panel); text as in the report
            return (lambda path, w, h: _tile_png(card, path, w, h, min_scale=TEXT_SCALE, **kw)) \
                if card is not None else None
        if want(fields, "chrom"):
            card = None
            if not st["auto"]:  # peaks integrated only for the report: drawn by the report
                card = getattr(frame.pda, "chrom", None) if st["src"] == "pda" else \
                    (getattr(ms_tab.views[0], "card", None) if getattr(ms_tab, "views", None) else None)
            panels.append(("chromp", st, "%s" % (st["label"] or ""), tile(card)))
            cap.append("%s chromatogram%s" % ("PDA" if st["src"] == "pda" else "MS",
                                              " at %.0f nm (bandwidth %s nm)" % (st["wl"], st["bw"]) if st["wl"]
                                              else ""))
        if st.get("uv") is not None and want(fields, "uv"):
            panels.append(("uv", st["uv"], "UV", tile(getattr(frame.pda, "spec_card", None))
                           if st.get("uv_src") == "window" else None))
            cap.append("UV/Vis spectrum (%s)" % esc(st.get("uv_desc", "")))
        for m in st["ms"][:2]:
            if want(fields, "ms%d" % m["e"]):
                panels.append(("ms", m["spec"], m["label"], tile(spec_cards.get(m["e"])) if m["src"] == "window"
                               else None))
                cap.append("%s mass spectrum (%s)" % (m["label"], esc(m["desc"])))
        pur = 100 * main["area"] / (st["total"] or 1.0)
        mz = "; ".join("<i>m/z</i> (%s) %s" % (m["label"].replace("ESI", "").strip("()"),
                                                 ", ".join("%.1f" % a for a, _ in m["tops"][:2])) for m in st["ms"])
        text = "<b>LC-MS (ESI):</b> t<sub>R</sub> = %.2f min, purity %.1f %% (%s); %s." % (
            main["rt"], pur, ("%.0f nm" % st["wl"]) if st["wl"] else "TIC", mz)
        tail = " (t<sub>R</sub> of the main peak %.2f min)" % main["rt"]
    else:
        st = _hrms_state(frame)
        spec, fi, dec = st["spec"], st["fi"], st["dec"]
        tab = st["tab"]

        def tile(card, **kw):  # the tile as shown (None: the report draws the panel); text as in the report
            return (lambda path, w, h: _tile_png(card, path, w, h, min_scale=TEXT_SCALE, **kw)) \
                if card is not None else None
        if st["traces"] and want(fields, "chrom"):
            tr = st["traces"][0]
            view = next((v for v in (getattr(tab, "views", None) or [])
                         if v.key == st["e"] and getattr(v, "kind", "main") == "main"), None)
            panels.append(("chrom", np.column_stack([tr["t"], tr["y"]]), tr.get("label") or "TIC",
                           tile(getattr(view, "card", None))))
            cap.append("%s chromatogram" % (tr.get("label") or "TIC"))
        if want(fields, "spec"):
            panels.append(("ms", spec, "ESI(%s)" % ("+" if st["data"].adduct_sign(st["e"]) > 0 else "−"),
                           tile(st["col"].get("spec_card"))))
            cap.append("mass spectrum (%s)" % esc(st["col"]["desc"] or ""))
        if want(fields, "iso"):
            panels.append(("zoom", (spec, fi), fi["name"] if fi else "isotope pattern"))
            cap.append("isotope pattern" + (", measured (blue) and calculated (orange)" if fi else ""))
        if dec and want(fields, "dec"):
            dp = getattr(tab, "dec", None)
            ent = getattr(dp, "active", None)
            mt = None
            if ent is not None and ent.res is dec and hasattr(dp, "report_images"):
                def mt(path, w, h):
                    return dp.report_images(ent, tmp, "d_mass", w, h, TEXT_SCALE, with_z=False,
                                            min_scale=TEXT_SCALE)["mass"]
            panels.append(("mass", np.asarray(dec["mass"], float), "deconvoluted", mt))
            cap.append("zero charge mass spectrum")
        if fi and fi.get("found") is not None and fi.get("ok") is False:
            text = "<b>HRMS (ESI):</b> <i>m/z</i> %s calcd for %s %.4f, not found (nearest peak %.4f%s, " \
                   "outside tolerance)." % (
                       esc(fi.get("ion", "")), esc(fi.get("ion_formula") or fi.get("formula", "")), fi["mz"],
                       fi["found"], ", Δ %+.1f ppm" % fi["err"] if fi.get("err") is not None else "")
        elif fi and fi.get("found") is not None:
            text = "<b>HRMS (ESI):</b> <i>m/z</i> %s calcd for %s %.4f, found %.4f (Δ %+.1f ppm)." % (
                esc(fi.get("ion", "")), esc(fi.get("ion_formula") or fi.get("formula", "")), fi["mz"], fi["found"],
                fi["err"])
        else:
            tops = local_max(spec[:, 0], spec[:, 1], n=2, sep=0.5)
            text = "<b>HRMS (ESI):</b> <i>m/z</i> found %s." % (
                ", ".join("%.4f" % a for a, _ in tops))
        tail = ""
    panels = panels[:4]
    n = len(panels)
    if not n:  # no figure chosen: the text block only
        rep.add("box", "Characterisation text", text)
        return rep
    ncol = 2 if n > 1 else 1
    pw = W_FULL if ncol == 1 else (W_FULL - 0.6) / 2.0  # panel width (cm), 0.6 cm between the columns
    ph = 5.5
    pngs = []
    for i, pan in enumerate(panels):
        kind, s_, lab = pan[:3]
        tile = pan[3] if len(pan) > 3 else None
        path = os.path.join(tmp, "d_panel%d.png" % i)
        img = None
        if tile is not None:  # the tile as shown in the window
            try:
                img = tile(path, pw, ph)
            except Exception as ex:
                print("report: no image of the tile for panel %d (%s); the report draws it" % (i, ex))
                img = None
        pngs.append(img or _si_panel_png(path, kind, s_, lab, pw, ph))
    figp = _compose(os.path.join(tmp, "d_sheet.png"), pngs, ncol, pw, ph, gap=0.6)
    rep.add("figure", figp, W_FULL)
    parts = ["(%s) %s" % ("abcd"[i], c) for i, c in enumerate(cap[:4])]
    rep.add("sicap", "<b>Figure S%s.</b> %s of %s. %s%s." % (
        esc(fields.get("fig_no") or "1"), "LC-MS analysis" if tail else "HRMS analysis", esc(name),
        "; ".join(parts), tail))
    rep.add("box", "Characterisation text", text)
    return rep


# ==========================================================================
# comparison report (Compare view of LCMS Analysis: several runs on top of each other)
# ==========================================================================
def _compare_tab(frame):
    """The Compare view of the window (None: not LCMS Analysis, or the view not made yet)."""
    try:
        tab = vars(frame).get("_compare_tab")
    except TypeError:
        tab = None
    if tab is None:
        fn = getattr(frame, "compare_tab", None)
        if callable(fn):
            try:
                tab = fn(create=False)
            except Exception:
                tab = None
    return tab


def _compare_drawn(frame):
    """(Compare view, the traces it draws)."""
    tab = _compare_tab(frame)
    return tab, (list(getattr(tab, "drawn", None) or []) if tab is not None else [])


def compare_shown(frame):
    """Is the Compare view on screen (instead of the views of a file)?"""
    tab = _compare_tab(frame)
    try:
        return bool(tab is not None and tab.IsShown())
    except Exception:
        return False


def _cmp_settings(tab):
    import lcms_compare_core as K
    s = K.default_settings()
    s.update(getattr(tab, "s", None) or {})
    return s


def _cmp_path(doc):
    return (getattr(doc, "attrs", None) or {}).get("path") or ""


def _cmp_info(tab, doc):
    try:
        return tab.info(doc) or {}
    except Exception:
        return {}


def _cmp_label(tab, en):
    try:
        return tab.label_of(en)
    except Exception:
        return en.get("label") or os.path.splitext(os.path.basename(_cmp_path(en["doc"]).rstrip("\\/")))[0]


def _cmp_shift(o):
    """Shift of a drawn trace on the time axis (alignment plus the shift set for its file), min."""
    return float(o.get("align", 0.0) or 0.0) + float(o["entry"].get("shift", 0.0) or 0.0)


def _cmp_signal(s):
    """(short, long) description of the signal compared."""
    k = s.get("signal", "pda")
    pol = "+" if s.get("polarity", "+") == "+" else "−"
    bw = float(s.get("bw", 4.0) or 0.0)
    if k == "pda":
        if s.get("own_wl"):
            return "PDA", "PDA chromatograms, each run at the wavelength of its PDA view (bandwidth %g nm)" % bw
        w = s.get("wl")
        ws = ("%g nm" % float(w)) if w is not None else "PDA"
        return ws, "PDA chromatograms at %s (bandwidth %g nm)" % (ws, bw)
    if k == "max":
        return "PDA max plot", "PDA max plots (highest absorbance at each time)"
    if k == "tic":
        return "TIC ESI(%s)" % pol, "ESI(%s) total ion chromatograms" % pol
    if k == "bpc":
        return "BPC ESI(%s)" % pol, "ESI(%s) base peak chromatograms" % pol
    if s.get("mz") is None:
        return "m/z", "mass chromatograms"
    win = "± %g%s" % (float(s.get("mz_win", 0.5)), " ppm" if s.get("mz_ppm") else "")
    return ("m/z %g ESI(%s)" % (float(s["mz"]), pol),
            "ESI(%s) mass chromatograms of m/z %g %s" % (pol, float(s["mz"]), win))


def _cmp_steps(tab, s, drawn):
    """The processing, as phrases that follow 'the traces were ...'."""
    out = []
    n = int(s.get("smooth") or 0)
    if n >= 3:
        out.append("smoothed over %d points" % n)
    b = s.get("baseline", "none")
    if b == "offset":
        out.append("corrected for a baseline offset (lowest point at zero)")
    elif b == "drift":
        out.append("corrected for baseline drift (straight line between the ends)")
    elif b == "rolling":
        out.append("corrected for the baseline (rolling minimum over %g min)" % float(s.get("rolling_min", 1.0)))
    t0, t1 = s.get("t0"), s.get("t1")
    if t0 is not None and t1 is not None:
        out.append("cut to %.2f to %.2f min" % (t0, t1))
    elif t0 is not None:
        out.append("cut to the time after %.2f min" % t0)
    elif t1 is not None:
        out.append("cut to the time before %.2f min" % t1)
    bdoc = getattr(tab, "blank_doc", None)
    if bdoc is not None:
        ben = next((en for en in (getattr(tab, "entries", None) or []) if en["doc"] is bdoc), None)
        out.append("corrected by subtracting the blank run %s" % (_cmp_label(tab, ben) if ben else ""))
    if s.get("align") is not None:
        out.append("aligned on the peak at %.2f min (highest point within ± %g min)" % (
            float(s["align"]), float(s.get("align_win", 0.3))))
    own = [o for o in drawn if abs(float(o["entry"].get("shift", 0.0) or 0.0)) > 1e-9]
    if own:
        out.append("shifted in time as set for %d run%s" % (len(own), "s" if len(own) > 1 else ""))
    sc = s.get("scale", "abs")
    if sc == "max":
        out.append("scaled each to its tallest peak (100 %)")
    elif sc == "ref" and s.get("ref_t") is not None:
        out.append("scaled each to its peak at %.2f min (100 %%)" % float(s["ref_t"]))
    return out


def _cmp_layout(s, n):
    import lcms_compare_core as K
    if n < 2:
        return "one trace"
    lay = s.get("layout", "stacked")
    if lay not in ("stacked", "offset"):
        return "overlaid"
    sp, sk = int(round(K.spacing_of(s))), int(round(K.skew_of(s)))
    if lay == "offset":
        t = "offset (each run %d %% of the tallest trace above" % sp
        if sk:
            t += " and %d %% of the time span to the right of" % sk
        t += " the previous one; the first run in front)" if s.get("reverse") else " the next one; the last run in front)"
        return t
    t = "stacked, spacing %d %% of the tallest trace" % sp
    if sk:
        t += ", each moved right by %d %% of the time span" % sk
    return t + (", first run at the bottom" if s.get("reverse") else ", first run at the top")


def _cmp_join(parts):
    parts = [p for p in parts if p]
    if len(parts) < 2:
        return "".join(parts)
    return ", ".join(parts[:-1]) + " and " + parts[-1]


def _cmp_unit(s, drawn):
    # % only for traces really put in % (the view's "scaled": not without a reference time)
    units = set("%" if o.get("scaled", s.get("scale", "abs") != "abs") else o.get("units", "") for o in drawn)
    return units.pop() if len(units) == 1 else ""


def _cmp_region_table(tab):
    """(header, rows) of the areas of the selected region as the view's table gives them, or None without a
    region."""
    if getattr(tab, "area_range", None) is None:
        return None
    try:
        tab._prepare_area_action()  # the region as typed in its fields (also without Enter)
    except Exception as ex:
        print("report: region areas not refreshed (%s)" % ex)
    tb = getattr(tab, "_area_table", None)
    if getattr(tab, "area_range", None) is None or not tb:
        return None
    return [str(h) for h in tb[0]], [[esc(str(v)) for v in row] for row in tb[1]]


def _cmp_peak_table(tab, drawn, s):
    """Peak table of the comparison (main peak of each trace, area % at the guide lines)."""
    tb = getattr(tab, "_peak_table", None) or getattr(tab, "_table", None)
    if tb and len(tb) == 2 and len(tb[1]) == len(drawn):
        hdr, rows = [str(h) for h in tb[0]], [list(r) for r in tb[1]]
    else:  # the view's table is not there (or for other traces): made the same way
        import lcms_compare_core as K
        hdr, rows, _pk = K.peak_rows([{"label": o["label"], "proc": o["proc"], "units": o.get("units", ""),
                                       "shift": _cmp_shift(o)} for o in drawn],
                                     s.get("guides") or [], s.get("integ_thr", 1.0))
    # "Main peak (min)" is the run's own retention time (as the labels of the plot); the height header
    # carries the unit, e.g. "Height (mAU)"
    keep = [i for i, h in enumerate(hdr) if h in ("Trace", "Main peak (min)", "Area %", "Area") or
            h.startswith("Height") or h.startswith("Area (") or
            h.lower().endswith("min area %")]
    ix = {h: i for i, h in enumerate(hdr)}
    shifts = [_cmp_shift(o) for o in drawn]
    run_col = any(abs(v) > 1e-9 for v in shifts) and "Main peak (min)" in ix
    out_h, out_r, rts, areas = [], [], [], []
    for i in keep:
        out_h.append(hdr[i])
        if hdr[i] == "Main peak (min)" and run_col:
            out_h.append("On the plot (min)")
    for r, sh in zip(rows, shifts):
        row = []
        rt = None
        for i in keep:
            v = str(r[i]) if i < len(r) else ""
            row.append(esc(v))
            if hdr[i] == "Main peak (min)":
                try:
                    rt = float(v)
                except ValueError:
                    rt = None
                if run_col:  # where the shifted trace shows it
                    row.append("%.3f" % (rt + sh) if rt is not None else "")
        try:
            areas.append(float(r[ix["Area %"]]) if "Area %" in ix else None)
        except (ValueError, IndexError):
            areas.append(None)
        rts.append(rt)
        out_r.append(row)
    return out_h, out_r, rts, areas


def _cmp_t_end(doc):
    a = getattr(doc, "attrs", None) or {}
    pda, ms = a.get("pda_data"), a.get("ms_data")
    try:
        return float(pda.times[-1]) if pda is not None else float(np.max(ms.rt))
    except Exception:
        return None


def _cmp_notes(tab, drawn, s):
    """Plain notes: runs left out (and why), open files not compared, different methods, run lengths,
    injection volumes."""
    import lcms_compare_core as K
    notes = []
    ents = list(getattr(tab, "entries", None) or [])
    bdoc = getattr(tab, "blank_doc", None)
    shown = {id(o["entry"]["doc"]) for o in drawn}
    for en in ents:
        if not en.get("include") or en["doc"] is bdoc or id(en["doc"]) in shown:
            continue
        why = "no data in the time window"
        try:
            r = tab.raw(en)
            if isinstance(r, Exception):
                why = str(r)
        except Exception:
            pass
        notes.append("%s is left out: %s." % (_cmp_label(tab, en), why))
    if bdoc is not None:
        ben = next((en for en in ents if en["doc"] is bdoc), None)
        try:
            r = tab.raw(ben) if ben is not None else None
            if isinstance(r, Exception):
                notes.append("The blank %s has %s: nothing was subtracted." % (_cmp_label(tab, ben), r))
        except Exception:
            pass
    off = [_cmp_label(tab, en) for en in ents if not en.get("include") and en["doc"] is not bdoc]
    if off:
        notes.append("Not ticked: %s." % ", ".join(off))
    infos = [{"label": o["label"], "t_end": _cmp_t_end(o["entry"]["doc"]),
              "method": _cmp_info(tab, o["entry"]["doc"]).get("method_file")} for o in drawn]
    notes += K.consistency_notes(infos)
    def _num(v):
        try:
            return (0, float(v), v)
        except ValueError:
            return (1, 0.0, v)
    vols = sorted(set(str(_cmp_info(tab, o["entry"]["doc"]).get("inj_vol") or "").strip() for o in drawn) - {""},
                  key=_num)
    if len(vols) > 1:
        notes.append("Different injection volumes: %s µL." % ", ".join(vols))
    if s.get("scale") == "ref" and s.get("ref_t") is None:
        notes.append("No reference time: the traces are on the same scale.")
    return notes


def _compare_fig(path, drawn, s, w, h):
    """The comparison drawn from the traces of the view (when the image of its tile is not available)."""
    from matplotlib.transforms import blended_transform_factory
    import lcms_compare_core as K
    f, ax = new_fig(w, h)
    st = K.offset_layout(s)
    tr = blended_transform_factory(ax.transAxes, ax.transData)
    lo = hi = None
    top = max([float(o["proc"].get("top", 0.0) or 0.0) for o in drawn] + [0.0])
    for n, o in enumerate(drawn):
        pr = o["proc"]
        x = np.asarray(pr["t"], float) + float(o.get("dx", 0.0))
        y = np.asarray(pr["y"], float) + float(o.get("dy", 0.0))
        col = o.get("colour") or PAL[n % len(PAL)]
        ax.plot(x, y, color=col, lw=float(s.get("lw", 0.8) or 0.8), zorder=3 + (len(drawn) - n) * 0.01)
        if s.get("fill"):
            ax.fill_between(x, float(o.get("dy", 0.0)), y, color=col, alpha=0.16, lw=0)
        if len(x):
            lo = float(x[0]) if lo is None else min(lo, float(x[0]))
            hi = float(x[-1]) if hi is None else max(hi, float(x[-1]))
        if (st or len(drawn) == 1) and s.get("labels", "right") in ("right", "left", "outside") and len(y):
            right = s.get("labels", "right") == "right"
            k = max(1, len(y) // 12)
            edge = y[-k:] if right else y[:k]
            ax.text(0.995 if right else 0.005, float(np.nanmax(edge)) + 0.04 * top, o["label"], transform=tr,
                    ha="right" if right else "left", va="bottom", fontsize=6.5, color=col)
    if lo is not None and hi > lo:
        ax.set_xlim(lo, hi)
    for g_ in s.get("guides") or []:
        ax.axvline(g_, color="#7B8594", lw=0.7, ls=(0, (4, 3)), zorder=2)
    if s.get("align") is not None:
        ax.axvline(float(s["align"]), color=ACCENT, lw=0.6, ls=(0, (1, 2)), zorder=2)
    unit = _cmp_unit(s, drawn)
    ylab = {"%": "Relative intensity (%)", "mAU": "Absorbance (mAU)", "counts": "Intensity (counts)"}.get(unit,
                                                                                                         "Intensity")
    stack_bar = K.yaxis_mode(s, len(drawn)) == "bar"
    style(ax, "Time (min)", "" if stack_bar else ylab)
    ys = [np.asarray(o["proc"]["y"], float) + float(o.get("dy", 0.0)) for o in drawn if len(o["proc"]["y"])]
    if ys:
        ax.set_ylim(min(0.0, min(float(np.nanmin(v)) for v in ys)), max(float(np.nanmax(v)) for v in ys) * 1.08
                    or 1.0)
    if K.yaxis_mode(s, len(drawn)) == "none":
        from matplotlib.ticker import NullLocator
        ax.set_yticks([])
        ax.yaxis.set_minor_locator(NullLocator())
    if stack_bar:
        from matplotlib.ticker import NullLocator
        ax.set_yticks([])
        ax.yaxis.set_minor_locator(NullLocator())
        bar = K.nice(0.4 * top)
        if bar > 0:
            ax.plot([-0.012, -0.012], [0.0, bar], transform=tr, color=INK, lw=1.4, solid_capstyle="butt",
                    clip_on=False)
            ax.text(-0.02, bar / 2.0, ("%g %s" % (bar, unit)).strip(), transform=tr, ha="right", va="center",
                    rotation=90, fontsize=6.5, clip_on=False)
    elif len(drawn) > 1 and s.get("labels", "right") != "none":
        from matplotlib.lines import Line2D
        ax.legend([Line2D([], [], color=o.get("colour") or PAL[n % len(PAL)], lw=1.4) for n, o in enumerate(drawn)],
                  [o["label"] for o in drawn], loc="upper right", fontsize=6.5, frameon=False, handlelength=1.4)
    return save_fig(f, path, left=1.0 / w, right=0.975, top=1 - 0.3 / h, bottom=0.95 / h)


def _cmp_mz(tab):
    """The peak picked with the m/z tool of the Compare view and its spectrum tiles (CompareTab.mz_report),
    or None."""
    try:
        fn = getattr(tab, "mz_report", None)
        return fn() if fn is not None else None
    except Exception as ex:
        print("report: spectra of the picked peak not available (%s)" % ex)
        return None


def _cmp_mz_caption(mz):
    """Caption of the spectra of the picked peak: run, time range, delay, averaging, background, and the
    suggested neutral mass worded as a suggestion."""
    pk = mz["pick"]
    pols = [p for p, _ in mz["cards"]]
    names = " and ".join("ESI+" if p == "+" else "ESI−" for p in pols)
    sps = [pk["spectra"][p] for p in pols]
    one = len(pols) == 1
    kind = pk.get("kind")
    what = "time range" if kind == "range" else "scan" if kind == "scan" else "peak"
    where = "picked in the Compare view" if not mz.get("plot", True) else "marked in the comparison plot" if \
        mz.get("marked", True) else \
        "picked in the Compare view (not marked in the comparison plot: its trace or the %s is not shown there)" % what
    cap = "Mass spectrum (%s) of the %s of %s %s: " % (names, what, mz["label"], where) if one else \
        "Mass spectra (%s) of the %s of %s %s: " % (names, what, mz["label"], where)
    rngs = [sp["range"] for sp in sps]
    own = [sp["bg"] for sp in sps if sp.get("bg") and sp.get("bg_manual")]
    if kind in ("range", "scan"):  # dragged across the trace, or Alt + click (one scan)
        r = rngs[0]
        if all(q[1] - q[0] < 1e-9 for q in rngs):
            cap += "the scan nearest to %.2f min (%s)" % (pk["ms_apex"], ", ".join(
                "%s at %.2f min" % ("ESI+" if p == "+" else "ESI−", q[0]) for p, q in zip(pols, rngs)))
        else:
            cap += "averaged %.2f to %.2f min (range)" % (r[0], r[1])
        if own:
            cap += ", background %.2f to %.2f min" % tuple(own[0])
        notes = []
        if not all(q[1] - q[0] < 1e-9 for q in rngs):
            notes.append(", ".join("%s %d scans" % ("ESI+" if p == "+" else "ESI−", sp["n"]) for p, sp in zip(pols, sps)))
        if pk.get("pda"):
            notes.append("MS detector delay %.2f min" % pk.get("delay", 0.0))
        if notes:
            cap += " (%s)" % "; ".join(notes)
        cap += "; bin width %g m/z" % pk.get("binw", 0.05)
        cap += ". Numbers: m/z of the strongest ions."
        if mz.get("suggestion"):
            txt = mz["suggestion"].replace("Suggested M", "M", 1)
            cap += " Suggested neutral mass (to be confirmed): %s." % txt
        return esc(cap)
    if pk.get("kind") == "none":
        cap += "no peak at %.2f min, the scans around that time averaged (%s)" % (pk["ms_apex"], ", ".join(
            "%.2f to %.2f min" % r for r in rngs[:1]))
    elif all(r[1] - r[0] < 1e-9 for r in rngs):
        cap += "the scan nearest to the peak top (%s)" % ", ".join(
            "%s at %.2f min" % ("ESI+" if p == "+" else "ESI−", r[0]) for p, r in zip(pols, rngs))
    else:
        r = rngs[0]
        cap += "the %d scans from the start to the end of the peak averaged (%.2f to %.2f min)" % (
            max(sp["n"] for sp in sps), r[0], r[1]) if len(set(sp["n"] for sp in sps)) == 1 else \
            "the scans from the start to the end of the peak averaged (%.2f to %.2f min; %s)" % (
                r[0], r[1], ", ".join("%s %d scans" % ("ESI+" if p == "+" else "ESI−", sp["n"])
                                      for p, sp in zip(pols, sps)))
    if pk.get("kind") != "none":
        cap += ", peak top at %.2f min" % pk["ms_apex"]
    if pk.get("pda"):
        cap += " (MS detector delay %.2f min)" % pk.get("delay", 0.0)
    cap += "; bin width %g m/z" % pk.get("binw", 0.05)
    bgs = [sp["bg"] for sp in sps if sp.get("bg")]
    if own:
        cap += "; background %.2f to %.2f min subtracted" % tuple(own[0])
    elif bgs:
        cap += "; the spectrum just before the peak (%.2f to %.2f min) subtracted" % bgs[0]
    cap += ". Numbers: m/z of the strongest ions."
    if mz.get("suggestion"):
        txt = mz["suggestion"].replace("Suggested M", "M", 1)
        cap += " Suggested neutral mass (to be confirmed): %s." % txt
    return esc(cap)


def build_compare(frame, fields, tmp):
    region = _cmp_region_table(_compare_tab(frame))  # (first: the view may be drawn again with the region typed)
    tab, drawn = _compare_drawn(frame)
    if not drawn:
        raise RuntimeError("no trace in the Compare view")
    s = _cmp_settings(tab)
    n = len(drawn)
    short, long_ = _cmp_signal(s)
    if n == 1:  # one run: one chromatogram
        long_ = long_.replace("chromatograms", "chromatogram").replace("max plots", "max plot")
    steps = _cmp_steps(tab, s, drawn)
    layout = _cmp_layout(s, n)
    docs = [o["entry"]["doc"] for o in drawn]
    infos = [_cmp_info(tab, d) for d in docs]
    names = [os.path.basename(_cmp_path(d).rstrip("\\/")) for d in docs]
    compound = fields.get("compound") or ""
    sample = compound or ("%d runs" % n if n > 1 else (infos[0].get("sample_name") or names[0]))
    import unidec_theme as T
    files_txt = names[0] + ("" if n == 1 else " and %d more" % (n - 1) if n > 2 else " and " + names[1])
    meta = {"kind": KIND_NAMES["compare"], "sample": sample, "version": getattr(T, "APP_VERSION", ""),
            "date": datetime.datetime.now().strftime("%d.%m.%Y %H:%M"),
            "footer": "Compare view of LCMS Analysis  ·  " + files_txt}
    rep = Report("compare", meta)
    rep.add("title", "Comparison report", "  ·  ".join(
        [esc(compound)] * bool(compound) + ["%d run%s" % (n, "s" if n > 1 else ""), esc(long_[:1].upper() + long_[1:])]))
    # information
    insts = []
    for p in [_cmp_path(d) for d in docs]:
        try:
            i_ = _instrument(type("F", (), {"path": p})())
        except Exception:
            i_ = ""
        if i_ and i_ not in insts:
            insts.append(i_)
    pairs = []
    if compound:
        pairs.append(("Compound", esc(compound)))
    pairs.append(("Runs", "%d: %s" % (n, esc(", ".join(names))) if n <= 4 else
                  "%d: %s and %d more%s" % (n, esc(", ".join(names[:3])), n - 3,
                                            " (table of the runs)" if want(fields, "files") else "")))
    if insts:
        pairs.append(("Instrument", esc(", ".join(insts))))
    acq = sorted(a for a in (i.get("acquired") for i in infos) if isinstance(a, datetime.datetime))
    if acq:
        pairs.append(("Acquired", fmt_date(acq[0]) if len(acq) == 1 or acq[0] == acq[-1] else
                      "%s to %s" % (fmt_date(acq[0]), fmt_date(acq[-1]))))
    pairs.append(("Signal", esc(long_[:1].upper() + long_[1:])))
    pairs.append(("Processing", esc((lambda t: t[:1].upper() + t[1:])("; ".join(steps))) if steps else
                  "none (the traces as recorded)"))
    pairs.append(("Display", esc(layout[:1].upper() + layout[1:])))
    ops = []
    for i in infos:
        if i.get("operator") and i["operator"] not in ops and i["operator"] != "System Administrator":
            ops.append(i["operator"])
    op = fields.get("operator") or ", ".join(ops)
    if op:
        pairs.append(("Operator", esc(op)))
    if fields.get("notes"):
        pairs.append(("Notes", esc(fields["notes"]).replace("\n", "<br/>")))
    rep.add("info", pairs)
    p_hdr, p_rows, rts, areas = _cmp_peak_table(tab, drawn, s)
    selected_region = region is not None
    items = [("%d" % n, "runs compared" if n > 1 else "run"), (short, "signal")]
    rr = [v for v in rts if v is not None]
    if rr:
        lo, hi = min(rr), max(rr)
        items.append(("%.2f min" % lo if hi - lo < 0.005 else "%.2f to %.2f min" % (lo, hi),
                       "retention time of the main peak" + (" in the runs" if n > 1 else "")))
    aa = [v for v in areas if v is not None]
    if aa:
        lo, hi = min(aa), max(aa)
        items.append(("%.1f %%" % lo if hi - lo < 0.05 else "%.1f to %.1f %%" % (lo, hi), "area % of the main peak"))
    if selected_region:
        ra = [r["relative"] for r in tab.area_results if r.get("relative") is not None]
        if ra:
            lo, hi = min(ra), max(ra)
            items.append(("%.1f %%" % lo if hi - lo < 0.05 else "%.1f to %.1f %%" % (lo, hi),
                          "relative area of the selected region"))
        a, b = tab.area_range
        reference = next((r["label"] for r in tab.area_results if r["reference"]), "none")
        rep.add("info", [("Selected region", "%.4f to %.4f min" % (a, b)),
                         ("100 % reference", esc(reference))])
    rep.add("result", items)
    notes = _cmp_notes(tab, drawn, s)
    # the comparison as shown
    tile = False
    if want(fields, "plot"):
        rep.add("section", "Comparison, %s" % short if short not in ("PDA",) else "Comparison")
        import lcms_compare_core as K
        st = s.get("layout", "stacked") == "stacked" and n > 1
        h = min(12.0, max(9.0, 7.5 + 0.8 * n)) if st else 8.5
        rng = (9.0, 12.0) if st else (7.0, 10.0)
        p = _tile_png(getattr(tab, "card", None), os.path.join(tmp, "e_compare.png"), W_FULL, h, h_range=rng,
                      min_scale=TEXT_SCALE)
        tile = p is not None
        if p is None:
            p = _compare_fig(os.path.join(tmp, "e_compare.png"), drawn, s, W_FULL, h)
        rep.add("figure", p, W_FULL)
        cap = "%s%s of %s" % ("As shown in the Compare view of LCMS Analysis: " if tile else "",
                              long_ if tile else long_[:1].upper() + long_[1:],
                              "%d runs, %s" % (n, layout) if n > 1 else "one run")
        if steps:
            cap += "; " + "; ".join(steps)
        cap += "."
        guides = sorted(s.get("guides") or [])
        if guides:
            cap += " Dashed lines: guide lines at %s min." % _cmp_join(
                ["%.2f" % v for v in guides])
        if s.get("align") is not None:
            cap += " Dotted line: the alignment time."
        mz_ = _cmp_mz(tab)
        if tile and mz_ is not None and want(fields, "mzspec") and mz_.get("marked", True):
            if mz_.get("kind") == "range":
                cap += " Shading on the trace of %s: the time range averaged for the mass spectra that follow." % esc(
                    mz_["label"])
            elif mz_.get("kind") == "scan":
                cap += " Ring on the trace of %s: the time of the scan whose mass spectra follow." % esc(mz_["label"])
            else:
                cap += " Ring and shading on the trace of %s: the peak whose mass spectra follow." % esc(mz_["label"])
        bgs = list(getattr(tab, "_mz_bgs", None) or []) if getattr(tab, "mz_on", False) else []
        if tile and bgs:  # background ranges of the m/z tool (Shift + drag), drawn on their traces
            cap += " Grey hatched band%s on the trace%s of %s: the background range%s subtracted from %s mass spectra." % (
                ("s", "s", _cmp_join([esc(v) for v in bgs]), "s", "their") if len(bgs) > 1 else
                ("", "", esc(bgs[0]), "", "its"))
        if K.yaxis_mode(s, n) == "bar":
            # the bar of the figure shown: the tile's (lcms_compare: half the tallest trace) or the one drawn
            # here (_compare_fig: 0.4 of it); the caption said 200 mAU under a tile with a 500 mAU bar
            bar = K.nice((0.5 if tile else 0.4) * max([float(o["proc"].get("top", 0.0) or 0.0) for o in drawn] +
                                                       [0.0]))
            if bar > 0:
                cap += " Bar at the left: %s." % esc(("%g %s" % (
                    bar, _cmp_unit(s, drawn))).strip())
        if s.get("rt_labels", "none") != "none":
            cap += " Numbers at the peaks: retention times%s." % (
                " in each run (before the shift of its trace)" if any(abs(_cmp_shift(o)) > 1e-9 for o in drawn)
                else "")
        rep.add("caption", cap)
    # the spectra of the peak picked with the m/z tool (as its two tiles show them)
    mz = _cmp_mz(tab)
    if mz is not None and want(fields, "mzspec"):
        figs = []
        for pol, card in mz["cards"]:
            w = (W_FULL - 0.4) / 2.0 if len(mz["cards"]) > 1 else 11.0
            p = _tile_png(card, os.path.join(tmp, "e_cmp_ms%s.png" % ("pos" if pol == "+" else "neg")), w, 6.2,
                          min_scale=TEXT_SCALE)
            if p is not None:
                figs.append((p, w))
        if figs:
            rep.add("section", "Mass spectra of the picked peak")
            if len(figs) == 2:
                rep.add("pair", [("figure", figs[0][0], figs[0][1])], [("figure", figs[1][0], figs[1][1])],
                        figs[0][1] + 0.1)
            else:
                rep.add("figure", figs[0][0], figs[0][1])
            rep.add("caption", _cmp_mz_caption(dict(mz, plot=tile)))  # (the mark is in the plot of the tile only)
    # the runs
    if want(fields, "files"):
        rep.add("section", "Runs compared")
        hdr = ["Label", "File", "Sample", "Sample ID", "Acquired", "Method", "Vial", "Injected (µL)",
               "Shift (min)"]
        rows = []
        for o, d, si, fn in zip(drawn, docs, infos, names):
            sh = _cmp_shift(o)
            acq_ = si.get("acquired")
            rows.append([esc(o["label"]), esc(fn), esc(si.get("sample_name") or ""), esc(si.get("sample_id") or ""),
                         esc(fmt_date(acq_)) if acq_ else "", esc(si.get("method_file") or ""),
                         esc(si.get("vial") or ""), esc(si.get("inj_vol") or ""),
                         "%+.3f" % sh if abs(sh) > 1e-9 else "0"])
        bdoc = getattr(tab, "blank_doc", None)
        if bdoc is not None:
            ben = next((en for en in (getattr(tab, "entries", None) or []) if en["doc"] is bdoc), None)
            si = _cmp_info(tab, bdoc)
            acq_ = si.get("acquired")
            rows.append(["Blank: " + esc(_cmp_label(tab, ben) if ben else ""),
                         esc(os.path.basename(_cmp_path(bdoc).rstrip("\\/"))), esc(si.get("sample_name") or ""),
                         esc(si.get("sample_id") or ""), esc(fmt_date(acq_)) if acq_ else "",
                         esc(si.get("method_file") or ""), esc(si.get("vial") or ""), esc(si.get("inj_vol") or ""),
                         "%+.3f" % float(ben.get("shift", 0.0) or 0.0) if ben and abs(float(ben.get("shift", 0.0) or
                                                                                         0.0)) > 1e-9 else "0"])
        rep.add("table", hdr, rows, [2.7, 3.1, 2.0, 1.4, 1.9, 3.1, 0.9, 1.3, 1.4], None)
        small = []
        if n > 1:
            small.append("In the order of the plot%s." % (
                "" if s.get("layout", "stacked") not in ("stacked", "offset") else
                (" (the first run at the bottom)" if s.get("layout") == "stacked" else " (the first run in front)")
                if s.get("reverse") else (" (the first run at the top)" if s.get("layout") == "stacked" else
                                          " (the first run at the back)")))
        if bdoc is not None:
            small.append("Blank: subtracted from every trace.")
        if small:
            rep.add("small", " ".join(small))
    if notes:
        rep.add("text", "<b>Note%s:</b> %s" % ("s" if len(notes) > 1 else "", esc(" ".join(notes))))
    # the peaks
    if want(fields, "peaks") and p_rows:
        rep.add("section", "Peaks of each trace")
        widths = []
        for h_ in p_hdr:
            widths.append(4.2 if h_ == "Trace" else 2.6 if h_.startswith("Height") or h_.startswith("Area (")
                          else 1.7 if h_ == "Area %" else 2.4)
        rep.add("table", [esc(h_) for h_ in p_hdr], p_rows, widths, None)
        txt = "Main peak: the tallest (threshold %g %% of it). Height and area unscaled." % float(
            s.get("integ_thr", 1.0) or 1.0)
        if "On the plot (min)" in p_hdr:
            txt += " Main peak: time in the run; on the plot: on the shifted time axis."
        rep.add("small", txt)
    if want(fields, "peaks") and region is not None:
        r_hdr, r_rows = region
        rep.add("section", "Selected region areas")
        rep.add("table", [esc(h_) for h_ in r_hdr], r_rows,
                [3.6, 2.4, 2.6, 2.8, 2.4, 1.4, 2.4][:len(r_hdr)] + [2.4] * max(0, len(r_hdr) - 7), None)
        txt = "Region peak: time in the run (before its shift)."
        if getattr(tab, "_area_x_error", None):
            txt += " X not used: %s." % esc(tab._area_x_error)
        rep.add("small", txt)
    # methods
    if want(fields, "methods"):
        meths = sorted(set(i.get("method_file") for i in infos if i.get("method_file")))
        ver = getattr(T, "APP_VERSION", "")
        if all(_cmp_path(d).lower().endswith(".lcd") for d in docs):
            txt = "The LC-MS runs were recorded on %s (LabSolutions%s). " % (
                esc(" and ".join("a " + i_.replace(" (LabSolutions)", "") for i_ in insts)) or "a Shimadzu LC-MS",
                ", method %s" % esc(meths[0]) if len(meths) == 1 else "")
        else:  # other vendors (lcms_sources)
            txt = "The runs were recorded on %s%s. " % (
                esc(" and ".join(insts)) or "an LC system", " (method %s)" % esc(meths[0]) if len(meths) == 1 else "")
        verbs = steps + ["stacked" if (s.get("layout", "stacked") == "stacked" and n > 1) else
                         "offset in a waterfall plot" if (s.get("layout") == "offset" and n > 1) else
                         "overlaid" if n > 1 else "drawn"]
        txt += "The %s of %s %s with MSpektra%s (LCMS Analysis, Compare view)." % (
            esc(long_), "the %d runs were" % n if n > 1 else "the run was", esc(_cmp_join(verbs)),
            (" " + esc(ver)) if ver else "")
        if want(fields, "peaks") and p_rows:
            txt += (" Peaks were integrated automatically in each trace (threshold %g %% of the tallest peak)." %
                    float(s.get("integ_thr", 1.0) or 1.0))
        rep.add("section", "Methods")
        rep.add("box", "Methods text", txt)
    return rep


BUILDERS = {"hrms": build_hrms, "lcms": build_lcms, "deconv": build_deconv, "si": build_si,
            "compare": build_compare}


# ==========================================================================
# PDF
# ==========================================================================
def render_pdf(rep, path):
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.units import mm, cm
    from reportlab.lib import colors
    from reportlab.lib.styles import ParagraphStyle
    from reportlab.platypus import (BaseDocTemplate, PageTemplate, Frame, Paragraph, Spacer, Table, TableStyle,
                                    Image, KeepTogether, Flowable)
    from reportlab.pdfgen import canvas as rl_canvas
    from PIL import Image as PILImage
    F = _register_fonts()
    sans, serif = F["sans"], F["serif"]
    bold = sans + "-Bold" if sans.startswith("R") else "Helvetica-Bold"
    C_NAVY, C_ACC, C_INK = colors.HexColor(NAVY), colors.HexColor(ACCENT), colors.HexColor(INK)
    C_MUTED, C_FAINT, C_LINE = colors.HexColor("#5B6475"), colors.HexColor("#8A94A6"), colors.HexColor("#D5DAE3")
    C_ZEBRA, C_HEAD, C_HL = colors.HexColor("#F4F6FA"), colors.HexColor("#E9EDF5"), colors.HexColor("#E3ECFB")
    PW, PH = A4
    LM = RM = 16 * mm
    CW = PW - LM - RM
    S = {
        "title": ParagraphStyle("t", fontName=bold, fontSize=16, leading=19, textColor=C_NAVY),
        "sub": ParagraphStyle("s", fontName=sans, fontSize=9, leading=12, textColor=C_MUTED),
        "body": ParagraphStyle("b", fontName=sans, fontSize=8.5, leading=11.5, textColor=C_INK),
        "small": ParagraphStyle("sm", fontName=sans, fontSize=7.3, leading=9.5, textColor=C_MUTED, spaceBefore=2),
        "cap": ParagraphStyle("c", fontName=sans, fontSize=7.6, leading=10, textColor=C_MUTED, spaceBefore=1,
                              spaceAfter=4),
        "cell": ParagraphStyle("ce", fontName=sans, fontSize=7.8, leading=9.8, textColor=C_INK),
        "key": ParagraphStyle("k", fontName=sans, fontSize=7.8, leading=9.8, textColor=C_MUTED),
        "si": ParagraphStyle("si", fontName=serif, fontSize=10, leading=14, textColor=C_INK),
    }

    def P(t, st="body"):
        return Paragraph(t.replace("<sup>", "<super>").replace("</sup>", "</super>"), S[st])

    class Section(Flowable):
        def __init__(self, n, text):
            Flowable.__init__(self)
            self.n, self.text, self.height = n, text, 17

        def wrap(self, aw, ah):
            self.width = aw
            return aw, self.height

        def draw(self):
            c = self.canv
            c.setFillColor(C_ACC)
            c.roundRect(0, 3, 15, 11, 2.5, stroke=0, fill=1)
            c.setFillColor(colors.white)
            c.setFont(bold, 7.5)
            c.drawCentredString(7.5, 6, str(self.n))
            c.setFillColor(C_NAVY)
            c.setFont(bold, 10.5)
            c.drawString(21, 5, self.text)
            c.setStrokeColor(C_LINE)
            c.setLineWidth(0.5)
            c.line(21 + c.stringWidth(self.text, bold, 10.5) + 6, 8.5, self.width, 8.5)

    class ResultBox(Flowable):
        def __init__(self, items):
            Flowable.__init__(self)
            self.items, self.width = items, CW
            w = self.width / max(1, len(items))
            # a label too long for one line at 6 pt is written on two (it ran into the next box)
            self.lines = [self._split(lab, w - 18) for _, lab in items]
            self.height = 44 if any(len(ln) > 1 for ln in self.lines) else 36

        @staticmethod
        def _split(lab, room):
            from reportlab.pdfbase.pdfmetrics import stringWidth
            if stringWidth(lab, sans, 6.0) <= room:
                return [lab]
            words, a = lab.split(" "), []
            while words and stringWidth(" ".join(a + [words[0]]), sans, 6.0) <= room:
                a.append(words.pop(0))
            if not a:
                a.append(words.pop(0))
            return [" ".join(a), " ".join(words)] if words else [" ".join(a)]

        def wrap(self, aw, ah):
            return self.width, self.height

        def draw(self):
            c = self.canv
            c.setFillColor(colors.HexColor("#F2F6FD"))
            c.setStrokeColor(colors.HexColor("#C9D7F2"))
            c.roundRect(0, 0, self.width, self.height, 5, stroke=1, fill=1)
            w = self.width / max(1, len(self.items))
            two = self.height > 36
            for i, (val, lab) in enumerate(self.items):
                x, room = i * w + 10, w - 18
                fs = 12.0
                while fs > 7 and c.stringWidth(val, bold, fs) > room:
                    fs -= 0.5
                c.setFillColor(C_INK)
                c.setFont(bold, fs)
                c.drawString(x, 25 if two else 17, val)
                lines = self.lines[i]
                c.setFillColor(C_MUTED)
                if len(lines) == 1:
                    fl = 7.0
                    while fl > 5.2 and c.stringWidth(lab, sans, fl) > room:
                        fl -= 0.3
                    c.setFont(sans, fl)
                    c.drawString(x, 13 if two else 7, lab)
                else:
                    for k, ln in enumerate(lines):
                        fl = 6.0
                        while fl > 4.6 and c.stringWidth(ln, sans, fl) > room:
                            fl -= 0.3
                        c.setFont(sans, fl)
                        c.drawString(x, 14 - 7.5 * k, ln)
                if i:
                    c.setStrokeColor(colors.HexColor("#C9D7F2"))
                    c.line(i * w, 6, i * w, self.height - 6)

    def image(path, wcm):
        w, h = PILImage.open(path).size
        width = min(wcm * cm, CW)
        return Image(path, width=width, height=width * h / w)

    def info(pairs):
        # a long value (Notes of a few thousand characters) made its row taller than the page, which a table
        # cannot split (LayoutError): such values follow the table as paragraphs of their own
        vw = CW / 2 - 24 * mm - 6
        long_ = [(k, v) for k, v in pairs if P(v, "cell").wrap(vw, PH)[1] > 12 * S["cell"].leading]
        pairs = [kv for kv in pairs if kv not in long_]
        out = [info_table(pairs)] if pairs else []
        for k, v in long_:
            out += [Spacer(1, 4), P(k, "key"), P(v, "cell")]
        return out

    def info_table(pairs):
        half = (len(pairs) + 1) // 2
        left, right = pairs[:half], pairs[half:]
        rows = []
        for i in range(half):
            k1, v1 = left[i]
            k2, v2 = right[i] if i < len(right) else ("", "")
            rows.append([P(k1, "key"), P(v1, "cell"), P(k2, "key"), P(v2, "cell")])
        t = Table(rows, colWidths=[24 * mm, CW / 2 - 24 * mm, 24 * mm, CW / 2 - 24 * mm])
        t.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"), ("TOPPADDING", (0, 0), (-1, -1), 2.2),
                               ("BOTTOMPADDING", (0, 0), (-1, -1), 2.2), ("LEFTPADDING", (0, 0), (-1, -1), 3),
                               ("LINEBELOW", (0, 0), (-1, -1), 0.3, C_LINE)]))
        return t

    def table(header, rows, widths, hl):
        body = [[P("<b>%s</b>" % h, "cell") for h in header]] + [[P(str(c), "cell") for c in r] for r in rows]
        tw = sum(widths) * cm
        k = min(1.0, CW / tw) if tw else 1.0
        t = Table(body, colWidths=[w * cm * k for w in widths], repeatRows=1, hAlign="LEFT")
        st = [("VALIGN", (0, 0), (-1, -1), "MIDDLE"), ("BACKGROUND", (0, 0), (-1, 0), C_HEAD),
              ("LINEBELOW", (0, 0), (-1, 0), 0.6, C_NAVY), ("LINEBELOW", (0, -1), (-1, -1), 0.6, C_NAVY),
              ("TOPPADDING", (0, 0), (-1, -1), 2), ("BOTTOMPADDING", (0, 0), (-1, -1), 2),
              ("LEFTPADDING", (0, 0), (-1, -1), 3), ("RIGHTPADDING", (0, 0), (-1, -1), 3)]
        for i in range(2, len(body), 2):
            st.append(("BACKGROUND", (0, i), (-1, i), C_ZEBRA))
        if hl is not None:
            st.append(("BACKGROUND", (0, hl + 1), (-1, hl + 1), C_HL))
        t.setStyle(TableStyle(st))
        return t

    def flow(b, n_sec):
        kind = b[0]
        if kind == "title":
            return [P(esc(b[1]), "title"), Spacer(1, 2), P(b[2], "sub"), Spacer(1, 6)]
        if kind == "info":
            return info(b[1]) + [Spacer(1, 7)]
        if kind == "result":
            return [ResultBox(b[1])]
        if kind == "section":
            n_sec[0] += 1
            return [Section(n_sec[0], b[1])]
        if kind == "figure":
            return [image(b[1], b[2])]
        if kind == "caption":
            return [P(b[1], "cap")]
        if kind == "table":
            return [table(b[1], b[2], b[3], b[4]), Spacer(1, 3)]
        if kind == "text":
            return [P(b[1], "body"), Spacer(1, 2)]
        if kind == "small":
            return [P(b[1], "small")]
        if kind == "sicap":
            return [Spacer(1, 4), P(b[1], "si"), Spacer(1, 10)]
        if kind == "box":
            inner = [P("<b>%s</b>" % esc(b[1]), "small"), Spacer(1, 3), P(b[2], "si")]
            t = Table([[inner]], colWidths=[CW])
            t.setStyle(TableStyle([("BOX", (0, 0), (-1, -1), 0.6, colors.HexColor("#C9D7F2")),
                                   ("BACKGROUND", (0, 0), (-1, -1), colors.HexColor("#F7F9FD")),
                                   ("LEFTPADDING", (0, 0), (-1, -1), 8), ("RIGHTPADDING", (0, 0), (-1, -1), 8),
                                   ("TOPPADDING", (0, 0), (-1, -1), 6), ("BOTTOMPADDING", (0, 0), (-1, -1), 8)]))
            return [t]
        if kind == "pair":
            left = sum((flow(x, n_sec) for x in b[1]), [])
            right = sum((flow(x, n_sec) for x in b[2]), [])
            t = Table([[left, right]], colWidths=[b[3] * cm, CW - b[3] * cm])
            t.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"), ("LEFTPADDING", (0, 0), (-1, -1), 0),
                                   ("LEFTPADDING", (1, 0), (1, 0), 6)]))
            return [t]
        return []

    story, n_sec = [], [0]
    blocks = rep.blocks
    i = 0
    while i < len(blocks):
        b = blocks[i]
        if b[0] == "section":  # a heading stays with what follows it
            grp = flow(b, n_sec)
            j = i + 1
            while j < len(blocks) and len(grp) < 4 and blocks[j][0] not in ("section",):
                grp += flow(blocks[j], n_sec)
                j += 1
                if blocks[j - 1][0] in ("figure", "table", "pair"):
                    break
            story.append(KeepTogether(grp))
            i = j
        else:
            story += flow(b, n_sec)
            i += 1
    meta = rep.meta

    class NC(rl_canvas.Canvas):
        def __init__(self, *a, **k):
            rl_canvas.Canvas.__init__(self, *a, **k)
            self._pg = []

        def showPage(self):
            self._pg.append(dict(self.__dict__))
            self._startPage()

        def _fit(self, text, font, size, room):
            """text shortened with an ellipsis to room points (a long sample name ran over the header)."""
            if self.stringWidth(text, font, size) <= room:
                return text
            lo, hi = 0, len(text)
            while lo < hi:
                mid = (lo + hi + 1) // 2
                if self.stringWidth(text[:mid].rstrip() + "…", font, size) <= room:
                    lo = mid
                else:
                    hi = mid - 1
            return text[:lo].rstrip() + "…"

        def save(self):
            n = len(self._pg)
            for s_ in self._pg:
                self.__dict__.update(s_)
                self.saveState()
                self.setFillColor(C_NAVY)
                self.setFont(bold, 8.5)
                self.drawString(LM, PH - 12 * mm, "MSpektra")
                self.setFillColor(C_MUTED)
                self.setFont(sans, 8.5)
                x0 = LM + self.stringWidth("MSpektra ", bold, 8.5)
                self.drawString(x0, PH - 12 * mm, "·  " + meta["kind"])
                tail = "   ·   " + meta["date"]
                room = (PW - RM) - (x0 + self.stringWidth("·  " + meta["kind"], sans, 8.5) + 12) - \
                    self.stringWidth(tail, sans, 8.5)
                self.drawRightString(PW - RM, PH - 12 * mm, self._fit(meta["sample"], sans, 8.5, room) + tail)
                self.setStrokeColor(C_ACC)
                self.setLineWidth(0.8)
                self.line(LM, PH - 14 * mm, PW - RM, PH - 14 * mm)
                self.setStrokeColor(C_LINE)
                self.setLineWidth(0.4)
                self.line(LM, 12.5 * mm, PW - RM, 12.5 * mm)
                self.setFont(sans, 6.8)
                self.setFillColor(C_FAINT)
                pg = "Page %d of %d" % (self._pageNumber, n)
                room = CW - self.stringWidth(pg, sans, 6.8) - 12
                self.drawString(LM, 8.5 * mm, self._fit("MSpektra %s  ·  %s" % (meta["version"], meta["footer"]),
                                                        sans, 6.8, room))
                self.drawRightString(PW - RM, 8.5 * mm, pg)
                self.restoreState()
                rl_canvas.Canvas.showPage(self)
            rl_canvas.Canvas.save(self)

    doc = BaseDocTemplate(path, pagesize=A4, leftMargin=LM, rightMargin=RM, topMargin=19 * mm, bottomMargin=17 * mm,
                          title="%s: %s" % (meta["kind"], meta["sample"]), author="MSpektra %s" % meta["version"])
    doc.addPageTemplates([PageTemplate(id="p", frames=[Frame(LM, 17 * mm, CW, PH - 36 * mm, id="f", leftPadding=0,
                                                             rightPadding=0, topPadding=0, bottomPadding=0)])])
    doc.build(story, canvasmaker=NC)
    return path


# ==========================================================================
# Word
# ==========================================================================
def _runs(par, markup, size=None, color=None, font=None):
    """Add simple markup (<b>, <i>, <sub>, <sup>/<super>, <br/>) as runs."""
    from docx.shared import Pt, RGBColor
    state = {"b": False, "i": False, "sub": False, "sup": False}
    for tok in re.split(r"(<[^>]+>)", markup or ""):
        if not tok:
            continue
        if tok.startswith("<"):
            t = tok.strip("<>/ ").lower()
            closing = tok.startswith("</")
            if t in ("br", "br/"):
                par.add_run().add_break()
                continue
            t = {"super": "sup", "strong": "b", "em": "i"}.get(t, t)
            if t in state:
                state[t] = not closing
            continue
        r = par.add_run(html.unescape(tok))
        r.bold = state["b"] or None
        r.italic = state["i"] or None
        if state["sub"]:
            r.font.subscript = True
        if state["sup"]:
            r.font.superscript = True
        if size:
            r.font.size = Pt(size)
        if color:
            r.font.color.rgb = RGBColor.from_string(color.lstrip("#"))
        if font:
            r.font.name = font
    return par


def _shade(cell, hexcol):
    from docx.oxml import OxmlElement
    from docx.oxml.ns import qn
    tcPr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:val"), "clear")
    shd.set(qn("w:color"), "auto")
    shd.set(qn("w:fill"), hexcol.lstrip("#"))
    tcPr.append(shd)


def render_docx(rep, path):
    import docx
    from docx.shared import Pt, Cm
    from docx.enum.table import WD_TABLE_ALIGNMENT
    d = docx.Document()
    sec = d.sections[0]
    sec.page_width, sec.page_height = Cm(21.0), Cm(29.7)
    sec.left_margin = sec.right_margin = Cm(1.6)
    sec.top_margin, sec.bottom_margin = Cm(1.8), Cm(1.6)
    st = d.styles["Normal"]
    st.font.name = "Arial"
    st.font.size = Pt(9)
    meta = rep.meta
    hp = sec.header.paragraphs[0]
    _runs(hp, "<b>MSpektra</b> · %s    %s · %s" % (esc(meta["kind"]), esc(meta["sample"]),
                                                                   esc(meta["date"])), 8, "#5B6475")
    fp = sec.footer.paragraphs[0]
    _runs(fp, "MSpektra %s · %s" % (esc(meta["version"]), esc(meta["footer"])), 7, "#8A94A6")
    n_sec = [0]

    def para(markup, size=9, color="#1B2330", space=3, font=None, align=None):
        p = d.add_paragraph()
        p.paragraph_format.space_after = Pt(space)
        p.paragraph_format.space_before = Pt(0)
        if align:
            p.alignment = align
        _runs(p, markup, size, color, font)
        return p

    def table(header, rows, widths, hl, target=None):
        t = (target or d).add_table(rows=1 + len(rows), cols=len(header))
        t.alignment = WD_TABLE_ALIGNMENT.LEFT
        t.style = "Table Grid" if "Table Grid" in [s.name for s in d.styles] else None
        for j, h in enumerate(header):
            c = t.rows[0].cells[j]
            c.paragraphs[0].text = ""
            _runs(c.paragraphs[0], "<b>%s</b>" % h, 8)
            _shade(c, "#E9EDF5")
        for i, r in enumerate(rows):
            for j, v in enumerate(r):
                c = t.rows[i + 1].cells[j]
                _runs(c.paragraphs[0], str(v), 8)
                if hl is not None and i == hl:
                    _shade(c, "#E3ECFB")
        for j, w in enumerate(widths[:len(header)]):
            t.columns[j].width = Cm(w)  # the grid too (LibreOffice reads only the grid)
            for row in t.rows:
                row.cells[j].width = Cm(w)
        return t

    def info(pairs):
        half = (len(pairs) + 1) // 2
        t = d.add_table(rows=half, cols=4)
        for i in range(half):
            for k, (key, val) in enumerate([pairs[i], pairs[half + i] if half + i < len(pairs) else ("", "")]):
                c1, c2 = t.rows[i].cells[2 * k], t.rows[i].cells[2 * k + 1]
                _runs(c1.paragraphs[0], esc(key) if "<" not in key else key, 8, "#5B6475")
                _runs(c2.paragraphs[0], val, 8)
        for row in t.rows:
            for j, w in enumerate((2.5, 6.4, 2.5, 6.4)):
                row.cells[j].width = Cm(w)
        para("", 4)

    for b in rep.blocks:
        kind = b[0]
        if kind == "title":
            para("<b>%s</b>" % esc(b[1]), 16, NAVY, 1)
            para(b[2], 9, "#5B6475", 6)
        elif kind == "info":
            info(b[1])
        elif kind == "result":
            t = d.add_table(rows=2, cols=len(b[1]))
            for j, (val, lab) in enumerate(b[1]):
                _runs(t.rows[0].cells[j].paragraphs[0], "<b>%s</b>" % esc(val), 12)
                _runs(t.rows[1].cells[j].paragraphs[0], esc(lab), 7, "#5B6475")
                _shade(t.rows[0].cells[j], "#F2F6FD")
                _shade(t.rows[1].cells[j], "#F2F6FD")
            para("", 4)
        elif kind == "section":
            n_sec[0] += 1
            p = para("<b>%d  %s</b>" % (n_sec[0], esc(b[1])), 11, NAVY, 3)
            p.paragraph_format.space_before = Pt(8)
            p.paragraph_format.keep_with_next = True
        elif kind == "figure":
            p = d.add_paragraph()
            p.paragraph_format.keep_with_next = True
            p.add_run().add_picture(b[1], width=Cm(min(b[2], 17.8)))
        elif kind == "caption":
            para(b[1], 7.5, "#5B6475", 5)
        elif kind == "table":
            table(b[1], b[2], b[3], b[4])
            para("", 3)
        elif kind == "text":
            para(b[1], 9)
        elif kind == "small":
            para(b[1], 7.5, "#5B6475")
        elif kind == "sicap":
            para(b[1], 10, "#1B2330", 10, "Times New Roman")
        elif kind == "box":
            t = d.add_table(rows=1, cols=1)
            c = t.rows[0].cells[0]
            _shade(c, "#F7F9FD")
            _runs(c.paragraphs[0], "<b>%s</b>" % esc(b[1]), 7.5, "#5B6475")
            p2 = c.add_paragraph()
            _runs(p2, b[2], 10, "#1B2330", "Times New Roman")
        elif kind == "pair" and all(x[0] == "figure" for x in b[1] + b[2]):
            p = d.add_paragraph()  # two figures side by side (e.g. a mass spectrum and its charge states)
            p.paragraph_format.keep_with_next = True
            for x in b[1] + b[2]:
                p.add_run().add_picture(x[1], width=Cm(x[2]))
        elif kind == "pair":
            for x in b[1] + b[2]:
                if x[0] == "figure":
                    d.add_paragraph().add_run().add_picture(x[1], width=Cm(x[2]))
                elif x[0] == "table":
                    table(x[1], x[2], x[3], x[4])
                elif x[0] == "small":
                    para(x[1], 7.5, "#5B6475")
    d.save(path)
    return path


# ==========================================================================
# entry point
# ==========================================================================
def create(frame, kind, fields, path, fmt="pdf"):
    """Build the report of the file shown in frame and write it to path."""
    ensure_libs()
    _APEX["on"] = bool(getattr(getattr(frame, "ms", None), "APEX_LABELS", False)
                       and getattr(getattr(frame, "ms_data", None), "has_profile", False))
    tmp = tempfile.mkdtemp(prefix="msreport_")
    try:
        rep = BUILDERS[kind](frame, fields, tmp)
        if fmt == "docx":
            render_docx(rep, path)
        else:
            render_pdf(rep, path)
    finally:
        try:
            import shutil
            shutil.rmtree(tmp, ignore_errors=True)
        except Exception:
            pass
    return path

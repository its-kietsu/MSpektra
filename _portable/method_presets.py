"""Method presets of LCMS Postrun and HRMS Postrun.

A method preset is a named snapshot of every setting of a window kind (the settings of its views, of
the deconvolution, of the Compare view or the calibration), not of the data of a file (ranges, peaks,
mass chromatogram lists, labels, names, colours and shifts of the runs) and not of the view (zoom,
tile heights, the session only Graph properties). Presets are kept in the settings file
(config/portable_settings.json, key "method_presets"), can be exported to a file and imported from
one (name.msmethod.json), and one preset per window kind can be the default, applied when a file is
opened in that window.

The windows have short get and set functions (MSTab.method_get/method_set, PDATab, HRMSTab
exact_get/finder_get, CalibPanel, CompareTab); everything else is here. apply_preset() is the one
entry point that applies a preset. The storage, the validation and the files need no wx (unit tests:
cpp/msengine/tests/check_presets.py); the menu, the dialogs and the window side import wx when used.
"""
import copy
import json
import math
import os
import time

import unidec_theme as T
import lcms_compare_core as K
import hrms_calib as HC
import ms_formula as MF
import ms_shifts as MSH

KEY = "method_presets"  # settings: {kind: {name: {"version", "created", "modified", "values"}}, "default": {kind: name}}
FORMAT = "MS Analysis method"
FORMAT_VERSION = 1
EXT = ".msmethod.json"
KINDS = {"lcms": "LCMS Postrun", "hrms": "HRMS Postrun"}
_FRAMES = {"LCMSFrame": "lcms", "HRMSFrame": "hrms"}
NAME_MAX = 60


class MethodError(Exception):
    """A preset or a method file that cannot be used (the reason in plain words)."""


# ==========================================================================
# what a preset contains: section, key, label, value type, program default
# ==========================================================================
def _bool():
    return ("bool",)


def _int(lo, hi):
    return ("int", lo, hi)


def _float(lo, hi, none=False):
    return ("float", lo, hi, none)


def _choice(*keys):
    return ("choice", tuple(keys))


def _floats(lo, hi, n):
    return ("floats", lo, hi, n)


NUMTEXT = ("numtext",)  # a number typed in a field of the deconvolution settings, or empty

SECTION_LABELS = {"ms": "Mass spectrometry view", "pda": "PDA view", "window": "MS and PDA link",
                  "deconvolution": "Deconvolution", "compare": "Compare view", "calibration": "Calibration",
                  "exact_mass": "Exact mass", "formula_finder": "Formula finder", "mass_shifts": "Mass shifts"}
# sections of each file (applied to the file or to every open file) and of the window (all files)
PER_FILE = {"lcms": ("ms", "pda"), "hrms": ("ms", "exact_mass", "formula_finder")}
WINDOW_WIDE = {"lcms": ("deconvolution", "mass_shifts", "compare", "window"),
               "hrms": ("deconvolution", "mass_shifts", "calibration")}

_INTEGRATION = [("integration_min_height", "Peak integration, minimum height (% of top)", _float(0, 100), 2.0),
                ("integration_min_width", "Peak integration, minimum width (s)", _float(0, 3600), 2.0)]


def _ms_fields(hr):
    f = [("trace", "Trace", _choice("tic", "bpc", "xic"), "tic"),
         ("xic_window", "Default mass chromatogram window", _float(1e-6, 1e5), 0.02 if hr else 0.5),
         ("xic_window_unit", "Unit of the mass chromatogram window", _choice("m/z", "ppm"), "m/z"),
         ("xic_shown", "Mass chromatograms shown", _choice("stacked", "overlaid"), "stacked"),
         ("smoothing", "Smoothing of the chromatograms (points)", _int(0, 51), 0),
         ("scale_100", "Scale each trace to 100 %", _bool(), False),
         ("link_time_axes", "Link the time axes", _bool(), True),
         ("subtract_background", "Subtract the background", _bool(), True)]
    if not hr:  # (HRMS spectra keep the profile m/z axis of the instrument: no bin width)
        f.append(("bin_width", "Bin width of the spectra (m/z)", _float(0.001, 10.0), 0.05))
    f.append(("spectrum_display", "Spectrum display", _choice("profile", "sticks"), "profile"))
    return f + _INTEGRATION


PDA_FIELDS = [("wavelength", "Wavelength of the slider (nm)", _float(100, 1100, True), None),
              ("bandwidth", "Bandwidth (nm)", _float(0, 200), 4.0),
              ("colour_scale", "Colour scale of the map", _choice("compressed", "linear"), "compressed"),
              ("more_wavelengths", "More wavelengths (nm)", _floats(100, 1100, 20), []),
              ("max_plot", "Max plot", _bool(), False),
              ("smoothing", "Smoothing of the chromatograms (points)", _int(0, 51), 0),
              ("scale_100", "Scale each trace to 100 %", _bool(), False),
              ("overlay_ms", "Overlay the MS trace", _bool(), False),
              ("ms_delay", "MS detector delay (min)", _float(-30, 30), 0.0),
              ("subtract_background", "Subtract the background", _bool(), True)] + _INTEGRATION

WINDOW_FIELDS = [("link_ms_pda", "Link MS and PDA times", _bool(), True)]


def _keys(pairs):
    return [k for k, _ in pairs]


_CS = K.default_settings()
COMPARE_SESSION = ("t0", "t1", "align", "ref_t", "guides")  # kept by the view for the session only (in a preset too)
COMPARE_FIELDS = [
    ("signal", "Signal", _choice(*_keys(K.SIGNALS)), _CS["signal"]),
    ("wl", "Wavelength (nm)", _float(100, 1100, True), _CS["wl"]),
    ("bw", "Bandwidth (nm)", _float(0, 200), _CS["bw"]),
    ("own_wl", "Each run at its own wavelength", _bool(), _CS["own_wl"]),
    ("polarity", "Polarity", _choice("+", "-"), _CS["polarity"]),
    ("mz", "m/z of the mass chromatogram", _float(1e-6, 1e6, True), _CS["mz"]),
    ("mz_win", "Window of the mass chromatogram", _float(1e-9, 1e6), _CS["mz_win"]),
    ("mz_ppm", "Window in ppm", _bool(), _CS["mz_ppm"]),
    ("smooth", "Smoothing (points)", _int(0, 51), _CS["smooth"]),
    ("baseline", "Baseline", _choice(*_keys(K.BASELINES)), _CS["baseline"]),
    ("rolling_min", "Rolling minimum window (min)", _float(0.01, 1000), _CS["rolling_min"]),
    ("t0", "Time window from (min)", _float(-1e4, 1e5, True), None),
    ("t1", "Time window to (min)", _float(-1e4, 1e5, True), None),
    ("align", "Align at (min)", _float(-1e4, 1e5, True), None),
    ("align_win", "Alignment search window (min)", _float(0.01, 100), _CS["align_win"]),
    ("layout", "Layout", _choice(*_keys(K.LAYOUTS)), _CS["layout"]),
    ("scale", "Scale", _choice(*_keys(K.SCALES)), _CS["scale"]),
    ("ref_t", "Reference time (min)", _float(-1e4, 1e5, True), None),
    ("spacing", "Spacing, stacked (%)", _int(0, 400), _CS["spacing"]),
    ("skew", "Skew, stacked (%)", _int(0, 30), _CS["skew"]),
    ("off_spacing", "Spacing, offset (%)", _int(0, 400), _CS["off_spacing"]),
    ("off_skew", "Skew, offset (%)", _int(0, 30), _CS["off_skew"]),
    ("reverse", "Reverse order", _bool(), _CS["reverse"]),
    ("colours", "Colours", _choice(*_keys(K.COLOURS)), _CS["colours"]),
    ("lw", "Line width (pt)", _float(0.2, 4.0), _CS["lw"]),
    ("fill", "Fill under the traces", _bool(), _CS["fill"]),
    ("labels", "Names of the runs", _choice(*_keys(K.LABELS)), _CS["labels"]),
    ("label_text", "Name shown", _choice(*_keys(K.LABEL_TEXT)), _CS["label_text"]),
    ("rt_labels", "Retention times", _choice(*_keys(K.RT_LABELS)), _CS["rt_labels"]),
    ("rt_min", "Retention time label threshold (%)", _float(0, 100), _CS["rt_min"]),
    ("guides", "Guide lines (min)", _floats(-1e4, 1e5, 50), []),
    ("yaxis", "Y axis", _choice(*_keys(K.YAXES)), _CS["yaxis"]),
    ("integ_thr", "Peak table threshold (%)", _float(0, 100), _CS["integ_thr"]),
    ("spec_avg", "Spectra of a peak", _choice("peak", "scan"), "peak"),
    ("spec_bg", "Spectra: background before the peak subtracted", _bool(), True),
    ("spec_labels", "Spectra: ions labelled", _int(0, 30), 8),
]

CAL_FIELDS = [("calibrant", "Reference list", _choice(*[c[0] for c in HC.CALIBRANTS]), "Sodium formate, positive"),
              ("model", "Model", _choice(*HC.MODELS), HC.MODELS[5]),
              ("order", "HPC order (0 = automatic)", _int(0, HC.HPC_MAX_ORDER), 0),
              ("tolerance_ppm", "Search window (ppm)", _float(1e-3, 1e4), 30.0),
              ("min_intensity", "Minimum intensity (% of the base peak)", _float(0, 100), 0.5),
              ("custom_list", "Custom reference list (m/z)", _floats(1e-3, 1e6, 500), []),
              ("auto", "Calibrate automatically on opening", _bool(), False)]

EXACT_FIELDS = [("ion", "Ion", _choice(*MF.ADDUCT_LABELS), "[M+H]+"),
                ("tolerance_ppm", "Tolerance (ppm)", _float(1e-6, 1000), 5.0)]
# the elements of the formula finder (hrms.FormulaDialog.ELEMENTS) with their limits
ELEMENTS = [("C", 0, 80), ("H", 0, 160), ("N", 0, 10), ("O", 0, 20), ("S", 0, 3), ("P", 0, 2), ("F", 0, 0),
            ("Cl", 0, 0), ("Br", 0, 0), ("Na", 0, 0), ("Si", 0, 0), ("B", 0, 0)]
FINDER_FIELDS = [("limits", "Element limits", ("limits", tuple(e for e, _, _ in ELEMENTS), 0, 500),
                  {e: [lo, hi] for e, lo, hi in ELEMENTS}),
                 ("rank_by_isotopes", "Rank by the isotope pattern", _bool(), True)]

# the mass shift finder of the deconvolution results (ms_shifts.default_settings; key "mass_shifts"), without its
# display choice (Show every matched pair)
_MSD = MSH.default_settings()
SHIFT_FIELDS = [("shifts", "Shift list (name, average and monoisotopic mass, in use)", ("shiftlist",), _MSD["shifts"]),
                ("tags", "Tags (name, average and monoisotopic mass)", ("taglist",), []),
                ("kmax", "At most k tags per molecule", _int(1, MSH.MAX_K), _MSD["kmax"]),
                ("max_extra", "Other shifts with the tags", _int(0, 2), _MSD["max_extra"]),
                ("tol_envelope", "Tolerance, isotope envelopes", ("tol",), None),
                ("tol_resolved", "Tolerance, isotope resolved results", ("tol",), None),
                ("min_rel", "Masses used (share of the tallest)", _float(0, 1), _MSD["min_rel"])]

# deconvolution: the settings of the deconvolution window (deconv_tab: _default_cfg, SETTINGS), stored as
# that window keeps them (choices as their number in the list, numbers as typed)
_DEC = [("method", "Method", _int(0, 2)), ("preset", "Preset", None), ("pol", "Polarity", _int(0, 1)),
        ("suggest_open", "Suggestions shown", _bool()),
        ("z0", "Charge from", NUMTEXT), ("z1", "Charge to", NUMTEXT), ("m0", "Mass from", NUMTEXT),
        ("m1", "Mass to", NUMTEXT), ("step", "Mass step", NUMTEXT), ("mz0", "m/z range from", NUMTEXT),
        ("mz1", "m/z range to", NUMTEXT), ("min_int", "Minimum intensity", NUMTEXT),
        ("min_int_unit", "Unit of the minimum intensity", _int(0, 1)),
        ("adduct", "Charge carrier (positive ions)", _int(0, 4)),
        ("adduct_neg", "Charge carrier (negative ions)", _int(0, 4)),
        ("adduct_custom", "Carrier mass", NUMTEXT), ("baseline", "Subtract the baseline", _bool()),
        ("bl_width", "Baseline window", NUMTEXT), ("rp", "Resolving power", NUMTEXT),
        ("width", "Peak width", NUMTEXT), ("psfun", "Peak shape", _int(0, 2)), ("iso", "Isotopes", _int(0, 1)),
        ("isotopemode", "Isotope mode", _int(0, 2)), ("disp_rp", "Show at resolving power", NUMTEXT),
        ("numit", "Iterations", NUMTEXT), ("zzsig", "Charge smoothing", NUMTEXT),
        ("psig", "Point smoothing", NUMTEXT), ("beta", "Suppression (beta)", NUMTEXT),
        ("msig", "Mass smoothing", NUMTEXT), ("poolflag", "m/z to mass", _int(0, 2)),
        ("binning", "Data reduction", _int(0, 3)), ("mzbins", "Bin size", NUMTEXT),
        ("smooth", "Data smoothing", NUMTEXT), ("rounds", "Rounds", NUMTEXT),
        ("min_rounds", "Minimum rounds", NUMTEXT), ("iterations", "Iterations per round", NUMTEXT),
        ("phaseres", "IsoDec model", _int(0, 1)), ("iso_matchtol", "Matching tolerance", NUMTEXT),
        ("iso_minpeaks", "Minimum isotope peaks", NUMTEXT), ("iso_css_thresh", "Similarity threshold", NUMTEXT),
        ("iso_knockdown_rounds", "Knock-down rounds", NUMTEXT), ("iso_datathreshold", "Data threshold", NUMTEXT),
        ("iso_window", "Peak window (points)", NUMTEXT), ("iso_thresh", "Peak threshold", NUMTEXT),
        ("win", "Peak window", NUMTEXT), ("thr", "Peak threshold", NUMTEXT), ("labels", "Labels", NUMTEXT)]
N_DECONV_PRESETS = {"lcms": 5, "hrms": 7}  # deconv_tab.PRESETS offered (HRMS adds two isotope resolved ones)


def _dec_fields(kind):
    return [(k, lab, spec if spec is not None else _int(0, N_DECONV_PRESETS[kind] - 1), None) for k, lab, spec in _DEC]


def schema(kind):
    """{section: [(key, label, value type, program default)]} of a window kind."""
    if kind == "lcms":
        return {"ms": _ms_fields(False), "pda": PDA_FIELDS, "window": WINDOW_FIELDS,
                "deconvolution": _dec_fields(kind), "mass_shifts": SHIFT_FIELDS, "compare": COMPARE_FIELDS}
    if kind == "hrms":
        return {"ms": _ms_fields(True), "exact_mass": EXACT_FIELDS, "formula_finder": FINDER_FIELDS,
                "calibration": CAL_FIELDS, "deconvolution": _dec_fields(kind), "mass_shifts": SHIFT_FIELDS}
    raise MethodError("Unknown window kind: %r" % (kind,))


def factory_values(kind):
    """The settings of a window opened for the first time (no deconvolution: its defaults are in
    deconv_tab._default_cfg)."""
    return {sec: {k: copy.deepcopy(d) for k, lab, spec, d in fields}
            for sec, fields in schema(kind).items() if sec != "deconvolution"}


def contents_text(kind):
    """What a preset of this window kind contains, one line per section (help, README)."""
    lines = []
    for sec, fields in schema(kind).items():
        if sec == "deconvolution":
            lines.append("Deconvolution: every setting of the deconvolution window")
            continue
        lines.append("%s: %s" % (SECTION_LABELS[sec], ", ".join(lab for k, lab, spec, d in fields)))
    return "\n".join(lines)


# ==========================================================================
# validation
# ==========================================================================
def _short(v):
    t = json.dumps(v, ensure_ascii=False) if not isinstance(v, str) else '"%s"' % v
    return t if len(t) <= 40 else t[:37] + "..."


def _is_num(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(float(v))


def _g(v):
    return "%g" % v


def check_value(spec, v):
    """(True, the value to use) or (False, why it cannot be used)."""
    kind = spec[0]
    if kind == "bool":
        if isinstance(v, bool):
            return True, v
        if isinstance(v, int) and v in (0, 1):
            return True, bool(v)
        return False, "is not true or false"
    if kind == "int":
        lo, hi = spec[1], spec[2]
        if _is_num(v) and float(v) == int(v) and lo <= int(v) <= hi:
            return True, int(v)
        return False, "is not a whole number from %d to %d" % (lo, hi)
    if kind == "float":
        lo, hi, none = spec[1], spec[2], spec[3]
        if v is None and none:
            return True, None
        if _is_num(v) and lo <= float(v) <= hi:
            return True, float(v)
        return False, "is not a number from %s to %s%s" % (_g(lo), _g(hi), " (or empty)" if none else "")
    if kind == "choice":
        if isinstance(v, str) and v in spec[1]:
            return True, v
        opts = ", ".join(spec[1][:8]) + (", ..." if len(spec[1]) > 8 else "")
        return False, "is not one of: %s" % opts
    if kind == "floats":
        lo, hi, n = spec[1], spec[2], spec[3]
        if isinstance(v, (list, tuple)) and len(v) <= n and all(_is_num(x) and lo <= float(x) <= hi for x in v):
            return True, [float(x) for x in v]
        return False, "is not a list of at most %d numbers from %s to %s" % (n, _g(lo), _g(hi))
    if kind == "numtext":
        if _is_num(v):
            v = "%.10g" % v
        if isinstance(v, str) and len(v.strip()) <= 40:
            t = v.strip()
            if not t:
                return True, ""
            try:
                x = float(t.replace(",", "."))
            except ValueError:
                x = None
            if x is not None and math.isfinite(x):
                return True, t
        return False, "is not a number (or empty)"
    if kind in ("shiftlist", "taglist"):
        key = "shifts" if kind == "shiftlist" else "tags"
        n_max = MSH.MAX_SHIFTS if kind == "shiftlist" else MSH.MAX_TAGS
        what = "shifts (name, average and monoisotopic mass)" if kind == "shiftlist" else \
            "tags (name, average and monoisotopic mass)"
        if isinstance(v, (list, tuple)) and len(v) <= n_max and all(isinstance(e, dict) for e in v):
            out = MSH.clean_settings({key: list(v)})[key]
            if len(out) == len(v):
                return True, out
        return False, "is not a list of at most %d %s" % (n_max, what)
    if kind == "tol":  # as the mass shift finder takes it (ms_shifts.clean_settings: up to 100 Da or 100000 ppm)
        if v is None:
            return True, None
        if isinstance(v, (list, tuple)) and len(v) == 2 and v[1] in ("Da", "ppm") and _is_num(v[0]):
            t = MSH.clean_settings({"tol_envelope": [float(v[0]), v[1]]})["tol_envelope"]
            if t is not None:
                return True, t
        return False, "is not empty (automatic) or a number above 0 in Da (at most 100) or ppm (at most 100000)"
    if kind == "limits":
        elements, lo, hi = spec[1], spec[2], spec[3]
        if not isinstance(v, dict):
            return False, "is not a list of elements with their minimum and maximum"
        out = {}
        for e in elements:
            if e not in v:
                continue
            p = v[e]
            if not (isinstance(p, (list, tuple)) and len(p) == 2 and all(_is_num(x) and float(x) == int(x) and
                                                                         lo <= int(x) <= hi for x in p)
                    and int(p[0]) <= int(p[1])):
                return False, "%s %s is not a minimum and maximum from %d to %d" % (e, _short(p), lo, hi)
            out[e] = [int(p[0]), int(p[1])]
        return True, out
    return False, "has a type this version does not know"


def validate(kind, values):
    """(values that can be used, problems): unknown sections and keys are left out quietly, values
    of the wrong type or out of range are left out with a line in problems (the current setting then
    stays when the preset is applied); missing keys are simply not in the result."""
    sch = schema(kind)
    clean, problems = {}, []
    if not isinstance(values, dict):
        return clean, ["the settings are not a list of named values"]
    for sec, fields in sch.items():
        part = values.get(sec)
        if part is None:
            continue
        if not isinstance(part, dict):
            problems.append("%s: not a list of named values" % SECTION_LABELS[sec])
            continue
        out = {}
        for key, label, spec, dflt in fields:
            if key not in part:
                continue
            ok, v = check_value(spec, part[key])
            if ok:
                out[key] = v
            else:
                shown = "" if spec[0] == "limits" and isinstance(part[key], dict) else _short(part[key]) + " "
                problems.append("%s, %s: %s%s" % (SECTION_LABELS[sec], label, shown, v))
        if out:
            clean[sec] = out
    return clean, problems


def _same(a, b):
    if _is_num(a) and _is_num(b):
        return abs(float(a) - float(b)) <= 1e-6 * max(1.0, abs(float(a)), abs(float(b)))
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        return len(a) == len(b) and all(_same(x, y) for x, y in zip(a, b))
    if isinstance(a, dict) and isinstance(b, dict):
        return set(a) == set(b) and all(_same(a[k], b[k]) for k in a)
    return a == b


TOLERANCE = {("pda", "wavelength"): 1.5}  # the slider takes the nearest wavelength recorded (1.2 nm apart)


def differences(kind, current, preset_values):
    """Labels of the settings of preset_values that current does not have (both validated)."""
    out = []
    labels = {sec: {k: lab for k, lab, spec, d in fields} for sec, fields in schema(kind).items()}
    for sec, part in preset_values.items():
        cur = current.get(sec, {})
        for k, v in part.items():
            tol = TOLERANCE.get((sec, k))
            if tol is not None and _is_num(cur.get(k)) and _is_num(v):
                if abs(cur[k] - v) > tol:
                    out.append("%s, %s" % (SECTION_LABELS[sec], labels[sec].get(k, k)))
                continue
            if k in cur and not _same(cur[k], v):
                out.append("%s, %s" % (SECTION_LABELS[sec], labels[sec].get(k, k)))
    return out


# ==========================================================================
# storage in the settings file
# ==========================================================================
def _now():
    return time.strftime("%Y-%m-%d %H:%M")


def clean_name(name):
    name = " ".join(str(name or "").split())
    if not name:
        raise MethodError("Enter a name for the method")
    if len(name) > NAME_MAX:
        raise MethodError("The name is too long (at most %d characters)" % NAME_MAX)
    return name


def _store():
    """The presets in the settings file, whatever an older version or a hand edit left there."""
    raw = T._load().get(KEY)
    raw = raw if isinstance(raw, dict) else {}
    out = {"default": {}}
    for kind in KINDS:
        d = raw.get(kind)
        out[kind] = {str(n): p for n, p in d.items() if isinstance(p, dict) and isinstance(p.get("values"), dict)} \
            if isinstance(d, dict) else {}
    dflt = raw.get("default")
    if isinstance(dflt, dict):
        out["default"] = {k: v for k, v in dflt.items() if k in KINDS and isinstance(v, str) and v in out[k]}
    return out


def _write(store):
    T._save({KEY: store})


def _check_kind(kind):
    if kind not in KINDS:
        raise MethodError("Unknown window kind: %r" % (kind,))


def presets(kind):
    """{name: preset} of a window kind, sorted by name."""
    _check_kind(kind)
    d = _store()[kind]
    return {n: d[n] for n in sorted(d, key=lambda s: s.lower())}


def names(kind):
    return list(presets(kind))


def get_preset(kind, name):
    p = presets(kind).get(name)
    if p is None:
        raise MethodError('There is no method "%s" in %s' % (name, KINDS[kind]))
    return p


def save_preset(kind, name, values, overwrite=True):
    """Stores the values (validated) under name; returns (name, problems)."""
    _check_kind(kind)
    name = clean_name(name)
    clean, problems = validate(kind, values)
    st = _store()
    old = st[kind].get(name)
    if old is not None and not overwrite:
        raise MethodError('A method "%s" exists already' % name)
    st[kind][name] = {"version": FORMAT_VERSION, "created": (old or {}).get("created") or _now(),
                      "modified": _now(), "values": clean}
    _write(st)
    return name, problems


def rename_preset(kind, old, new):
    _check_kind(kind)
    new = clean_name(new)
    st = _store()
    if old not in st[kind]:
        raise MethodError('There is no method "%s"' % old)
    if new == old:
        return new
    if new in st[kind]:
        raise MethodError('A method "%s" exists already: choose another name' % new)
    st[kind][new] = st[kind].pop(old)
    if st["default"].get(kind) == old:
        st["default"][kind] = new
    _write(st)
    return new


def delete_preset(kind, name):
    _check_kind(kind)
    st = _store()
    if name not in st[kind]:
        raise MethodError('There is no method "%s"' % name)
    del st[kind][name]
    if st["default"].get(kind) == name:
        del st["default"][kind]
    _write(st)


def default_name(kind):
    """The preset applied when a file is opened in this window kind, or None (as before)."""
    return _store()["default"].get(kind)


def set_default(kind, name):
    _check_kind(kind)
    st = _store()
    if name is None:
        st["default"].pop(kind, None)
    else:
        if name not in st[kind]:
            raise MethodError('There is no method "%s"' % name)
        st["default"][kind] = name
    _write(st)


# ==========================================================================
# method files (name.msmethod.json)
# ==========================================================================
def file_text(kind, name, preset=None):
    p = preset or get_preset(kind, name)
    clean, _ = validate(kind, p.get("values"))
    return json.dumps({"format": FORMAT, "format_version": FORMAT_VERSION, "kind": kind, "window": KINDS[kind],
                       "name": name, "created": p.get("created") or _now(), "modified": p.get("modified") or "",
                       "program": "MS Analysis %s" % getattr(T, "APP_VERSION", ""), "values": clean},
                      indent=1, ensure_ascii=False)


def export_file(path, kind, name):
    text = file_text(kind, name)
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write(text)
    os.replace(tmp, path)
    return path


def read_file(path, kind):
    """(name, values, notes) of a method file for this window kind; raises MethodError with the reason
    when it is not one. notes: values left out (wrong type or range) and other remarks."""
    _check_kind(kind)
    base = os.path.basename(path)
    try:
        with open(path, "rb") as fh:
            raw = fh.read(5_000_001)
    except OSError as ex:
        raise MethodError("%s could not be read: %s" % (base, ex.strerror or ex))
    if len(raw) > 5_000_000:
        raise MethodError("%s is not a method file (too large)" % base)
    try:
        d = json.loads(raw.decode("utf-8-sig"))
    except (UnicodeDecodeError, ValueError) as ex:
        raise MethodError("%s is not a method file of MS Analysis (not readable as JSON: %s)" % (base, ex))
    if not isinstance(d, dict) or d.get("format") != FORMAT or not isinstance(d.get("values"), dict):
        raise MethodError("%s is not a method file of MS Analysis (exported with Method > Export to a file)" % base)
    k = d.get("kind")
    if k != kind:
        if k in KINDS:
            raise MethodError("%s is a method of %s; import it there (this window is %s)" % (base, KINDS[k],
                                                                                           KINDS[kind]))
        raise MethodError("%s is a method for a window this version does not know (%r)" % (base, k))
    notes = []
    ver = d.get("format_version")
    if not _is_num(ver):
        notes.append("the file does not say its format version; read as version %d" % FORMAT_VERSION)
    elif ver > FORMAT_VERSION:
        notes.append("the file was written by a newer version (format %s); the settings this version knows "
                     "were read" % _g(ver))
    try:
        name = clean_name(d.get("name") or os.path.splitext(os.path.splitext(base)[0])[0])
    except MethodError:
        name = clean_name(os.path.splitext(base)[0][:NAME_MAX] or "Imported method")
    clean, problems = validate(kind, d["values"])
    if not clean:
        raise MethodError("%s has no setting this version can use%s" % (
            base, (":\n" + "\n".join(problems[:12])) if problems else ""))
    return name, clean, problems + notes


def free_name(kind, name):
    """name, or name (2), (3), ... when a method of that name exists."""
    have = set(names(kind))
    if name not in have:
        return name
    i = 2
    while True:
        n = "%s (%d)" % (name[:NAME_MAX - 6], i)
        if n not in have:
            return n
        i += 1


# ==========================================================================
# the window side (wx)
# ==========================================================================
def kind_of(window):
    """'lcms' or 'hrms' for a window (or one of its files), else None."""
    fr = _frame(window)
    for klass in type(fr).__mro__:
        if klass.__name__ in _FRAMES:
            return _FRAMES[klass.__name__]
    return None


def _frame(window):
    try:
        return window.window()  # PostrunFrame.window(): the frame also for one of its files
    except Exception:
        return window


def _doc_of(window):
    fr = _frame(window)
    return window if window is not fr else getattr(fr, "active", None)


def _attrs(doc):
    return object.__getattribute__(doc, "attrs")


def _open_docs(fr):
    return [d for d in fr.docs if _attrs(d).get("path") or _attrs(d).get("loading")]


def _log(*a):
    try:
        print("[method presets]", *a)
    except Exception:
        pass


def _dec_panel(doc):
    ms = _attrs(doc).get("ms")
    return getattr(ms, "dec", None)


def _dec_key(kind):
    return "deconv_hr_v2" if kind == "hrms" else "deconv_lc_v2"


def _dec_get(doc, kind):
    pn = _dec_panel(doc)
    if pn is not None:
        return dict(pn._saved_cfg() or pn.cfg)
    import deconv_tab
    saved = T._load().get(_dec_key(kind))
    cfg = deconv_tab._default_cfg(kind == "hrms")
    if isinstance(saved, dict):
        cfg.update({k: v for k, v in saved.items() if k in cfg})
    return cfg


def _dec_set(fr, kind, v):
    """The deconvolution settings of the window kind (every file and window of that kind reads them
    when its deconvolution window opens)."""
    import deconv_tab
    key = _dec_key(kind)
    pn = next((_dec_panel(d) for d in fr.docs if _dec_panel(d) is not None), None)
    cfg = (pn._saved_cfg() if pn is not None else None) or deconv_tab._default_cfg(kind == "hrms")
    cfg.update(v)
    T._save({key: cfg})
    for p in list(deconv_tab.DeconvPanel.instances):
        if getattr(p, "key", None) == key:
            p.cfg = dict(cfg)


def _shifts_get():
    """The settings of the mass shift finder (shared by every file and window), without the display choice."""
    s = MSH.clean_settings(T._load().get("mass_shifts"))
    s.pop("show_all", None)
    return s


def _shifts_set(v):
    """The settings of the mass shift finder; the results shown with mass shifts are drawn again."""
    s = MSH.clean_settings(T._load().get("mass_shifts"))
    s.update({k: copy.deepcopy(x) for k, x in v.items() if k != "show_all"})
    s = MSH.clean_settings(s)
    T._save({"mass_shifts": s})
    MSH.SETTINGS_REV[0] += 1  # (undo: a change made in this program)
    try:
        import deconv_tab
        pn = next(iter(list(deconv_tab.DeconvPanel.instances)), None)
        if pn is not None and getattr(pn, "shifts", None) is not None:
            pn.shifts.redraw_all()  # (every panel)
    except Exception as ex:
        _log("mass shifts not drawn again: %s" % ex)


def _compare_tab(fr):
    return fr.__dict__.get("_compare_tab")


def _compare_get(fr):
    tab = _compare_tab(fr)
    if tab is not None:
        return tab.method_get()
    # the view is made when first shown: the settings it would start with
    s = {k: copy.deepcopy(d) for k, lab, spec, d in COMPARE_FIELDS}
    saved = T._load().get("compare")
    if isinstance(saved, dict):
        s.update({k: v for k, v in K.migrate(saved).items() if k in s and k not in COMPARE_SESSION})
    s.update(fr.__dict__.get("_method_compare") or {})
    return s


def _compare_set(fr, v):
    tab = _compare_tab(fr)
    if tab is not None:
        tab.method_set(v)
        return
    # not made yet: saved as the view saves its settings, the rest when it is made (LCMSFrame.compare_tab)
    pend = dict(fr.__dict__.get("_method_compare") or {})
    pend.update(v)
    fr.__dict__["_method_compare"] = pend
    saved = T._load().get("compare")
    saved = dict(saved) if isinstance(saved, dict) else {}
    saved.update({k: x for k, x in v.items() if k not in COMPARE_SESSION})
    T._save({"compare": saved})


def _section_get(fr, doc, kind, sec):
    a = _attrs(doc) if doc is not None else {}
    if sec == "ms":
        return a["ms"].method_get()
    if sec == "pda":
        return a["pda"].method_get()
    if sec == "exact_mass":
        return a["ms"].exact_get()
    if sec == "formula_finder":
        return a["ms"].finder_get()
    if sec == "calibration":
        return a["ms"].cal.method_get()
    if sec == "deconvolution":
        return _dec_get(doc, kind)
    if sec == "compare":
        return _compare_get(fr)
    if sec == "window":
        return {"link_ms_pda": bool(T._load().get("link_ms_pda", True))}
    if sec == "mass_shifts":
        return _shifts_get()
    return {}


def capture(window):
    """(values, problems): the current settings of the window (the file shown, or the file given, and
    the settings of the window), validated."""
    fr = _frame(window)
    kind = kind_of(fr)
    doc = _doc_of(window)
    raw = {}
    for sec in schema(kind):
        try:
            raw[sec] = _section_get(fr, doc, kind, sec)
        except Exception as ex:
            _log("%s not read: %s" % (sec, ex))
    return validate(kind, raw)


def _apply_file(doc, kind, sec, v):
    a = _attrs(doc)
    if sec == "ms":
        a["ms"].method_set(v)
    elif sec == "pda":
        a["pda"].method_set(v)
    elif sec == "exact_mass":
        a["ms"].exact_set(v)
    elif sec == "formula_finder":
        a["ms"].finder_set(v)


def _apply_window(fr, kind, sec, v):
    if sec == "deconvolution":
        _dec_set(fr, kind, v)
    elif sec == "compare":
        _compare_set(fr, v)
    elif sec == "window":
        if "link_ms_pda" in v:
            T._save({"link_ms_pda": bool(v["link_ms_pda"])})
    elif sec == "mass_shifts":
        _shifts_set(v)
    elif sec == "calibration":
        cals = [getattr(_attrs(d).get("ms"), "cal", None) for d in fr.docs]
        for cal in [c for c in cals if c is not None]:
            cal.method_set(v)


def _resolve(kind, preset):
    if isinstance(preset, dict):
        return preset.get("name") or "", preset.get("values") or {}
    return preset, get_preset(kind, preset).get("values") or {}


def _set_chip(doc, name):
    try:
        _attrs(doc)["method_name"] = name
        chips = _attrs(doc).get("chips") or {}
        if "method" in chips:
            chips["method"].SetValue(name or "")
    except Exception as ex:
        _log("chip: %s" % ex)


def apply_preset(window, preset, scope="file"):
    """Applies a method preset (see _apply_preset); one undo step "apply method NAME" unless the file is being
    opened (scope "open": the opened state is the start of the undo history)."""
    if scope == "open":
        return _apply_preset(window, preset, scope)
    import undo
    label = 'apply method "%s"' % (preset.get("name") or "unnamed" if isinstance(preset, dict) else preset)
    if scope == "all":
        label += " to every open file"
    with undo.action(_frame(window), label, full=True):
        return _apply_preset(window, preset, scope)


def _apply_preset(window, preset, scope="file"):
    """Applies a method preset: the one entry point (menu, default on opening, tests).

    window: the window, or one of its files; preset: the name of a stored preset, or a dict with
    "values" (and "name"); scope: "file" (the file given, or the file shown), "all" (every open file),
    "open" (the file given is being opened: called by the window before the file is read, with the
    default preset; the settings of the window are applied too when no other file is open).
    The settings of the files go through the set functions of their views (controls show the values,
    chromatograms and spectra are computed again as after a change in the panel); the settings of the
    window (deconvolution, Compare view, calibration) apply to every file. Returns the message shown in
    the status bar."""
    fr = _frame(window)
    kind = kind_of(fr)
    if kind is None:
        raise MethodError("Method presets are for LCMS Postrun and HRMS Postrun")
    name, raw = _resolve(kind, preset)
    vals, problems = validate(kind, raw)
    doc = _doc_of(window)
    if scope == "all":
        docs = _open_docs(fr) or ([doc] if doc is not None else [])
    elif scope in ("file", "open"):
        docs = [doc] if doc is not None else []
    else:
        raise MethodError("Unknown scope %r (file, all or open)" % (scope,))
    window_parts = scope != "open" or not [d for d in _open_docs(fr) if d is not doc]
    done, failed = [], []
    for sec in PER_FILE[kind]:
        if sec not in vals:
            continue
        for d in docs:
            try:
                _apply_file(d, kind, sec, vals[sec])
            except Exception as ex:
                import traceback
                traceback.print_exc()
                failed.append("%s (%s)" % (SECTION_LABELS[sec], ex))
        done.append(SECTION_LABELS[sec])
    if window_parts:
        for sec in WINDOW_WIDE[kind]:
            if sec not in vals:
                continue
            try:
                _apply_window(fr, kind, sec, vals[sec])
                done.append(SECTION_LABELS[sec])
            except Exception as ex:
                import traceback
                traceback.print_exc()
                failed.append("%s (%s)" % (SECTION_LABELS[sec], ex))
    for d in docs:
        _set_chip(d, name)
    files = [_attrs(d).get("file_name") or "" for d in docs]
    files = [f for f in files if f]
    where = ("every open file (%d)" % len(files)) if scope == "all" and len(files) > 1 else \
        (files[0] if files else "this window")
    msg = 'Method "%s" applied to %s (%s)' % (name or "unnamed", where, ", ".join(done) or "nothing to apply")
    if problems:
        msg += "; %d value%s not valid, left as they were (%s)" % (
            len(problems), "s" if len(problems) > 1 else "", "; ".join(problems[:3]) + ("; ..." if len(problems) > 3
                                                                                      else ""))
    if failed:
        msg += "; not applied: " + "; ".join(failed)
    _log(msg)
    if scope != "open":
        try:
            fr.SetStatusText(msg, 0)
        except Exception:
            pass
    return msg


def file_opening(window, doc):
    """Called by the window when a file is about to be read into doc: the default method of the window
    kind, if one is set (none: nothing changes)."""
    try:
        kind = kind_of(window)
        name = default_name(kind) if kind else None
        if not name:
            return None
        msg = apply_preset(doc, name, "open")
        try:
            _frame(window).SetStatusText("Reading %s \u2026 (default method \"%s\")" % (
                _attrs(doc).get("file_name") or "", name), 0)
        except Exception:
            pass
        return msg
    except Exception as ex:
        import traceback
        traceback.print_exc()
        _log("default method not applied: %s" % ex)
        return None


def current_name(window):
    doc = _doc_of(window)
    return _attrs(doc).get("method_name") if doc is not None else None


def _set_current_everywhere(fr, old, new):
    for d in fr.docs:
        if _attrs(d).get("method_name") == old:
            _set_chip(d, new)


# ----------------------------------------------------------------- actions (menu)
def _status(fr, text):
    try:
        fr.SetStatusText(text, 0)
    except Exception:
        pass
    _log(text)


def _message(fr, text, icon="info"):
    import wx
    wx.MessageBox(text, "Method", {"info": wx.ICON_INFORMATION, "warn": wx.ICON_WARNING,
                                   "error": wx.ICON_ERROR}[icon], fr)


def _ask_name(fr, title, value, note, ok):
    import unilcms as U
    v = U.ask_form(fr, title, [dict(key="name", label="Name", value=value, width=260)], note=note, ok=ok,
                   colour=T.GROUP["blue"])
    return None if v is None else v["name"]


def _confirm(fr, text, yes, no="Cancel"):
    import wx
    dlg = wx.MessageDialog(fr, text, "Method", wx.YES_NO | wx.ICON_QUESTION)
    try:
        dlg.SetYesNoLabels(yes, no)
    except Exception:
        pass
    try:
        return dlg.ShowModal() == wx.ID_YES
    finally:
        dlg.Destroy()


def save_as(window, name=None, replace=None):
    """Save the current settings as a new method (name None: asks). Returns the name or None."""
    fr = _frame(window)
    kind = kind_of(fr)
    values, problems = capture(window)
    while True:
        if name is None:
            name = _ask_name(fr, "Save the current settings as a method", current_name(window) or "",
                             "Saves the settings of %s, not the data of the files." % KINDS[kind], "Save")
            if name is None:
                return None
        try:
            name = clean_name(name)
        except MethodError as ex:
            _message(fr, str(ex), "warn")
            name = None
            continue
        if name in names(kind):
            if replace is None:
                if not _confirm(fr, 'A method "%s" exists already. Replace its settings with the current ones?'
                                % name, "Replace"):
                    name = None
                    continue
            elif not replace:
                return None
        break
    name, _ = save_preset(kind, name, values)
    doc = _doc_of(window)
    if doc is not None:
        _set_chip(doc, name)
    msg = 'Settings saved as the method "%s"' % name
    if problems:
        msg += "; left out (field without a valid value): " + "; ".join(problems[:3])
    _status(fr, msg)
    return name


def update_current(window, name=None):
    fr = _frame(window)
    kind = kind_of(fr)
    name = name or current_name(window)
    if not name or name not in names(kind):
        _status(fr, "No method to update: save the settings as a method first")
        return None
    values, problems = capture(window)
    save_preset(kind, name, values)
    _status(fr, 'Method "%s" updated with the current settings%s' % (
        name, ("; left out: " + "; ".join(problems[:3])) if problems else ""))
    return name


def rename(window, old, new=None):
    fr = _frame(window)
    kind = kind_of(fr)
    while True:
        if new is None:
            new = _ask_name(fr, 'Rename the method "%s"' % old, old, None, "Rename")
            if new is None:
                return None
        try:
            new = rename_preset(kind, old, new)
            break
        except MethodError as ex:
            _message(fr, str(ex), "warn")
            new = None
    _set_current_everywhere(fr, old, new)
    _status(fr, 'Method "%s" renamed to "%s"' % (old, new))
    return new


def delete(window, name, confirm=True):
    fr = _frame(window)
    kind = kind_of(fr)
    if confirm and not _confirm(fr, 'Delete the method "%s"?' % name, "Delete"):
        return False
    delete_preset(kind, name)
    _set_current_everywhere(fr, name, None)
    _status(fr, 'Method "%s" deleted' % name)
    return True


def choose_default(window, name):
    fr = _frame(window)
    kind = kind_of(fr)
    set_default(kind, name)
    _status(fr, ('Default method: "%s"' % name) if name else "No default method")


def export(window, name, path=None):
    import unilcms as U
    fr = _frame(window)
    kind = kind_of(fr)
    if path is None:
        path = U.ask_save_file(fr, 'Export the method "%s"' % name, U.safe_file_name(name) + EXT,
                               "MS Analysis method (*%s)|*%s" % (EXT, EXT))
        if not path:
            return None
    try:
        export_file(path, kind, name)
    except (OSError, MethodError) as ex:
        _message(fr, "The method could not be exported to %s:\n%s" % (path, ex), "error")
        return None
    _status(fr, 'Method "%s" exported to %s' % (name, path))
    return path


def import_(window, path=None, replace=None):
    """Imports a method file (path None: asks). replace: what to do when a method of that name exists
    (None: ask; True: replace it; False: keep both, the imported one gets a new name)."""
    import wx
    fr = _frame(window)
    kind = kind_of(fr)
    if path is None:
        dlg = wx.FileDialog(fr, "Import a method of %s" % KINDS[kind],
                            wildcard="MS Analysis method (*%s;*.json)|*%s;*.json|All files (*.*)|*.*" % (EXT, EXT),
                            style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST)
        try:
            if dlg.ShowModal() != wx.ID_OK:
                return None
            path = dlg.GetPath()
        finally:
            dlg.Destroy()
    try:
        name, values, notes = read_file(path, kind)
    except MethodError as ex:
        _message(fr, str(ex), "warn")
        _status(fr, "Not imported: " + str(ex).split("\n")[0])
        return None
    if name in names(kind):
        if replace is None:
            dlg = wx.MessageDialog(fr, 'A method "%s" exists already.' % name, "Import a method",
                                   wx.YES_NO | wx.CANCEL | wx.ICON_QUESTION)
            try:
                dlg.SetYesNoCancelLabels("Replace it", "Keep both", "Cancel")
            except Exception:
                pass
            try:
                r = dlg.ShowModal()
            finally:
                dlg.Destroy()
            if r == wx.ID_CANCEL:
                return None
            replace = r == wx.ID_YES
        if not replace:
            name = free_name(kind, name)
    save_preset(kind, name, values)
    msg = 'Method "%s" imported from %s' % (name, os.path.basename(path))
    if notes:
        msg += "; %d note%s" % (len(notes), "s" if len(notes) > 1 else "")
        _message(fr, "%s.\n\nLeft out or noted:\n%s" % (msg, "\n".join(notes[:15]) + ("\n..." if len(notes) > 15
                                                                                      else "")), "warn")
    _status(fr, msg)
    return name


def menu_items(window):
    """The entries of the Method menu (unilcms.build_menu format)."""
    fr = _frame(window)
    kind = kind_of(fr)
    have = names(kind)
    cur = current_name(window)
    if cur not in have:
        cur = None
    dflt = default_name(kind)
    head = "Method: %s" % (cur or "none")
    if cur:
        try:
            now, _ = capture(window)
            pv, _ = validate(kind, get_preset(kind, cur)["values"])
            if differences(kind, now, pv):
                head += " (settings changed since)"
        except Exception as ex:
            _log("compare with the preset: %s" % ex)

    def listed(fn, mark=True):
        if not have:
            return [("(no method saved yet)", None)]
        return [(n + ("   (default)" if n == dflt else ""), (lambda n=n: fn(n))) + ((n == cur,) if mark else ())
                for n in have]

    items = [(head, None), (None, None),
             ("Apply to this file", listed(lambda n: apply_preset(fr, n, "file"))),
             ("Apply to every open file", listed(lambda n: apply_preset(fr, n, "all"))),
             (None, None),
             ("Save the current settings as\u2026", lambda: save_as(fr)),
             (('Update "%s" with the current settings' % cur) if cur else "Update this method with the current "
                                                                         "settings",
              (lambda: update_current(fr)) if cur else None),
             (None, None),
             ("Rename", listed(lambda n: rename(fr, n), mark=False)),
             ("Delete", listed(lambda n: delete(fr, n), mark=False)),
             (None, None),
             ("Default when a file is opened",
              [("None", lambda: choose_default(fr, None), dflt is None)] +
              ([(None, None)] + [(n, (lambda n=n: choose_default(fr, n)), n == dflt) for n in have] if have else [])),
             (None, None),
             ("Export to a file", listed(lambda n: export(fr, n), mark=False)),
             ("Import from a file\u2026", lambda: import_(fr)),
             (None, None),
             ("What a method contains\u2026", lambda: _message(fr, "A method of %s contains:\n\n%s\n\nNot in a "
                                                                  "method: the data of the files, the zoom and the "
                                                                  "Graph properties." % (KINDS[kind],
                                                                                         contents_text(kind))))]
    return items


def show_menu(window, anchor=None):
    import wx
    import unilcms as U
    fr = _frame(window)
    menu = U.build_menu(fr, menu_items(fr))
    try:
        if anchor is not None:
            U.popup_menu(anchor, menu, wx.Point(0, anchor.GetSize()[1]))
        else:
            U.popup_menu(fr, menu)
    finally:
        menu.Destroy()


_CLS = {}


def _chip_class():
    if "chip" in _CLS:
        return _CLS["chip"]
    import wx
    from unidec_theme import C, crisp
    Base = T._cls("Chip")

    class MethodChip(Base):
        """Method in the top bar: the method of the file shown; a click opens the Method menu."""
        TIP = "Click for the method presets"

        def __init__(self, parent, frame):
            Base.__init__(self, parent, "Method")
            self.frame = frame
            self.SetCursor(wx.Cursor(wx.CURSOR_HAND))
            self.Bind(wx.EVT_LEFT_UP, lambda e: show_menu(self.frame, self))
            self.SetValue("")

        def SetValue(self, v):
            v = (v or "").strip() or "none"
            self.full_value = v
            Base.SetValue(self, v if len(v) <= 22 else v[:21].rstrip() + "\u2026")  # (the top bar is full)
            self.SetToolTip("Method: %s\n\n%s" % (v, self.TIP))
            self.GetParent().Refresh()  # the chips left of it move when its width changes: the whole bar again

        def _fit(self):
            Base._fit(self)
            s = self.GetMinSize()
            s = wx.Size(s.width + self.FromDIP(12), s.height)
            self.SetMinSize(s)
            self.SetInitialSize(s)

        def _on_paint(self, e):
            dc = wx.AutoBufferedPaintDC(self)
            dc.SetBackground(wx.Brush(self.GetParent().GetBackgroundColour()))
            dc.Clear()
            gc = crisp(dc)
            w, h = self.GetClientSize()
            gc.SetBrush(wx.Brush(wx.Colour(C["accent_bg"])))
            gc.SetPen(wx.Pen(wx.Colour(C["line2"]), 1))
            gc.DrawRoundedRectangle(0.5, 0.5, w - 1, h - 1, (h - 1) / 2.0)
            gc.SetFont(self.GetFont(), wx.Colour(C["muted"]))
            t1 = self._label + "  "
            w1, th = gc.GetTextExtent(t1)
            x = self.FromDIP(11)
            gc.DrawText(t1, x, (h - th) / 2.0)
            gc.SetFont(self._bold, wx.Colour(C["accent_text"]))
            gc.DrawText(self._value, x + w1, (h - th) / 2.0)
            cx, cy, r = w - self.FromDIP(14), h / 2.0, self.FromDIP(7) / 2.0
            path = gc.CreatePath()
            path.MoveToPoint(cx - r, cy - r / 2.0)
            path.AddLineToPoint(cx + r, cy - r / 2.0)
            path.AddLineToPoint(cx, cy + r * 0.7)
            path.CloseSubpath()
            gc.SetPen(wx.TRANSPARENT_PEN)
            gc.SetBrush(wx.Brush(wx.Colour(C["muted"])))
            gc.FillPath(path)

    _CLS["chip"] = MethodChip
    return MethodChip


def add_chip(frame, toolbar):
    """The Method chip in the top bar of LCMS and HRMS Postrun (one chip value per file, as the other
    chips); nothing for other windows."""
    if kind_of(frame) is None:
        return None
    try:
        chip = _chip_class()(toolbar, frame)
    except Exception as ex:
        import traceback
        traceback.print_exc()
        _log("Method chip not made: %s" % ex)
        return None
    frame.chips["method"] = chip
    toolbar.add(chip, 3)
    return chip

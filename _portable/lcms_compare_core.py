# -*- coding: utf-8 -*-
"""Data side of the Compare view of LCMS Analysis: one trace per open file (PDA chromatogram at a
wavelength, PDA max plot, TIC, base peak chromatogram or mass chromatogram), processed the same way
for every file (smoothing, baseline, time window, alignment on a peak, scale) and placed on top of
each other. No wx here: the window (lcms_compare.py) draws what these functions return, and they
are tested without it (tests/check_compare.py). The calculations run in the C++ library (msengine,
src/compare.cpp) when it is there; the numpy code below is the fallback and the reference
(tests/check_cmp_cpp.py compares the two)."""
import csv
import math
import os
import re
import zipfile
from xml.sax.saxutils import escape

import numpy as np

SIGNALS = [("pda", "PDA at a wavelength"), ("max", "PDA max plot"), ("tic", "TIC"),
           ("bpc", "Base peak chromatogram"), ("xic", "Mass chromatogram (m/z)")]
BASELINES = [("none", "None"), ("offset", "Offset"), ("drift", "Drift"), ("rolling", "Rolling minimum")]
SCALES = [("abs", "Same scale"), ("max", "Each to its tallest peak (100 %)"),
          ("ref", "Each to its peak at the reference time (100 %)")]
LAYOUTS = [("stacked", "Stacked"), ("offset", "Offset"), ("overlay", "Overlay")]
YAXES = [("auto", "Automatic"), ("bar", "Scale bar"), ("axis", "Axis with values"), ("none", "None")]
COLOURS = [("palette", "Palette"), ("gradient", "Gradient"), ("black", "Black"), ("dark", "Dark blue")]
LABELS = [("right", "At the right end"), ("outside", "Right of the frame"), ("left", "At the left end"),
          ("legend", "Legend"), ("none", "None")]
LABEL_TEXT = [("sample", "Sample name"), ("file", "File name"), ("both", "Sample name and file name")]
RT_LABELS = [("none", "None"), ("main", "Tallest peak of each trace"), ("peaks", "Peaks above the threshold")]
COMMON_WL = [214, 220, 230, 254, 260, 280, 320, 350, 400]


def default_settings():
    return {"signal": "pda", "wl": None, "bw": 4.0, "own_wl": False, "polarity": "+", "mz": None, "mz_win": 0.5,
            "mz_ppm": False, "smooth": 0, "baseline": "none", "rolling_min": 1.0, "t0": None, "t1": None,
            "align": None, "align_win": 0.3, "layout": "stacked", "scale": "abs", "ref_t": None, "spacing": 100,
            "skew": 0, "reverse": False, "colours": "palette", "lw": 0.8, "fill": False, "labels": "right",
            "label_text": "sample", "rt_labels": "none", "rt_min": 10.0, "guides": [], "yaxis": "auto",
            "off_spacing": 15, "off_skew": 4, "integ_thr": 1.0}


def migrate(saved):
    """Settings saved by older versions in the keys of this one: the scale bar check box (3.3) is the
    y axis choice (unticked: the axis with its values)."""
    out = dict(saved)
    if "yaxis" not in out and "scalebar" in out:
        out["yaxis"] = "auto" if out["scalebar"] else "axis"
    out.pop("scalebar", None)
    return out


def offset_layout(s):
    """Stacked or offset: the traces are moved apart (overlay: all on one baseline)."""
    return s.get("layout", "stacked") in ("stacked", "offset")


def spacing_of(s):
    """Spacing (% of the tallest trace) of the layout chosen: the offset layout keeps its own (small)
    spacing and skew, so that switching between the layouts keeps both looks."""
    return float(s.get("off_spacing", 15) if s.get("layout") == "offset" else s.get("spacing", 100))


def skew_of(s):
    return float(s.get("off_skew", 4) if s.get("layout") == "offset" else s.get("skew", 0))


def yaxis_mode(s, n=2):
    """The y axis drawn: "bar" (a scale bar instead of the axis: stacked traces apart from each other,
    where the axis reads true for the lowest only), "axis" (with its values: overlay, offset, a single
    trace) or "none"."""
    k = s.get("yaxis", "auto")
    if k == "auto":
        return "bar" if (s.get("layout", "stacked") == "stacked" and n > 1) else "axis"
    return k if k in ("bar", "axis", "none") else "axis"


def _lib():
    """msengine_py when the C++ library is in use, else None (the numpy code runs)."""
    try:
        import msengine_py
        return msengine_py if msengine_py.available() else None
    except Exception:
        return None


class TraceError(Exception):
    """A file without the data asked for (the reason in plain words)."""


def _event_of(ms, polarity):
    """The scan event of a polarity ('+' or '-'), else None."""
    evs = getattr(ms, "events", None) or []
    for e, ev in enumerate(evs):
        if (ev.get("polarity") or "+") == polarity:
            return e
    return None


def signal_text(s, wl=None):
    """Short description of the signal: '254 nm (bandwidth 4 nm)', 'TIC (+)', 'm/z 445.1 ± 0.5 (+)'."""
    pol = "+" if s.get("polarity", "+") == "+" else "−"
    k = s.get("signal", "pda")
    if k == "pda":
        w = wl if wl is not None else s.get("wl")
        return ("%g nm" % w if w is not None else "PDA") + " (bandwidth %g nm)" % float(s.get("bw", 4.0))
    if k == "max":
        return "PDA max plot"
    if k == "tic":
        return "TIC (%s)" % pol
    if k == "bpc":
        return "base peak chromatogram (%s)" % pol
    if s.get("mz") is None:
        return "mass chromatogram"
    return "m/z %g ± %g%s (%s)" % (float(s["mz"]), float(s.get("mz_win", 0.5)),
                                       " ppm" if s.get("mz_ppm") else "", pol)


def raw_trace(pda, ms, s, own_wl=None):
    """(t, y, units, description) of one file for the settings s; own_wl: the wavelength set in the
    PDA view of that file (used with s['own_wl']). Raises TraceError with the reason when the file
    has no such data."""
    k = s.get("signal", "pda")
    if k in ("pda", "max"):
        if pda is None:
            raise TraceError("no PDA data")
        if k == "max":
            t, y = pda.max_plot()
            return np.asarray(t, float), np.asarray(y, float), "mAU", "PDA max plot"
        wl = own_wl if (s.get("own_wl") and own_wl is not None) else s.get("wl")
        if wl is None:
            raise TraceError("no wavelength chosen")
        w = np.asarray(pda.wavelengths, float)
        if not (w[0] - 2 <= float(wl) <= w[-1] + 2):
            raise TraceError("%g nm is outside its PDA range (%.0f to %.0f nm)" % (float(wl), w[0], w[-1]))
        bw = max(float(s.get("bw", 4.0) or 0.0), 0.0)
        t, y = pda.chromatogram(float(wl), bw)
        return np.asarray(t, float), np.asarray(y, float), "mAU", signal_text(s, float(wl))
    if ms is None:
        raise TraceError("no MS data")
    pol = s.get("polarity", "+")
    e = _event_of(ms, pol)
    if e is None:
        raise TraceError("no %s ion scans" % ("positive" if pol == "+" else "negative"))
    if k in ("tic", "bpc"):
        t, y = ms.chromatogram(e, k)
    else:
        if s.get("mz") is None:
            raise TraceError("no m/z given")
        mz = float(s["mz"])
        tol = mz * float(s.get("mz_win", 0.5)) * 1e-6 if s.get("mz_ppm") else float(s.get("mz_win", 0.5))
        t, y = ms.chromatogram(e, "xic", mz, abs(tol))
    return np.asarray(t, float), np.asarray(y, float), "counts", signal_text(s)


def subtract_blank(t, y, bt, by):
    """y minus a blank run (bt, by) at the same times (interpolated; nothing subtracted outside the
    blank's time range)."""
    t, y = np.asarray(t, float), np.asarray(y, float)
    bt, by = np.asarray(bt, float), np.asarray(by, float)
    E = _lib()
    if E is not None:
        r = E.cmp_subtract_blank(t, y, bt, by)
        if r is not None:
            return r
    if len(bt) < 2:
        return y
    b = np.interp(t, bt, by, left=np.nan, right=np.nan)
    return y - np.where(np.isfinite(b), b, 0.0)


def smooth(y, n):
    n = int(n or 0)
    if n < 3:
        return np.asarray(y, float)
    try:
        import lcms_integrate as LI
        return np.asarray(LI.smooth(y, n), float)
    except Exception:
        k = np.ones(n) / n
        return np.convolve(np.asarray(y, float), k, mode="same")


def find_apex(t, y, tc, win):
    """Time of the highest peak top of y within tc +- win, refined by a parabola through it and its
    neighbours; None without a peak there. A peak top is the highest point within about 0.03 min
    on both sides, not at the edge of the window (the rising flank of a peak just outside the
    window was taken for a peak and the trace moved by the whole window), and it rises above the
    lowest points of the window on both sides by 10 x the noise of the trace (a wiggle of the
    baseline is not a peak)."""
    t, y = np.asarray(t, float), np.asarray(y, float)
    E = _lib()
    if E is not None:
        r = E.cmp_find_apex(t, y, tc, win, 0)
        if r is not False:
            return r
    m = np.where((t >= tc - win) & (t <= tc + win))[0]
    if len(m) < 3:
        return None
    dt = float(np.median(np.diff(t))) if len(t) > 1 else 0.0
    k = max(2, int(round(0.03 / dt))) if dt > 0 else 2
    noise = 1.4826 * float(np.median(np.abs(np.diff(y)))) / math.sqrt(2.0) if len(y) > 2 else 0.0
    cand = []
    for i in m[1:-1]:
        lo, hi = max(0, i - k), min(len(y), i + k + 1)
        if y[i] >= np.max(y[lo:hi]) and y[i] > np.min(y[lo:hi]):
            rise = y[i] - max(float(np.min(y[m[0]:i + 1])), float(np.min(y[i:m[-1] + 1])))
            if rise > 10.0 * noise:
                cand.append(int(i))
    if not cand:
        return None
    i = max(cand, key=lambda j: y[j])
    if 0 < i < len(t) - 1:
        y0, y1, y2 = y[i - 1], y[i], y[i + 1]
        d = y0 - 2 * y1 + y2
        if d < 0:
            off = 0.5 * (y0 - y2) / d
            if abs(off) <= 1:
                return float(t[i] + off * (t[i + 1] - t[i - 1]) / 2.0)
    return float(t[i])


def lambda_max(pda, s, shift=0.0):
    """Absorption maximum (nm, above 210) of the tallest peak of a PDA run: the peak of the PDA trace shown
    (the max plot for the other signals) in the time window (without one: after the first 8 % of the run,
    the injection peak), its UV spectrum minus the spectrum at the lowest point of the max plot 0.15 to
    1 min before it. Returns (nm, time of the peak, background subtracted) or None."""
    wl0 = float(s["wl"]) if (s.get("signal") == "pda" and s.get("wl") is not None) else None
    bw = max(float(s.get("bw", 4.0)), 0.0)
    E = _lib()
    if E is not None and getattr(pda, "A", None) is not None:
        r = E.cmp_lambda_max(pda.A, pda.wavelengths, pda.times, getattr(pda, "interval_s", 0.0), wl0, bw,
                             s.get("smooth", 0), shift, s.get("t0"), s.get("t1"))
        if r is not False:
            return r
    sh = float(shift or 0.0)
    t, y = pda.chromatogram(wl0, bw) if wl0 is not None else pda.max_plot()
    t, y = np.asarray(t, float), smooth(np.asarray(y, float), max(5, int(s.get("smooth", 0) or 0)))
    m = window_mask(t + sh, s.get("t0"), s.get("t1"))
    if s.get("t0") is None and len(t):
        m &= t >= t[0] + 0.08 * (t[-1] - t[0])
    if not np.any(m):
        return None
    k = int(np.flatnonzero(m)[np.argmax(y[m])])
    ta = float(t[k])
    tm, ym = pda.max_plot()
    tm, ym = np.asarray(tm, float), np.asarray(ym, float)
    before = (tm >= ta - 1.0) & (tm <= ta - 0.15)
    bg = None
    if np.any(before):
        tb = float(tm[np.flatnonzero(before)[np.argmin(ym[before])]])
        bg = (tb - 0.02, tb + 0.02)
    try:
        wl, a = pda.spectrum(ta, None, bg=bg)
    except TypeError:
        wl, a = pda.spectrum(ta)
        bg = None
    wl, a = np.asarray(wl, float), np.asarray(a, float)
    ok = wl >= 210.0
    if np.sum(ok) < 3:
        return None
    ww, aa = wl[ok], a[ok]
    loc = [i for i in range(1, len(aa) - 1) if aa[i] >= aa[i - 1] and aa[i] > aa[i + 1]]
    j = max(loc, key=lambda i: aa[i]) if loc else int(np.argmax(aa))
    return float(ww[j]), ta, bg is not None


def window_mask(t, t0, t1):
    m = np.ones(len(t), bool)
    if t0 is not None:
        m &= t >= float(t0)
    if t1 is not None:
        m &= t <= float(t1)
    return m


def baseline(t, y, kind, rolling_min=1.0):
    """y minus a baseline: 'offset' (lowest point at zero), 'drift' (straight line between the
    medians of the first and last 3 % of the points), 'rolling' (moving minimum over rolling_min
    minutes, smoothed; below zero clipped)."""
    y = np.asarray(y, float)
    E = _lib()
    if E is not None and len(y) == len(t):
        r = E.cmp_baseline(np.asarray(t, float), y, kind, rolling_min)
        if r is not None:
            return r
    if kind == "offset" and len(y):
        return y - float(np.nanmin(y))
    if kind == "drift" and len(y) > 10:
        k = max(2, len(y) // 33)
        a, b = float(np.median(y[:k])), float(np.median(y[-k:]))
        ta, tb = float(np.median(t[:k])), float(np.median(t[-k:]))
        line = a + (b - a) * (np.asarray(t, float) - ta) / ((tb - ta) or 1.0)
        return y - line
    if kind == "rolling" and len(y) > 20:
        try:
            import ms_deconv
            out = ms_deconv.subtract_baseline(np.column_stack([t, y]), width=max(float(rolling_min), 0.01))
            return np.asarray(out[:, 1], float)
        except Exception:
            return y - float(np.nanmin(y))
    return y


def process(t, y, s, shift=0.0):
    """One trace processed as the settings say: smoothing, time shift (alignment and the shift of
    the file), time window, baseline, scale. Returns {"t", "y" (as drawn, before the stacking
    offset), "factor" (y drawn = factor x signal), "top" (highest point in the window), "ok",
    "ref_missing", "scaled" (True when the trace was really put in %: not when it has no top above 0, or
    no reference time is given)}."""
    E = _lib()
    if E is not None:
        r = E.cmp_process(np.asarray(t, float), np.asarray(y, float), s, shift)
        if r is not None:
            return _mark_scaled(r, s)
    y = smooth(y, s.get("smooth", 0))
    t = np.asarray(t, float) + float(shift or 0.0)
    m = window_mask(t, s.get("t0"), s.get("t1"))
    t, y = t[m], y[m]
    y = baseline(t, y, s.get("baseline", "none"), s.get("rolling_min", 1.0))
    factor = 1.0
    sc = s.get("scale", "abs")
    if sc == "max" and len(y):
        top = float(np.nanmax(y))
        factor = 100.0 / top if top > 0 else 1.0
    ref_missing = False
    if sc == "ref" and len(y) and s.get("ref_t") is not None:
        win = max(float(s.get("align_win", 0.3)), 0.02)
        mm = (t >= float(s["ref_t"]) - win) & (t <= float(s["ref_t"]) + win)
        top = float(np.nanmax(y[mm])) if np.any(mm) else 0.0
        factor = 100.0 / top if top > 0 else 1.0
        ref_missing = not top > 0  # no signal there: the trace cannot be scaled (it is not in %)
    y = y * factor
    top = float(np.nanmax(y)) if len(y) else 0.0
    return _mark_scaled({"t": t, "y": y, "factor": factor, "top": top, "ok": len(t) > 1,
                         "ref_missing": ref_missing}, s)


def _mark_scaled(r, s):
    """r["scaled"]: the trace was really put in % (its unit is then %). Not when its tallest point is not
    above 0 ("max": factor 1), nor without a reference time ("ref"). The same test for the library and the
    reference above: a factor is applied exactly when the top after it is 100 (above 0), or when a
    reference time has a signal there."""
    sc = s.get("scale", "abs")
    r["scaled"] = bool(len(r["t"])) and ((sc == "max" and r["top"] > 0) or
                                         (sc == "ref" and s.get("ref_t") is not None and not r["ref_missing"]))
    return r


def alignment_apexes(raws, target, win, smooth_n=0):
    """Time of the peak top of each trace within target +- win (find_apex), None where there is no
    peak. raws: list of (t, y)."""
    E = _lib()
    out = []
    for t, y in raws:
        r = E.cmp_find_apex(np.asarray(t, float), np.asarray(y, float), float(target), float(win), smooth_n) \
            if E is not None else False
        out.append(find_apex(t, smooth(y, smooth_n), float(target), float(win)) if r is False else r)
    return out


def alignment_shifts(raws, target, win, smooth_n=0):
    """Shift of each trace (minutes) that puts its peak top within target +- win on target; 0 for a
    trace without a peak there. raws: list of (t, y)."""
    return [0.0 if a is None else float(target) - a for a in alignment_apexes(raws, target, win, smooth_n)]


def nice(v):
    """1, 2 or 5 x a power of ten, at most v (for the scale bar)."""
    if not v or v <= 0 or not math.isfinite(v):
        return 0.0
    p = 10 ** math.floor(math.log10(v))
    for m in (5, 2, 1):
        if m * p <= v:
            return m * p
    return p


def bar_text(v):
    """Number of the scale bar: plain from 0.001 to 99999, else m x 10^e (matplotlib mathtext),
    e.g. '200', '0.002', '5 \u00d7 10$^{6}$' (instead of '5e+06')."""
    if not v or not math.isfinite(v):
        return "0"
    if 1e-3 <= abs(v) < 1e5:
        return "%g" % v
    e = int(math.floor(math.log10(abs(v)) + 1e-9))
    m = v / 10.0 ** e
    return ("%g \u00d7 10$^{%d}$" % (m, e)) if abs(m - 1.0) > 1e-9 else ("10$^{%d}$" % e)


def stack(procs, spacing=100.0, reverse=False, skew=0.0, layout="stacked"):
    """Offsets of the traces: [(dx, dy)] in the order of procs. Stacked (and offset): each trace dy
    above the next one, dy = spacing % of the tallest trace; the first trace on top (reverse: at the
    bottom). skew: each step also moves the trace right by skew % of the time span (a waterfall)."""
    n = len(procs)
    E = _lib()
    if E is not None and n:
        r = E.cmp_stack([p["top"] for p in procs], [float(p["t"][-1] - p["t"][0]) if p["ok"] else 0.0 for p in procs],
                        [p["ok"] for p in procs], spacing, reverse, skew, layout in ("stacked", "offset"))
        if r is not None:
            return r
    if layout not in ("stacked", "offset") or n == 0:
        return [(0.0, 0.0)] * n
    h = max([p["top"] for p in procs if p["ok"]] + [0.0])
    step = h * float(spacing) / 100.0
    spans = [float(p["t"][-1] - p["t"][0]) for p in procs if p["ok"]]
    dxs = (max(spans) if spans else 0.0) * float(skew) / 100.0
    out = []
    for i in range(n):
        level = i if reverse else (n - 1 - i)
        out.append((level * dxs, level * step))
    return out


def integrate(t, y, thr_pct=1.0, t0=None, t1=None):
    """Peaks of a processed trace (lcms_integrate.auto_integrate in the time window)."""
    try:
        import lcms_integrate as LI
        rng = None if (t0 is None and t1 is None) else (float(t0 if t0 is not None else t[0]),
                                                         float(t1 if t1 is not None else t[-1]))
        return LI.auto_integrate(t, y, threshold_pct=float(thr_pct), min_width_s=2.0, t_range=rng)
    except Exception:
        return []


def number(v):
    """A value for a table cell: whole numbers from 1000 up (counts: 1017342, not 1.017e+06), four
    significant digits below."""
    v = float(v)
    if not math.isfinite(v):
        return ""
    return "%.0f" % v if abs(v) >= 1000 else "%.4g" % v


def peak_rows(entries, guides=(), thr_pct=1.0):
    """One row per trace: tallest peak (RT, height and actual area in the signal's units, area %) and the area % of
    the peak at each guide time. entries: dicts with "label", "proc", "units" and "shift" (min,
    optional: the time shift of the trace, alignment included; the table gives the retention time
    of the run itself, as the labels in the plot do, while the guide lines are at the times drawn).
    Returns (header, rows, peaks per entry)."""
    units = sorted(set(en.get("units", "") for en in entries))
    unit = units[0] if len(units) == 1 and units[0] else ""
    hdr = ["Trace", "Main peak (min)", "Height (%s)" % unit if unit else "Height",
           "Area (%s·s)" % unit if unit else "Area", "Area %"] + ["%.2f min area %%" % g for g in guides]
    rows, allp = [], []
    E = _lib()
    for en in entries:
        pr = en["proc"]
        r = E.cmp_peaks(pr["t"], pr["y"], pr["factor"], thr_pct, list(guides)) if (E is not None and pr["ok"]) else None
        if r is not None:
            pk, mi, mpct, gpct = r
            mp = pk[mi] if mi >= 0 else None
            try:
                import lcms_integrate as LI
                pk = LI.finish(pk)  # area % of each peak, as integrate() gives them
            except Exception:
                pass
            allp.append(pk)
            if mp is not None:
                row = [en["label"], "%.3f" % (mp["rt"] - float(en.get("shift", 0.0) or 0.0)),
                       number(mp["height"]) if len(units) == 1 else "%s %s" % (number(mp["height"]), en.get("units", "")),
                       number(mp["area"]), "%.1f" % mpct if mpct == mpct else ""]
            else:
                row = [en["label"], "", "", "", ""]
            row += ["%.1f" % v if v == v else "" for v in gpct]
            rows.append(row)
            continue
        y = pr["y"] / (pr["factor"] or 1.0)  # the signal itself (before the scale of the plot)
        pk = integrate(pr["t"], y, thr_pct) if pr["ok"] else []
        allp.append(pk)
        tot = sum(p["area"] for p in pk) or 0.0
        if pk:
            mp = max(pk, key=lambda p: p["height"])
            r = [en["label"], "%.3f" % (mp["rt"] - float(en.get("shift", 0.0) or 0.0)),
                 number(mp["height"]) if len(units) == 1 else "%s %s" % (number(mp["height"]), en.get("units", "")),
                 number(mp["area"]), "%.1f" % (100.0 * mp["area"] / tot) if tot > 0 else ""]
        else:
            r = [en["label"], "", "", "", ""]
        for g in guides:
            hit = [p for p in pk if p["t0"] <= g <= p["t1"]]
            if not hit and pk:
                near = min(pk, key=lambda p: abs(p["rt"] - g))
                hit = [near] if abs(near["rt"] - g) <= 0.1 else []
            r.append("%.1f" % (100.0 * hit[0]["area"] / tot) if (hit and tot > 0) else "")
        rows.append(r)
    return hdr, rows, allp


def _region_peak(t, y, a, b):
    """Trapezoidal positive area above the endpoint baseline, with exact limits.

    Interpolate the boundary samples rather than snapping each run's selection
    to a different acquisition point. Time is in minutes; area is signal * s.
    """
    t, y = np.asarray(t, float), np.asarray(y, float)
    if len(t) < 2 or len(t) != len(y) or not np.isfinite(t).all() or not np.isfinite(y).all():
        return None
    if np.any(np.diff(t) <= 0):
        return None
    lo, hi = max(a, float(t[0])), min(b, float(t[-1]))
    if not (math.isfinite(lo) and math.isfinite(hi)) or hi <= lo:
        return None
    inside = (t > lo) & (t < hi)
    ts = np.concatenate(([lo], t[inside], [hi]))
    ys = np.interp(ts, t, y)
    base = ys[0] + (ys[-1] - ys[0]) * (ts - lo) / (hi - lo)
    delta = ys - base
    # Insert zero crossings before clipping so partial trapezoids are exact.
    cross = np.flatnonzero(delta[:-1] * delta[1:] < 0)
    if len(cross):
        zeros = ts[cross] - delta[cross] * np.diff(ts)[cross] / np.diff(delta)[cross]
        ts = np.sort(np.concatenate((ts, zeros)))
        ys = np.interp(ts, t, y)
        base = ys[0] + (ys[-1] - ys[0]) * (ts - lo) / (hi - lo)
    d = np.maximum(ys - base, 0.0)
    k = int(np.argmax(d))
    return {"t0": lo, "t1": hi, "rt": float(ts[k]), "height": float(d[k]),
            "area": float(np.sum((d[1:] + d[:-1]) * np.diff(ts)) * 30.0),
            "b0": float(base[0]), "b1": float(base[-1]),
            "t": ts, "y": ys, "baseline": base}


def region_areas(entries, t0, t1, reference=None, x_values=None):
    """Integrate one displayed time region in every processed trace.

    The returned areas use the processed signal before display normalisation, in
    signal-units x seconds.  Every trace is integrated over the same processed
    time t0 to t1 (alignment and time shifts are already in ``proc["t"]``);
    ``dx``, the waterfall skew of the drawing, only moves where the fill is
    drawn, never the integrated window.  ``reference`` is a row index;
    None makes the largest positive area 100 %.  ``x_values`` is optional and is
    kept numeric for an area-versus-variable plot and spreadsheet export.
    """
    a, b = sorted((float(t0), float(t1)))
    if not np.isfinite([a, b]).all() or b <= a:
        return []
    vals = list(x_values or [])
    if vals:  # (None: no X value for that trace)
        if len(vals) != len(entries):
            raise ValueError("Give one X value per plotted trace (%d values required)" % len(entries))
        try:
            vals = [None if v is None else float(v) for v in vals]
        except (TypeError, ValueError):
            raise ValueError("X values must be numbers") from None
        if not all(v is None or math.isfinite(v) for v in vals):
            raise ValueError("X values must be finite numbers")
    out = []
    E = _lib()
    for i, en in enumerate(entries):
        pr = en.get("proc") or {}
        dx = float(en.get("dx", 0.0) or 0.0)
        factor = float(pr.get("factor", 1.0) or 1.0)
        p = None
        if pr.get("ok") and len(pr.get("t", [])) == len(pr.get("y", [])):
            y = np.asarray(pr["y"], float) / factor
            p = E.cmp_region(pr["t"], y, a, b) if E is not None else False
            if p is False:
                p = _region_peak(pr["t"], y, a, b)
        xv = vals[i] if i < len(vals) else float(i + 1)
        # no data there (the region outside the trace, a value not finite, time not rising): area,
        # height and relative area None, not 0
        out.append({"index": i, "label": en.get("label", "Trace %d" % (i + 1)),
                    "units": en.get("units", ""), "x": xv, "dx": dx, "peak": p,
                    "area": float(p["area"]) if p is not None else None,
                    "height": float(p["height"]) if p is not None else None,
                    "rt": float(p["rt"]) if p is not None else None,
                    "t0": float(p["t0"]) if p is not None else a,
                    "t1": float(p["t1"]) if p is not None else b})
    ref = None
    try:
        j = int(reference) if reference is not None else -1
        if 0 <= j < len(out) and (out[j]["area"] or 0.0) > 0:
            ref = out[j]
    except (TypeError, ValueError):
        pass
    if reference is None and out:
        ref = max(out, key=lambda r: r["area"] or 0.0)
        if not (ref["area"] or 0.0) > 0:
            ref = None
    denom = ref["area"] if ref is not None else 0.0
    for r in out:
        r["relative"] = 100.0 * r["area"] / denom if (denom > 0 and r["area"] is not None) else None
        r["reference"] = bool(ref is r)
    return out


def split_x(text):
    """The X values typed in a field: words separated by commas, semicolons or new lines ([] when empty)."""
    text = (text or "").strip()
    return [w.strip() for w in re.split(r"\r\n|[,;\r\n]", text)] if text else []


def parse_x_values(text):
    """X values typed (split_x) as numbers; ValueError naming the first value that is not a finite number."""
    out = []
    for i, w in enumerate(split_x(text)):
        try:
            v = float(w)
        except ValueError:
            v = None
        if v is None or not math.isfinite(v):
            raise ValueError('X value %d ("%s") is not a number' % (i + 1, w) if w else "X value %d is empty" % (i + 1))
        out.append(v)
    return out


def _xlsx_col(n):
    s = ""
    n = int(n) + 1
    while n:
        n, r = divmod(n - 1, 26)
        s = chr(65 + r) + s
    return s


def _xlsx_text(v):
    # XML 1.0 rejects most control characters; instrument/sample names can
    # contain them when copied from vendor metadata.
    s = "" if v is None else str(v)
    s = "".join(c for c in s if c in "\t\n\r" or ord(c) >= 32)
    return escape(s)


def write_area_xlsx(path, header, rows, meta=()):
    """Write a dependency-free Excel workbook containing actual area results.

    Cells that are numbers remain numbers in Excel (rather than CSV text), so
    the sheet is immediately usable for kinetic plots and further calculations.
    """
    path = os.fspath(path)
    grid = [[str(k), v] for k, v in meta]
    if grid:
        grid.append([])
    header_row = len(grid) + 1
    grid.append(list(header))
    grid.extend([list(r) for r in rows])
    xml_rows = []
    for ri, row in enumerate(grid, 1):
        cells = []
        for ci, value in enumerate(row):
            ref = "%s%d" % (_xlsx_col(ci), ri)
            style = ' s="1"' if ri == header_row else ""
            if isinstance(value, (int, float, np.integer, np.floating)) and math.isfinite(float(value)):
                cells.append('<c r="%s"%s><v>%.15g</v></c>' % (ref, style, float(value)))
            else:
                cells.append('<c r="%s" t="inlineStr"%s><is><t xml:space="preserve">%s</t></is></c>' %
                             (ref, style, _xlsx_text(value)))
        xml_rows.append('<row r="%d">%s</row>' % (ri, "".join(cells)))
    ncol = max([len(r) for r in grid] + [1])
    nrow = max(len(grid), 1)
    widths = []
    for ci in range(ncol):
        width = max([len(str(r[ci])) for r in grid if ci < len(r)] + [8])
        widths.append('<col min="%d" max="%d" width="%.1f" customWidth="1"/>' %
                      (ci + 1, ci + 1, min(max(width + 2, 10), 42)))
    sheet = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
             '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
             '<dimension ref="A1:%s%d"/><sheetViews><sheetView workbookViewId="0">'
             '<pane ySplit="%d" topLeftCell="A%d" activePane="bottomLeft" state="frozen"/>'
             '</sheetView></sheetViews><cols>%s</cols><sheetData>%s</sheetData>'
             '<autoFilter ref="A%d:%s%d"/></worksheet>') % (
                 _xlsx_col(ncol - 1), nrow, header_row, header_row + 1, "".join(widths), "".join(xml_rows),
                 header_row, _xlsx_col(ncol - 1), nrow)
    content_types = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                     '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
                     '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
                     '<Default Extension="xml" ContentType="application/xml"/>'
                     '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
                     '<Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
                     '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
                     '</Types>')
    root_rels = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                 '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                 '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
                 '</Relationships>')
    workbook = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
                'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
                '<sheets><sheet name="Peak areas" sheetId="1" r:id="rId1"/></sheets></workbook>')
    wb_rels = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
               '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
               '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>'
               '<Relationship Id="rId2" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>'
               '</Relationships>')
    styles = ('<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
              '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
              '<fonts count="2"><font><sz val="11"/><name val="Calibri"/></font>'
              '<font><b/><sz val="11"/><name val="Calibri"/></font></fonts>'
              '<fills count="2"><fill><patternFill patternType="none"/></fill>'
              '<fill><patternFill patternType="gray125"/></fill></fills>'
              '<borders count="1"><border/></borders><cellStyleXfs count="1"><xf/></cellStyleXfs>'
              '<cellXfs count="2"><xf fontId="0" fillId="0" borderId="0" xfId="0"/>'
              '<xf fontId="1" fillId="0" borderId="0" xfId="0" applyFont="1"/></cellXfs>'
              '</styleSheet>')
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        for name, data in (("[Content_Types].xml", content_types), ("_rels/.rels", root_rels),
                           ("xl/workbook.xml", workbook), ("xl/_rels/workbook.xml.rels", wb_rels),
                           ("xl/styles.xml", styles), ("xl/worksheets/sheet1.xml", sheet)):
            z.writestr(name, data.encode("utf-8"))
    return path


def consistency_notes(infos):
    """Plain warnings when the runs differ: run length, method file. infos: dicts with "label",
    "t_end", "method"."""
    out = []
    ends = [(i["label"], i["t_end"]) for i in infos if i.get("t_end")]
    if len(ends) > 1:
        lo, hi = min(e for _, e in ends), max(e for _, e in ends)
        if hi > 1.1 * lo:
            out.append("The runs have different lengths (%.1f to %.1f min)." % (lo, hi))
    meths = sorted(set(i["method"] for i in infos if i.get("method")))
    if len(meths) > 1:
        out.append("Different methods: %s." % ", ".join(meths))
    return out


def write_traces_csv(path, entries):
    """Every trace as drawn (processed and scaled, without the stacking offset): a time and a value
    column per trace."""
    cols = []
    for en in entries:
        pr = en["proc"]
        unit = "%" if en.get("scaled") else en.get("units", "")
        cols.append(("%s time (min)" % en["label"], pr["t"]))
        cols.append(("%s (%s)" % (en["label"], unit), pr["y"]))
    n = max([len(c[1]) for c in cols] + [0])
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        w.writerow([c[0] for c in cols])
        fmts = ["%.5f", "%.9g"] * len(entries)  # times to 0.0006 s, values with every digit that counts
        for i in range(n):
            w.writerow([(f % c[1][i]) if i < len(c[1]) else "" for f, c in zip(fmts, cols)])
    return path


def write_table_csv(path, header, rows, meta=()):
    with open(path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.writer(fh)
        for k, v in meta:
            w.writerow([k, v])
        if meta:
            w.writerow([])
        w.writerow(header)
        for r in rows:
            w.writerow(r)
    return path


# ------------------------------------------------------------------ m/z tool of the view (3.33): texts only
def adduct_text(name):
    """An adduct name as shown: '[M-H]-' with typographic minus signs ('[M−H]−')."""
    return str(name or "").replace("-", "−")


def halogen_hint(ratio):
    """What an M+2 partner of this intensity ratio (M+2 : M) suggests: one Cl about 0.32, two Cl 0.64, one Br
    0.97, Br and Cl 1.3, two Br 1.95 (natural isotope abundances); '' outside these bands."""
    try:
        r = float(ratio)
    except (TypeError, ValueError):
        return ""
    for lo, hi, txt in ((0.24, 0.45, "one Cl"), (0.52, 0.78, "two Cl"), (0.8, 1.17, "one Br"),
                        (1.17, 1.5, "one Br and one Cl"), (1.65, 2.3, "two Br")):
        if lo <= r < hi:
            return txt
    return ""


def mass_suggestion(masses, more=2, info=None):
    """Text of the neutral masses suggested by ms_adducts.neutral_masses (a list, best first):
    'Suggested M ≈ 444.1 ([M+H]+ 445.1, [M+Na]+ 467.1, [M−H]− 443.1), or 466.1, or 888.2' (the ions of
    the best one, positive first; up to more other masses); '' when there is none. With an M+2 partner
    (Br, Cl): 'M+2 446.1 at 0.98 (one Br?)'. A mobile phase additive (TFA ...) is named as such. info
    (ms_adducts.neutral_masses_info): the base peaks the best mass does not explain are named, so that a
    suggestion resting on weaker ions is not taken for the main compound."""
    ok = []
    for m in masses or []:
        try:
            if isinstance(m, dict) and m.get("mass") is not None and math.isfinite(float(m["mass"])):
                ok.append(m)
        except (TypeError, ValueError):
            pass
    if not ok:
        return ""
    best = ok[0]
    ions = [i for i in (best.get("ions") or []) if isinstance(i, dict) and i.get("mz") is not None]
    ions = [i for i in ions if i.get("polarity", "+") != "-"] + [i for i in ions if i.get("polarity", "+") == "-"]
    if best.get("additive"):  # only the mobile phase found (it was put after every other candidate)
        txt = "Only %s of the mobile phase (M ≈ %.1f) found, no compound" % (best["additive"], float(best["mass"]))
    elif best.get("assumed") and len(ions) == 1:  # no pair of ions agrees: one ion read as the protonated or
        i = ions[0]                               # deprotonated molecule
        txt = "M ≈ %.1f if %.1f is %s (one ion, no adduct pair)" % (
            float(best["mass"]), float(i["mz"]), adduct_text(str(i.get("adduct") or "").replace("assumed ", "")))
    else:
        # a weak case (its ions a few % of the base peaks; score below 0.5) is called so
        try:
            weak = float(best.get("score", 1.0)) < 0.5
        except (TypeError, ValueError):
            weak = False
        txt = ("Weak suggestion: M ≈ %.1f" if weak else "Suggested M ≈ %.1f") % float(best["mass"])
        if ions:
            txt += " (%s)" % ", ".join("%s %.1f" % (adduct_text(i.get("adduct")), float(i["mz"])) for i in ions)
    try:
        if best.get("m2") is not None and math.isfinite(float(best["m2"])):
            hint = halogen_hint(best.get("m2_ratio"))
            txt += "; M+2 %.1f at %.2f%s" % (float(best["m2"]), float(best.get("m2_ratio") or 0.0),
                                             (" (%s?)" % hint) if hint else "")
    except (TypeError, ValueError):
        pass
    if isinstance(info, dict) and not best.get("additive"):
        base, expl = info.get("base") or {}, info.get("explained") or {}
        # a base peak that is an ion of the mobile phase (TFA ...) is explained by it
        mob = [(i.get("polarity", "+"), float(i["mz"])) for m in ok if m.get("additive")
               for i in (m.get("ions") or []) if isinstance(i, dict) and i.get("mz") is not None]
        miss = []
        for pol, name in (("+", "+"), ("-", "\u2212")):
            try:
                b = base.get(pol)
                if b is None or not math.isfinite(float(b)) or expl.get(pol, True):
                    continue
                if any(q == pol and abs(v - float(b)) <= 0.3 for q, v in mob):
                    continue
                miss.append("%.1f (%s)" % (float(b), name))
            except (TypeError, ValueError):
                pass
        if miss:
            txt += "; base peak%s %s not explained" % ("s" if len(miss) > 1 else "", " and ".join(miss))
    # other masses: not the mobile phase, at least half as well supported as the first
    try:
        s0 = float(best.get("score", 0.0))
    except (TypeError, ValueError):
        s0 = 0.0
    alt = []
    for m in ok[1:]:
        try:
            if not m.get("additive") and float(m.get("score", 0.0)) >= 0.5 * s0:
                alt.append(m)
        except (TypeError, ValueError):
            pass
    for m in alt[:max(int(more), 0)]:
        txt += ", or %.1f" % float(m["mass"])
    return txt

"""
Chromatogram peak integration for LCMS Postrun and HRMS Postrun (MS Analysis).

Automatic integration follows the usual LC approach:
  * baseline level from a rolling lower envelope (follows gradient drift),
  * peaks found on the smoothed signal above that envelope, with a minimum
    height (percent of the largest peak) and a minimum width,
  * each peak runs from where the signal returns to the baseline, or to the
    valley shared with a neighbouring peak,
  * touching peaks share one straight baseline with vertical drop lines at
    the valleys; an isolated peak has its own straight baseline.
Areas are in signal units x seconds, as in LabSolutions.
"""
import numpy as np


def _lib():
    """msengine_py when the C++ library is in use: smoothing and integration run there (no
    scipy import, which costs seconds at the first integration under Windows) and the Python
    code is the fallback (the reference of tests/check_p.py)."""
    try:
        import msengine_py
        return msengine_py if msengine_py.available() else None
    except Exception:
        return None


def smooth(y, points):
    """Savitzky-Golay smoothing (points = window length, 0 or 1 = off). A
    quadratic over 3 points passes through all of them (no smoothing at all),
    so 3 points use a straight line: the 3 point moving average."""
    y = np.asarray(y, dtype=float)
    points = int(points or 0)
    if points < 3 or len(y) < points + 2:
        return y
    if points % 2 == 0:
        points += 1
    E = _lib()
    if E is not None:
        out = E.smooth(y, points)
        if out is not None:
            return out
    try:
        from scipy.signal import savgol_filter
        return savgol_filter(y, points, 2 if points >= 5 else 1)
    except Exception:
        k = np.ones(points) / points
        return np.convolve(y, k, mode="same")


def _envelope(y, window):
    """Rolling minimum followed by a rolling mean: a lower envelope that
    tracks slow baseline drift but not peaks narrower than the window."""
    n = len(y)
    w = int(max(3, min(window, n)))
    try:
        from scipy.ndimage import minimum_filter1d, uniform_filter1d
        return uniform_filter1d(minimum_filter1d(y, w, mode="nearest"), w, mode="nearest")
    except Exception:
        pad = np.pad(y, (w // 2, w - w // 2 - 1), mode="edge")
        m = np.array([pad[i:i + w].min() for i in range(n)])
        return np.convolve(np.pad(m, (w // 2, w - w // 2 - 1), mode="edge"), np.ones(w) / w, mode="valid")


def _noise(y):
    d = np.diff(y)
    if len(d) < 4:
        return 0.0
    return 1.4826 * np.median(np.abs(d - np.median(d))) / np.sqrt(2)


def _make_peak(t, y, i0, i1, b0, b1):
    """Peak between indices i0..i1 with baseline values b0 (at i0), b1 (at i1)."""
    ts = t[i0:i1 + 1]
    ys = y[i0:i1 + 1]
    if len(ts) < 2:
        return None
    base = b0 + (b1 - b0) * (ts - ts[0]) / max(ts[-1] - ts[0], 1e-12)
    d = np.clip(ys - base, 0, None)
    k = int(np.argmax(d))
    tsec = ts * 60.0
    area = float(np.sum((d[1:] + d[:-1]) * np.diff(tsec)) / 2.0)
    return {"rt": float(ts[k]), "t0": float(ts[0]), "t1": float(ts[-1]), "i0": int(i0), "i1": int(i1),
            "height": float(d[k]), "area": area, "b0": float(b0), "b1": float(b1), "apex_y": float(ys[k])}


def auto_integrate(t, y, threshold_pct=2.0, min_width_s=2.0, smooth_points=0, t_range=None,
                   baseline_window_min=1.5):
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(t) < 8:
        return []
    E = _lib()
    if E is not None and len(t) == len(y):
        peaks = E.lc_integrate(t, y, threshold_pct, min_width_s, smooth_points, t_range, baseline_window_min)
        if peaks is not None:
            return finish(peaks)
    ys = smooth(y, smooth_points if smooth_points else 5)
    dt_s = np.median(np.diff(t)) * 60.0 if len(t) > 1 else 1.0
    env = _envelope(ys, baseline_window_min * 60.0 / max(dt_s, 1e-6))
    d = ys - env
    lo_i, hi_i = 0, len(t) - 1
    if t_range is not None:
        a, b = min(t_range), max(t_range)
        idx = np.where((t >= a) & (t <= b))[0]
        if len(idx) < 8:
            return []
        lo_i, hi_i = int(idx[0]), int(idx[-1])
    seg = d[lo_i:hi_i + 1]
    top = float(seg.max()) if len(seg) else 0.0
    if top <= 0:
        return []
    noise = _noise(y[lo_i:hi_i + 1])
    thr = max(top * float(threshold_pct) / 100.0, 5.0 * noise)
    try:
        from scipy.signal import find_peaks
        wmin = max(1.0, float(min_width_s) / dt_s)
        pk = find_peaks(seg, height=thr, prominence=thr, width=wmin)[0]
        pk = pk + lo_i
    except Exception:
        pk = np.array([i for i in range(lo_i + 1, hi_i) if d[i] >= thr and d[i] >= d[i - 1] and d[i] > d[i + 1]])
    if len(pk) == 0:
        return []
    stop_level = max(noise * 2.0, 0.0)
    bounds = []
    for p in pk:
        lim = max(stop_level, 0.005 * d[p])
        i = p
        while i > lo_i:
            if d[i] <= lim:
                break  # back at the baseline
            if d[i - 1] > d[i] and d[max(i - 3, lo_i)] > d[i] + noise:
                break  # valley before the next peak
            i -= 1
        j = p
        while j < hi_i:
            if d[j] <= lim:
                break
            if d[j + 1] > d[j] and d[min(j + 3, hi_i)] > d[j] + noise:
                break
            j += 1
        bounds.append([i, j])
    # neighbouring peaks: meet at the valley
    for n in range(1, len(bounds)):
        a, b = bounds[n - 1], bounds[n]
        if b[0] <= a[1]:
            v = pk[n - 1] + int(np.argmin(ys[pk[n - 1]:pk[n] + 1]))
            a[1] = v
            b[0] = v
    # clusters of touching peaks share a baseline (drop lines at valleys)
    peaks = []
    n = 0
    while n < len(bounds):
        m = n
        while m + 1 < len(bounds) and bounds[m + 1][0] <= bounds[m][1]:
            m += 1
        c0, c1 = bounds[n][0], bounds[m][1]
        y0, y1 = ys[c0], ys[c1]
        for q in range(n, m + 1):
            i0, i1 = bounds[q]
            if i1 - i0 < 2:
                continue
            f0 = (t[i0] - t[c0]) / max(t[c1] - t[c0], 1e-12)
            f1 = (t[i1] - t[c0]) / max(t[c1] - t[c0], 1e-12)
            pkd = _make_peak(t, y, i0, i1, y0 + (y1 - y0) * f0, y0 + (y1 - y0) * f1)
            if pkd is not None and pkd["height"] > 0:
                peaks.append(pkd)
        n = m + 1
    return finish(peaks)


def manual_peak(t, y, t0, t1):
    """Peak over a dragged range; the baseline joins the signal at both ends."""
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    E = _lib()
    if E is not None and len(t) == len(y):
        p = E.lc_manual(t, y, t0, t1)
        if p is not None:
            return p or None
    a, b = min(t0, t1), max(t0, t1)
    idx = np.where((t >= a) & (t <= b))[0]
    if len(idx) < 3:
        return None
    i0, i1 = int(idx[0]), int(idx[-1])
    ys = smooth(y, 5)
    return _make_peak(t, y, i0, i1, ys[i0], ys[i1])


def peak_at(t, y, x, search_s=4.0, baseline_window_min=1.5):
    """The peak whose top is nearest to time x (within +- search_s seconds);
    its limits are found as in the automatic integration."""
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(t)
    if n < 8:
        return None
    E = _lib()
    if E is not None and len(t) == len(y):
        p = E.lc_peak_at(t, y, x, search_s, baseline_window_min)
        if p is not None:
            return p or None
    ys = smooth(y, 5)
    dt_s = np.median(np.diff(t)) * 60.0
    env = _envelope(ys, baseline_window_min * 60.0 / max(dt_s, 1e-6))
    d = ys - env
    noise = _noise(y)
    i = int(np.argmin(np.abs(t - x)))
    w = max(2, int(round(search_s / max(dt_s, 1e-6))))
    lo, hi = max(0, i - w), min(n - 1, i + w)
    p = lo + int(np.argmax(d[lo:hi + 1]))
    spread = 1.4826 * float(np.median(np.abs(d - np.median(d))))
    if d[p] <= max(4.0 * noise, 3.0 * spread) or d[p] <= 0:
        return None  # nothing that stands out from the baseline noise
    lim = max(2.0 * noise, 0.005 * d[p])
    i0 = p
    while i0 > 0:
        if d[i0] <= lim:
            break
        if d[i0 - 1] > d[i0] and d[max(i0 - 3, 0)] > d[i0] + noise:
            break
        i0 -= 1
    i1 = p
    while i1 < n - 1:
        if d[i1] <= lim:
            break
        if d[i1 + 1] > d[i1] and d[min(i1 + 3, n - 1)] > d[i1] + noise:
            break
        i1 += 1
    if i1 - i0 < 2:
        return None
    return _make_peak(t, y, i0, i1, ys[i0], ys[i1])


def split_peak(t, y, p, x):
    """Split peak p with a vertical drop line at time x; both parts keep the
    original baseline."""
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    E = _lib()
    if E is not None and len(t) == len(y) and 0 <= p["i0"] <= p["i1"] < len(t):
        r = E.lc_split(t, y, p["i0"], p["i1"], p["b0"], p["b1"], x)
        if r is not None:
            return r or None
    k = int(np.argmin(np.abs(t - x)))
    if not (p["i0"] + 1 < k < p["i1"] - 1):
        return None
    ta, tb = t[p["i0"]], t[p["i1"]]
    bm = p["b0"] + (p["b1"] - p["b0"]) * (t[k] - ta) / max(tb - ta, 1e-12)
    a = _make_peak(t, y, p["i0"], k, p["b0"], bm)
    b = _make_peak(t, y, k, p["i1"], bm, p["b1"])
    if a is None or b is None:
        return None
    return a, b


def finish(peaks):
    """Sort by trace, then retention time; area % within each trace."""
    peaks = sorted(peaks, key=lambda p: (p.get("trace", ""), p["rt"]))
    totals = {}
    for p in peaks:
        totals[p.get("trace", "")] = totals.get(p.get("trace", ""), 0.0) + p["area"]
    for p in peaks:
        tot = totals[p.get("trace", "")]
        p["area_pct"] = 100.0 * p["area"] / tot if tot > 0 else 0.0
    return peaks


def write_table(path, peaks, units):
    import csv
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["Peak", "Trace", "RT (min)", "Start (min)", "End (min)", "Height (%s)" % units,
                    "Area (%s*s)" % units, "Area % (within trace)"])
        for n, p in enumerate(peaks, 1):
            w.writerow([n, p.get("trace", ""), "%.3f" % p["rt"], "%.3f" % p["t0"], "%.3f" % p["t1"],
                        "%.6g" % p["height"], "%.6g" % p["area"], "%.2f" % p["area_pct"]])

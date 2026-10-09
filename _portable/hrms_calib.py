"""
Internal m/z calibration for HRMS Analysis (MSpektra).

Reference lists (m/z of the singly charged ions, electron mass included):
  sodium formate clusters  [Na(HCOONa)n]+  and  [HCOO(HCOONa)n]-
  ESI-L low concentration tune mix (Agilent G1969-85000), positive and negative
  caesium iodide clusters  [Cs(CsI)n]+
  or any list typed in or loaded from a text file.

Models (measured m/z -> corrected m/z):
  Single point (ppm shift)     m' = k m
  Linear                       m' = a + b m
  Quadratic                    m' = a + b m + c m^2
  Cubic                        m' = a + b m + c m^2 + d m^3
  TOF (square root)            sqrt(m') = a + b sqrt(m) + c m
  HPC (quadratic + correction) quadratic, then a polynomial correction of
                               its remaining errors inside the calibrated range
The TOF model follows the time of flight relation (t proportional to
sqrt(m/z)) with one correction term; it extrapolates more gently than the
cubic polynomial outside the calibrated range.

HPC follows the idea of Bruker's High Precision Calibration ("Quadratic +
HPC" in DataAnalysis; Gobom et al., Anal. Chem. 2002, 74, 3915): after the
ordinary calibration the errors left at the calibrant peaks are systematic,
so a higher order polynomial is fitted to them and subtracted. Bruker does
not publish its exact algorithm; here the order of the correction is chosen
by leave-one-out cross-validation (each calibrant predicted from the others),
so it only goes up when it really predicts better, and outside the
calibrated range the correction is held at its value at the edge.

Every calibration also reports a cross-validated error: the RMS of each
calibrant predicted by a fit without it. It shows how well the model will do
for peaks between the calibrants; the ordinary RMS always improves with more
terms, the cross-validated one does not.
"""
import numpy as np

E = 0.00054858
NA = 22.98976928
H = 1.00782503
C = 12.0
O = 15.99491462
CS = 132.90545196
I = 126.904473
HCOONA = H + C + 2 * O + NA

CALIBRANTS = [
    ("Sodium formate, positive", "+", [NA - E + n * HCOONA for n in range(1, 25)]),
    ("Sodium formate, negative", "-", [H + C + 2 * O + E + n * HCOONA for n in range(1, 25)]),
    ("ESI-L tune mix, positive", "+", [118.086255, 322.048121, 622.028960, 922.009798, 1221.990637, 1521.971475,
                                       1821.952313, 2121.933152, 2421.913990, 2721.894829]),
    ("ESI-L tune mix, negative", "-", [112.985587, 301.998139, 601.978977, 1033.988109, 1333.968947, 1633.949786,
                                       1933.930624, 2233.911463, 2533.892301, 2833.873139]),
    ("Caesium iodide clusters, positive", "+", [CS - E + n * (CS + I) for n in range(0, 12)]),
    ("Custom list", "", []),
]

MODELS = ["Single point (ppm shift)", "Linear", "Quadratic", "Cubic", "TOF (square root)",
          "HPC (quadratic + correction)"]
N_PARAMS = {0: 1, 1: 2, 2: 3, 3: 4, 4: 3, 5: 5}
HPC_MAX_ORDER = 8
_C_MODEL = {0: 5, 1: 0, 2: 1, 3: 2, 4: 3, 5: 4}  # model numbers of the library (msengine.h)


def _lib():
    """msengine_py when the C++ library is in use: the searches and fits below run there and
    the Python code is the fallback (the reference of tests/check_p.py)."""
    try:
        import msengine_py
        return msengine_py if msengine_py.available() else None
    except Exception:
        return None


_SLICE = 2048


def _window(it, i0, i1, t):
    """Bounds (lo, hi) of a slice around the points i0 .. i1 - 1 that the apex fit of a peak there
    (its points of at least t on a falling flank, _apex_at) cannot reach the cut of, so that the
    library is given the slice and not a copy of a whole spectrum; None: give it all."""
    n = len(it)
    lo, hi = max(0, i0 - _SLICE), min(n, i1 + _SLICE)
    if lo > 0 and it[lo] >= t and not np.any(np.diff(it[lo:i0 + 1]) < 0):
        return None
    if hi < n and it[hi - 1] >= t and not np.any(np.diff(it[max(i1 - 1, 0):hi]) > 0):
        return None
    return lo, hi


def calibrant_names(polarity=None):
    return [n for n, p, _ in CALIBRANTS if not polarity or not p or p == polarity]


def calibrant_list(name):
    for n, p, v in CALIBRANTS:
        if n == name:
            return list(v)
    return []


def read_list(path):
    """m/z values from a text or csv file (first number on each line; columns
    separated by commas, semicolons, tabs or spaces)."""
    import re
    out = []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            for part in re.split(r"[,;\s]+", line.strip()):
                try:
                    out.append(float(part))
                    break
                except ValueError:
                    continue
    return sorted(v for v in out if v > 0)


def _peak_position(mz, it, lo, hi, profile):
    E = _lib()
    if E is not None:
        r = E.peak_position(mz, it, lo, hi, profile)
        if r is not None:
            return r
    m = (mz >= lo) & (mz <= hi)
    if not np.any(m):
        return None, 0.0
    idx = np.where(m)[0]
    k = idx[np.argmax(it[idx])]
    return _apex_at(mz, it, k, profile)


def _apex_at(mz, it, k, profile):
    """Position and height of the peak whose highest point is k."""
    top = float(it[k])
    if top <= 0:
        return None, 0.0
    if not profile:
        return float(mz[k]), top
    # Gaussian fit (parabola through the logarithm) of the points above 30 %
    # of the apex: unbiased for a symmetric peak, however it is sampled
    a = k
    while a > 0 and it[a - 1] >= 0.3 * top and it[a - 1] <= it[a]:
        a -= 1
    b = k
    while b < len(it) - 1 and it[b + 1] >= 0.3 * top and it[b + 1] <= it[b]:
        b += 1
    if b - a >= 2:
        x = mz[a:b + 1] - mz[k]
        y = np.log(np.maximum(it[a:b + 1], 1e-30))
        w = it[a:b + 1] ** 2  # Caruana weighting
        try:
            c2, c1, c0 = np.polyfit(x, y, 2, w=np.sqrt(w))
            if c2 < 0:
                pos = -c1 / (2 * c2)
                if abs(pos) <= (mz[b] - mz[a]):
                    return float(mz[k] + pos), top
        except Exception:
            pass
        ww = it[a:b + 1]
        return float(np.sum(mz[a:b + 1] * ww) / np.sum(ww)), top
    # too few points: Gaussian through the three top points
    if 0 < k < len(it) - 1 and it[k - 1] > 0 and it[k + 1] > 0:
        la, lb, lc = np.log(it[k - 1]), np.log(it[k]), np.log(it[k + 1])
        den = la - 2 * lb + lc
        if den < 0:
            off = 0.5 * (la - lc) / den
            return float(mz[k] + off * (mz[k + 1] - mz[k])), top
    return float(mz[k]), top


def profile_apex(mz, it, i):
    """m/z of the apex of the profile peak whose highest point is i (the same
    Gaussian fit as the calibrant search, the formula check and the clicked
    peaks, so that a label and a measured m/z of one peak agree)."""
    mz = np.asarray(mz, float)
    it = np.asarray(it, float)
    i = int(i)
    if not 0 <= i < len(mz):
        return float("nan")
    E = _lib()
    if E is not None:
        w = _window(it, i, i + 1, 0.3 * it[i])
        r = E.profile_apex(mz, it, i) if w is None else E.profile_apex(mz[w[0]:w[1]], it[w[0]:w[1]], i - w[0])
        if r is not None:
            return r
    m, _ = _apex_at(mz, it, i, True)
    return float(mz[i]) if m is None else float(m)


def peak_near(mz, it, target, d, profile=True, min_rel=0.05):
    """Position and height of the peak nearest to target within +- d: of the
    peak maxima in the window at least min_rel of the tallest there, the one
    closest to target (a taller neighbour a few ppm away does not take the
    place of the peak looked for). (None, 0.0) when there is no peak."""
    mz = np.asarray(mz, float)
    it = np.asarray(it, float)
    E = _lib()
    if E is not None:
        a, b = np.searchsorted(mz, target - d, "left"), np.searchsorted(mz, target + d, "right")
        if b <= a:
            return None, 0.0
        seg = it[a:b]
        w = _window(it, a, b, 0.3 * float(seg[seg > 0].min())) if np.any(seg > 0) else (max(0, a - 1), min(len(it), b + 1))
        r = E.peak_near(mz, it, target, d, profile, min_rel) if w is None else \
            E.peak_near(mz[w[0]:w[1]], it[w[0]:w[1]], target, d, profile, min_rel)
        if r is not None:
            return r
    a, b = np.searchsorted(mz, target - d, "left"), np.searchsorted(mz, target + d, "right")
    if b <= a:
        return None, 0.0
    if not profile:
        k = np.arange(a, b)
    else:
        lo, hi = max(a - 1, 0), min(b + 1, len(mz))
        y = it[lo:hi]
        if len(y) < 3:
            k = np.arange(a, b)
        else:
            k = np.where((y[1:-1] > y[:-2]) & (y[1:-1] >= y[2:]))[0] + 1 + lo
            k = k[(k >= a) & (k < b)]
            if not len(k):  # the window on the flank of a peak: its highest point
                k = np.array([a + int(np.argmax(it[a:b]))])
    k = k[it[k] > 0]
    if not len(k):
        return None, 0.0
    k = k[it[k] >= min_rel * it[k].max()]
    j = int(k[np.argmin(np.abs(mz[k] - target))])
    m, h = _apex_at(mz, it, j, profile)
    return (float(m), float(h)) if m is not None else (None, 0.0)


def find_calibrants(spec, refs, tol_ppm=30.0, min_rel=0.001, profile=True):
    """For each reference m/z the most intense peak within +- tol_ppm.
    Returns a list of dicts: ref, found (None if absent), intensity, rel,
    err (ppm, measured - reference), use."""
    spec = np.asarray(spec, float)
    mz, it = spec[:, 0], spec[:, 1]
    top = float(it.max()) if len(it) else 0.0
    out = []
    rs = [r for r in sorted(set(float(v) for v in refs)) if not (len(mz) == 0 or r < mz[0] or r > mz[-1])]
    E = _lib()
    res = E.find_calibrants(mz, it, rs, tol_ppm, min_rel, profile) if (E is not None and rs) else None
    for k, r in enumerate(rs if res is not None else ()):
        ok = bool(res[3][k])
        m, inten = float(res[0][k]), float(res[1][k])
        out.append({"ref": r, "found": m if ok else None, "intensity": inten if ok else 0.0,
                    "rel": float(res[2][k]) if ok else 0.0, "err": (m - r) / r * 1e6 if ok else None, "use": ok})
    for r in (rs if res is None else ()):
        d = r * tol_ppm * 1e-6
        m, inten = _peak_position(mz, it, r - d, r + d, profile)
        ok = m is not None and top > 0 and inten >= min_rel * top
        out.append({"ref": float(r), "found": m if ok else None, "intensity": inten if ok else 0.0,
                    "rel": 100.0 * inten / top if (ok and top) else 0.0,
                    "err": (m - r) / r * 1e6 if ok else None, "use": bool(ok)})
    # two reference ions closer than the search window can find the same peak: it
    # belongs to the nearer one (the other one would enter the fit with the
    # distance between them as its error)
    groups = {}
    for row in out:
        if row["found"] is not None:
            groups.setdefault(round(row["found"], 9), []).append(row)
    for rows in groups.values():
        rows.sort(key=lambda q: abs(q["err"]))
        for far in rows[1:]:
            far.update(found=None, intensity=0.0, rel=0.0, err=None, use=False)
    return out


class Calibration(object):
    """Fitted correction; call it with measured m/z (array) to get corrected m/z."""

    def __init__(self, model, measured, reference, cv=True, hpc_order=None):
        self.model = int(model)
        m = np.asarray(measured, float)
        r = np.asarray(reference, float)
        need = N_PARAMS[self.model]
        if len(m) < need:
            raise ValueError("%s needs at least %d calibrant peak%s (found %d)" % (
                MODELS[self.model], need, "s" if need > 1 else "", len(m)))
        self.lo, self.hi = float(m.min()), float(m.max())
        self.x0 = float(m.mean())
        self.s = float(max(m.max() - m.min(), 1.0))
        if self._fit_lib(m, r, cv, hpc_order, need):
            return
        # least squares on relative errors (ppm), so every calibrant counts alike
        if self.model == 0:
            q = m / r
            self.coef = np.array([float(np.sum(q) / np.sum(q * q))])
        elif self.model == 4:
            A = np.column_stack([np.ones_like(m), np.sqrt(m), m]) / np.sqrt(r)[:, None]
            self.coef = np.linalg.lstsq(A, np.ones_like(r), rcond=None)[0]
        elif self.model == 5:
            self.base = Calibration(2, m, r, cv=False)
            self.hpc_order, self.hpc = _hpc_fit(m, r, self.base, self.lo, self.hi, order=hpc_order)
            self.coef = np.concatenate([self.base.coef, self.hpc])
        else:
            x = (m - self.x0) / self.s
            A = np.column_stack([x ** k for k in range(self.model + 1)]) / r[:, None]
            self.coef = np.linalg.lstsq(A, (r - m) / r, rcond=None)[0]  # fit the correction, not m itself
        self.measured, self.reference = m, r
        self.after = (self(m) - r) / r * 1e6
        self.before = (m - r) / r * 1e6
        self.rms = float(np.sqrt(np.mean(self.after ** 2))) if len(m) else 0.0
        self.rms_before = float(np.sqrt(np.mean(self.before ** 2))) if len(m) else 0.0
        # cross-validated error: each calibrant predicted by the fit without it
        self.cv_rms, self.cv_err = None, None
        if cv and len(m) > need:
            errs = []
            for i in range(len(m)):
                keep = np.arange(len(m)) != i
                try:
                    c = Calibration(self.model, m[keep], r[keep], cv=False,
                                    hpc_order=getattr(self, "hpc_order", None))
                    errs.append((float(c(m[i:i + 1])[0]) - r[i]) / r[i] * 1e6)
                except (ValueError, np.linalg.LinAlgError):
                    errs.append(np.nan)
            errs = np.array(errs)
            if np.all(np.isfinite(errs)):
                self.cv_err = errs
                self.cv_rms = float(np.sqrt(np.mean(errs ** 2)))

    def _fit_lib(self, m, r, cv, hpc_order, need):
        """The fit, its errors and the cross validation in the library (False: not available)."""
        E = _lib()
        if E is None:
            return False
        res = E.calibration_fit(_C_MODEL[self.model], -1 if hpc_order is None else int(hpc_order), m, r)
        if res is None:
            return False
        coef, after, rms, cv_err, cv_rms = res
        if self.model == 5:
            self.base = Calibration(2, m, r, cv=False)
            self.hpc_order = len(coef) - 6
            self.hpc = np.array(coef[5:], float)
            self.coef = np.concatenate([self.base.coef, self.hpc])
        else:
            self.coef = np.array(coef[2:], float)
        self.measured, self.reference = m, r
        self.after = after
        self.before = (m - r) / r * 1e6
        self.rms = rms
        self.rms_before = float(np.sqrt(np.mean(self.before ** 2))) if len(m) else 0.0
        self.cv_rms, self.cv_err = (cv_rms, cv_err) if (cv and len(m) > need and cv_rms is not None) else (None, None)
        return True

    def _lib_coef(self):
        """Coefficients in the layout of the library: x0, s, then those of the model."""
        if self.model == 5:
            return np.concatenate([[self.base.x0, self.base.s], self.base.coef, self.hpc])
        return np.concatenate([[self.x0, self.s], self.coef])

    def __call__(self, mz):
        E = _lib()
        if E is not None:
            y = E.calibration_apply(_C_MODEL[self.model], self._lib_coef(), self.lo, self.hi, mz)
            if y is not None:
                return y
        mz = np.asarray(mz, float)
        if self.model == 0:
            return mz * self.coef[0]
        if self.model == 4:
            a, b, c = self.coef
            return (a + b * np.sqrt(np.maximum(mz, 0)) + c * mz) ** 2
        if self.model == 5:
            y = self.base(mz)
            return y * (1.0 + _hpc_eval(self.hpc, mz, self.lo, self.hi))
        x = (mz - self.x0) / self.s
        corr = np.zeros_like(mz)
        for k, cf in enumerate(self.coef):
            corr = corr + cf * x ** k
        return mz + corr

    def name(self):
        if self.model == 5:
            return "HPC (quadratic + order %d correction)" % self.hpc_order
        return MODELS[self.model]

    def short(self):
        return "HPC" if self.model == 5 else MODELS[self.model].split(" ")[0].lower()

    def describe(self):
        cv = ", cross-validated %.2f ppm" % self.cv_rms if self.cv_rms is not None else ""
        return "%s calibration: RMS error %.2f ppm%s (before %.2f ppm), %d calibrants, m/z %.2f to %.2f" % (
            self.name(), self.rms, cv, self.rms_before, len(self.measured), self.lo, self.hi)

    def report(self, title, calibrant, rows):
        lines = ["MSpektra, HRMS: internal m/z calibration", "Data: " + title,
                 "Calibrant: " + calibrant, "Model: " + self.name(),
                 "RMS error: %.3f ppm (before: %.3f ppm)" % (self.rms, self.rms_before),
                 "Cross-validated RMS error (each calibrant predicted without it): %s" % (
                     "%.3f ppm" % self.cv_rms if self.cv_rms is not None else "n/a"),
                 "Calibrated range: m/z %.4f to %.4f" % (self.lo, self.hi)]
        if self.model == 0:
            lines.append("m' = k m, k = %.10f (%.3f ppm)" % (self.coef[0], (self.coef[0] - 1) * 1e6))
        elif self.model == 4:
            lines.append("sqrt(m') = a + b sqrt(m) + c m;  a = %.10g, b = %.10g, c = %.10g" % tuple(self.coef))
        elif self.model == 5:
            b = self.base
            lines.append("quadratic: m1 = m + sum_k c_k x^k with x = (m - %.6f) / %.6f" % (b.x0, b.s))
            lines.append("c = " + ", ".join("%.10g" % c for c in b.coef))
            lines.append("HPC correction: m' = m1 (1 + sum_j h_j P_j(t)), P_j Legendre polynomials, "
                         "t = 2 (m - %.6f) / %.6f - 1 held within -1..1" % (self.lo, self.hi - self.lo))
            lines.append("h = " + ", ".join("%.10g" % c for c in self.hpc))
        else:
            lines.append("m' = m + sum_k c_k x^k with x = (m - %.6f) / %.6f" % (self.x0, self.s))
            lines.append("c = " + ", ".join("%.10g" % c for c in self.coef))
        lines += ["", "reference m/z\tmeasured m/z\trelative intensity %\terror before (ppm)\terror after (ppm)\tused"]
        for r in rows:
            if r["found"] is None:
                lines.append("%.6f\t\t\t\t\tnot found" % r["ref"])
                continue
            after = (float(self(np.array([r["found"]]))[0]) - r["ref"]) / r["ref"] * 1e6
            lines.append("%.6f\t%.6f\t%.2f\t%.3f\t%.3f\t%s" % (r["ref"], r["found"], r["rel"], r["err"], after,
                                                             "yes" if r["use"] else "no"))
        return "\n".join(lines) + "\n"


def _hpc_t(mz, lo, hi):
    return np.clip(2.0 * (np.asarray(mz, float) - lo) / max(hi - lo, 1e-9) - 1.0, -1.0, 1.0)


def _hpc_eval(c, mz, lo, hi):
    from numpy.polynomial import legendre
    return legendre.legval(_hpc_t(mz, lo, hi), c)


def _hpc_coef(m, r, base, lo, hi, p):
    from numpy.polynomial import legendre
    res = (r - base(m)) / base(m)  # relative error left by the quadratic
    return legendre.legfit(_hpc_t(m, lo, hi), res, p)


def hpc_max_order(n):
    """Highest HPC correction order for n calibrants: the quadratic and the
    correction together leave at least 4 calibrants of freedom (n - 8; with
    n - 4 the fit went through every calibrant, and the leave-one-out choice
    then often took order 6 or 7 on 12 calibrants, which follows the noise:
    larger errors between the calibrants than without the correction)."""
    return int(max(0, min(HPC_MAX_ORDER, n - 8)))


def _hpc_fit(m, r, base, lo, hi, order=None):
    """Order and coefficients of the HPC correction. Without a given order,
    the order with the smallest leave-one-out error is taken (a higher order
    only when it is at least 5 % better), so noise is not fitted."""
    n = len(m)
    top = hpc_max_order(n)
    if order is not None:
        p = int(max(0, min(order, top)))
        return p, _hpc_coef(m, r, base, lo, hi, p)
    best_p, best_e = 0, None
    for p in range(0, max(top, 0) + 1):
        errs = []
        for i in range(n):
            k = np.arange(n) != i
            b = Calibration(2, m[k], r[k], cv=False)
            c = _hpc_coef(m[k], r[k], b, float(m[k].min()), float(m[k].max()), p)
            pred = float(b(m[i:i + 1])[0]) * (1.0 + float(_hpc_eval(c, m[i:i + 1], float(m[k].min()),
                                                                     float(m[k].max()))[0]))
            errs.append((pred - r[i]) / r[i])
        e = float(np.sqrt(np.mean(np.square(errs))))
        if best_e is None or e < 0.95 * best_e:
            best_p, best_e = p, e
    return best_p, _hpc_coef(m, r, base, lo, hi, best_p)


def fit(rows, model):
    used = [r for r in rows if r["use"] and r["found"] is not None]
    return Calibration(model, [r["found"] for r in used], [r["ref"] for r in used])


def _segment_counts_xic(data, refs, event, tol_ppm):
    """Reference ions per scan from mass chromatograms (readers without the
    centroids of every scan in memory, e.g. the C++ core): an ion counts in a
    scan when its intensity within tol_ppm reaches 3 % of the scan's base
    peak; a scan counts only when these ions carry 15 % of its signal (the
    sum of all its centroids). None when the chromatograms would be slow."""
    if not getattr(data, "xic_fast", True) or not hasattr(data, "chromatogram"):
        return None
    idx = data.event_scans(event)
    ev = data.events[event] if hasattr(data, "events") else {}
    lo, hi = ev.get("mz_low") or 0.0, ev.get("mz_high") or float("inf")
    cal, data.calib = data.calib, None  # the m/z as recorded
    try:
        # all chromatograms in one pass over the scans when the library can (ms_xic_multi)
        sel = [float(r) for r in refs if lo <= r <= hi]
        E = _lib()
        Y = E.xic_multi(data, event, [1e6] + sel, [1e7] + [r * tol_ppm * 1e-6 for r in sel]) if E is not None else None
        if Y is not None:
            total = Y[0]
        else:
            _, total = data.chromatogram(event, "xic", 1e6, 1e7)  # every centroid of the scan
        total = np.asarray(total, float)
        base = np.asarray(data.bpc, float)[idx]
        counts = np.zeros(len(idx), int)
        share = np.zeros(len(idx))
        for j, r in enumerate(sel):
            if Y is not None:
                y = Y[1 + j]
            else:
                _, y = data.chromatogram(event, "xic", float(r), float(r) * tol_ppm * 1e-6)
            y = np.asarray(y, float)
            hit = (y > 0) & (y >= 0.03 * base)
            counts += hit
            share += np.where(hit, y, 0.0)
    finally:
        data.calib = cal
    counts[share < 0.15 * total] = 0
    return counts


def find_calibrant_segment(data, refs, event=0, tol_ppm=30.0, min_hits=5):
    """Time range (min) of the scans that contain the calibrant: in each scan
    of the event the reference ions present (centroids within tol_ppm, above
    3 % of the scan's base peak, carrying at least 15 % of its signal) are
    counted; the longest run of scans with at least half the best count is
    returned as (t0, t1, hits), or None."""
    refs = np.sort(np.asarray(refs, float))
    idx = data.event_scans(event)
    lines = getattr(data, "lines", None)
    if lines is None:
        counts = _segment_counts_xic(data, refs, event, tol_ppm)
        if counts is None:
            return None
    else:
        counts = np.zeros(len(idx), int)
    for n, i in enumerate(idx if lines is not None else ()):
        a, b = lines.offsets[i], lines.offsets[i + 1]
        m = lines.mz[a:b]
        it = lines.inten[a:b]
        if not len(m):
            continue
        keep = it >= 0.03 * it.max()
        order = np.argsort(m[keep])
        mk, ik = m[keep][order], it[keep][order]
        if not len(mk):
            continue
        pos = np.searchsorted(mk, refs)
        hit, share = 0, 0.0
        for r, p in zip(refs, pos):
            best = None
            for q in (p - 1, p):
                if 0 <= q < len(mk) and abs(mk[q] - r) <= r * tol_ppm * 1e-6:
                    if best is None or ik[q] > ik[best]:
                        best = q
            if best is not None:
                hit += 1
                share += ik[best]
        # the calibrant ions must carry a real part of the signal (not noise)
        if share < 0.15 * float(it.sum()):
            hit = 0
        counts[n] = hit
    if not len(counts) or counts.max() < min(min_hits, max(3, len(refs) // 3)):
        return None
    need = max(min(min_hits, max(3, len(refs) // 3)), int(np.ceil(0.5 * counts.max())))
    good = counts >= need
    best, cur, start, best_rng = 0, 0, 0, None
    for n, g in enumerate(list(good) + [False]):
        if g:
            if cur == 0:
                start = n
            cur += 1
        else:
            if cur > best:
                best, best_rng = cur, (start, n - 1)
            cur = 0
    if best_rng is None:
        return None
    a, b = best_rng
    t = data.rt[idx]
    return float(t[a]), float(t[b]), int(counts[a:b + 1].max())

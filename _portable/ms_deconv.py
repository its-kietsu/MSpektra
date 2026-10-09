"""
Deconvolution back ends for MS Analysis (LCMS Postrun and HRMS Postrun).

unidec_deconvolute   the UniDec engine (Marty et al., Anal. Chem. 2015,
                     DOI 10.1021/acs.analchem.5b00140), run in this process.
maxent_deconvolute   maximum entropy charge state deconvolution, written for
                     MS Analysis from the published principle (not a copy of
                     any vendor code): a positive mass x charge distribution
                     is fitted to the spectrum while its entropy relative to
                     a flat default model is kept as high as possible; the
                     weight of the entropy term is lowered step by step until
                     the fit reaches the noise level of the data.
isodec_deconvolute   isotope resolved deconvolution with UniDec's IsoDec
                     (monoisotopic masses from isotope patterns).

Each function returns a dict:
  method, mass (N x 2 zero charge mass spectrum), fit (N x 2, fitted m/z
  spectrum), data (N x 2, m/z data that were fitted), peaks (list of dicts:
  mass, height, area, z for IsoDec), r2, notes (text), zdist (charge
  distribution, if available).
"""
import os

import numpy as np

PROTON = 1.007276467


def carrier_mass(p):
    """Mass added per charge, with its sign (UniDec convention: m/z =
    (M + z * a) / z): +1.00728 for protonation, -1.00728 for the loss of a
    proton (negative ions), +22.9892 for Na+, +34.9694 for Cl-, ..."""
    a = p.get("adduct_mass")
    try:
        a = float(a) if a not in (None, "") else 0.0
    except (TypeError, ValueError):
        a = 0.0
    if a == 0.0:
        return PROTON * (-1 if (p.get("sign", 1) or 1) < 0 else 1)
    return a


def _adduct(x):
    """A carrier mass, or +1 / -1 (polarity) for the proton."""
    x = float(x)
    return PROTON * x if x in (1.0, -1.0) else x


def _restrict(spec, mz_range):
    spec = np.asarray(spec, dtype=float)
    if spec.ndim != 2 or spec.shape[1] != 2:
        raise ValueError("A spectrum must have two columns: m/z and intensity")
    E = _lib31()
    if E is not None:
        rng = mz_range and mz_range[0] is not None and mz_range[1] is not None and mz_range[1] > mz_range[0]
        r = E.restrict(spec, mz_range[0] if rng else 0.0, mz_range[1] if rng else 0.0)
        if r is not None:
            return r
    spec = spec[np.isfinite(spec).all(axis=1)]
    spec = spec[spec[:, 0] > 0]
    if np.any(spec[:, 1] < 0):
        raise ValueError("Spectrum intensities must not be negative")
    spec = spec[np.argsort(spec[:, 0])]
    if mz_range and mz_range[0] is not None and mz_range[1] is not None and mz_range[1] > mz_range[0]:
        m = (spec[:, 0] >= mz_range[0]) & (spec[:, 0] <= mz_range[1])
        spec = spec[m]
    return spec


def _lib31():
    """msengine_py when the C++ library is in use: the label, species and data helpers of this
    module (isotope_peaks, refine_peak_mass, refine_species, group_isotopes, peak_mask,
    subtract_baseline, _restrict) run there and the Python code is the fallback (the reference
    of tests/check_p.py); no scipy import then."""
    try:
        import msengine_py
        return msengine_py if msengine_py.available() else None
    except Exception:
        return None


def subtract_baseline(spec, width=None, factor=15.0):
    """Spectrum minus a smooth baseline under its peaks (morphological
    opening: moving minimum, then moving maximum, then smoothed), with a
    window of `factor` peak widths (or `width` in m/z). Removes broad humps
    of chemical noise that a deconvolution would otherwise try to explain
    with masses at low charges."""
    spec = np.asarray(spec, float)
    x, y = spec[:, 0], spec[:, 1]
    if len(x) < 20:
        return spec
    E = _lib31()
    if E is not None:
        b = E.subtract_baseline(spec, width, factor)
        if b is not None:
            out = spec.copy()
            out[:, 1] = b
            return out
    from scipy.ndimage import minimum_filter1d, maximum_filter1d, uniform_filter1d
    if not width:
        width = factor * estimate_peak_width(spec)
    dx = float(np.median(np.diff(x)))
    k = int(max(5, min(len(x) // 4, round(width / max(dx, 1e-12)))))
    base = minimum_filter1d(y, k, mode="nearest")
    base = maximum_filter1d(base, k, mode="nearest")
    base = uniform_filter1d(base, max(3, k // 2), mode="nearest")
    out = spec.copy()
    out[:, 1] = np.clip(y - np.minimum(base, y), 0, None)
    return out


def pick_mass_peaks(mass, window, threshold):
    """Local maxima of a mass spectrum: within +- window (Da), above
    threshold x the largest value; area = sum over +- window / 2."""
    x, y = mass[:, 0], mass[:, 1]
    if len(x) < 3 or y.max() <= 0:
        return []
    top = y.max()
    dx = np.median(np.diff(x)) if len(x) > 1 else 1.0
    w = max(1, int(round(window / max(dx, 1e-9))))
    peaks = []
    order = np.argsort(y)[::-1]
    taken = np.zeros(len(y), bool)
    for i in order:
        if y[i] <= 0 or y[i] < threshold * top:
            break
        if taken[i]:
            continue
        lo, hi = max(0, i - w), min(len(y), i + w + 1)
        if y[i] < y[lo:hi].max():
            continue
        taken[lo:hi] = True
        a0, a1 = max(0, i - w // 2), min(len(y), i + w // 2 + 1)
        # centroid of the top of the peak (above half height) for a precise mass
        seg = slice(a0, a1)
        yy = y[seg]
        m = yy >= 0.5 * y[i]
        cen = float(np.sum(x[seg][m] * yy[m]) / np.sum(yy[m])) if np.any(m) else float(x[i])
        peaks.append({"mass": cen, "apex": float(x[i]), "height": float(y[i]),
                      "area": float(np.sum(yy) * dx)})
    peaks.sort(key=lambda p: p["mass"])
    tot = sum(p["height"] for p in peaks) or 1.0
    for p in peaks:
        p["rel"] = 100.0 * p["height"] / max(q["height"] for q in peaks)
        p["frac"] = 100.0 * p["height"] / tot
    return peaks


# ============================================================================
# UniDec
# ============================================================================
# What the bundled engine (UniDec 8.2.1, unidec/src) can and cannot do, from
# its source, for the configuration written here; the same rules as in the
# C++ driver (cpp/msengine/src/unidec.cpp).
#
# The 0xC0000409 stops of the engine (2.3 to 2.7) all had one cause: its
# setup removes every mass x charge point that lacks an allowed neighbour at
# the next charge state (MakeSparseBlur), and when nothing is left it writes
# its outputs from an empty fit array (fwrite of a null buffer in WriteDecon),
# which Windows ends with 0xC0000409. This happens with charge smoothing below
# 1 (no neighbour is looked at), with a single charge state, and with data
# that have gaps (zero runs removed without data reduction, isolated peaks
# left by a minimum intensity) when no peak has the same mass at the next
# charge state. engine_points() predicts it; unidec_deconvolute then puts
# gapped data on bins of their own spacing or refuses with the reason.


def engine_exit_text(code):
    """Exit code of unidec.exe in plain words (Windows exception codes of a
    crash and the codes of the engine's own exit() calls). Every text starts
    with "The UniDec engine" (run_method: no second run in Python)."""
    try:
        code = int(code)
    except (TypeError, ValueError):
        return "The UniDec engine stopped (%s)." % code
    u = code & 0xFFFFFFFF
    w = " (Windows code 0x%08X" % u
    advice = "Use a larger mass step, narrower mass, charge or m/z ranges, or Data reduction: Automatic."
    if u == 0xC0000409:
        return ("The UniDec engine stopped itself" + w + ": one of its internal checks failed). With these settings it "
                "most likely found no data point it could use, or it could not get enough memory. Check the m/z, mass "
                "and charge ranges and the minimum intensity, or use Data reduction: Automatic.")
    if u == 0xC0000005:
        return ("The UniDec engine crashed" + w + ", memory access violation), most likely because it could not get the "
                "memory for this mass and charge grid. " + advice)
    if u in (0xC0000017, 0xC000012D):
        return "The UniDec engine ran out of memory" + w + "). " + advice
    if u == 0xC00000FD:
        return "The UniDec engine crashed" + w + ", stack overflow). " + advice
    if u == 0xC0000374:
        return "The UniDec engine crashed" + w + ", its memory was corrupted). " + advice
    if u == 0xC000013A:
        return "The UniDec engine was interrupted" + w + ")."
    if u == 1:
        return ("The UniDec engine could not read or write its files, or could not get memory (engine code 1). Check "
                "that the output folder can be written to and has free space.")
    if u == 2:
        return "The UniDec engine could not read a number in its input file (engine code 2)."
    if u == 10:
        return ("The UniDec engine found no data point that gives a mass inside the mass range at the chosen charges "
                "(engine code 10). Check the m/z, mass and charge ranges.")
    if u == 11:
        return "The UniDec engine could not get the memory for its peak shape table (engine code 11). " + advice
    if u in (14, 103):
        return "The UniDec engine stopped because the peak shape or peak width is not valid (engine code %d)." % u
    if u == 100:
        return "The UniDec engine stopped because a charge of 0 is in the charge range (engine code 100)."
    if u == 104:
        return "The UniDec engine stopped because two m/z values of its input are identical (engine code 104)."
    if -64 < code < 0:
        return "The UniDec engine was ended by signal %d." % -code
    if u >= 0xC0000000:
        return "The UniDec engine crashed" + w + "). " + advice
    return "The UniDec engine stopped with code %d." % code


def _engine_blur(zzsig, msig):
    """(zlength, mlength) of the engine's neighbourhood (SetUpBlur)."""
    zs, ms = np.float32(zzsig), np.float32(msig)
    if zs >= 0 and ms >= 0:
        return 1 + 2 * int(zs), 1 + 2 * int(ms)
    zl = 1 + 2 * int(3 * abs(zs) + np.float32(0.5)) if zs != 0 else 1
    ml = 1 + 2 * int(3 * abs(ms) + np.float32(0.5)) if ms != 0 else 1
    return zl, ml


def engine_memory(L, Z, M, C, W, transient):
    """Peak memory (bytes) of the engine for L input points, Z charges, M
    masses of its output axis, C neighbours per point and W points in the
    widest peak shape window (1 for binned data). From its allocations: mass
    table (float) and two allowed flags (char) per m/z x charge cell, the
    neighbour table (int and float per neighbour and cell), blur, newblur and
    oldblur (float per cell), the peak shape (float per point and window
    point), small per point arrays, and either the temporary copy of point
    smoothing or suppression (float per cell) or the mass x charge output
    (float per cell)."""
    cells = float(L) * float(Z)
    base = cells * (18.0 + 8.0 * C) + 4.0 * float(L) * max(float(W), 1.0) + 40.0 * float(L)
    outputs = 4.0 * float(M) * float(Z) + 16.0 * float(M)
    return base + max(outputs, 4.0 * cells if transient else 0.0) + 80e6  # + the program (70 MB measured)


def engine_index_problem(L, Z, M, C, W, speedy):
    """What exceeds the engine's 32 bit index arithmetic, or ""."""
    lim = 2147483647.0
    if C * float(L) * Z > lim:
        return "the engine's neighbour table (%d entries)" % int(C * float(L) * Z)
    if float(L) * Z > lim:
        return "the engine's m/z x charge grid"
    if float(M) * Z > lim:
        return "the engine's mass x charge grid"
    if not speedy and float(L) * max(float(W), 1.0) > lim:
        return "the engine's peak shape table"
    return ""


def engine_window_points(x, mzsig, psfun):
    """Widest peak shape window of the engine for data that are not binned
    (SetStartsEnds): from the point nearest to x - 6 sigma to the one nearest
    to x + 6 sigma (sigma = FWHM / 2.35482 for the Gaussian, the FWHM for the
    other shapes), reflected at the low end of the data as the engine does."""
    sig = np.float32(mzsig)
    if int(psfun) == 0:
        sig = np.float32(sig / np.float32(2.35482))
    d = np.asarray(x, float).astype(np.float32)
    if len(d) < 2:
        return 1
    win = np.float32(np.float32(6.0) * abs(sig))
    lo = (d - win).astype(np.float32)
    hi = (d + win).astype(np.float32)
    start = np.where(lo < d[0], -_nearfast(d, (np.float32(2.0) * d[0] - lo).astype(np.float32)), _nearfast(d, lo))
    end = np.where(hi > d[-1], len(d) - 1 + _nearfast(d, (np.float32(2.0) * d[0] - hi).astype(np.float32)),
                   _nearfast(d, hi))
    return int(max(1, int(np.max(end - start))))


def _nearfast(d, p):
    """The engine's nearest point (udtools.c nearfast) of each value of p in
    the sorted float32 array d: the nearer neighbour, the upper one on ties."""
    j = np.clip(np.searchsorted(d, p, side="left"), 1, len(d) - 1)
    lo, hi = d[j - 1], d[j]
    take_hi = np.abs(p - lo) >= np.abs(p - hi)
    return np.where(take_hi, j, j - 1)


def engine_points(x, y, startz, numz, adduct, masslb, massub, zzsig, msig, mzsig, psfun):
    """Which points of the mass x charge grid the engine keeps (SetLimits,
    MakeSparseBlur, KillB). Returns (in_range, sure, maybe): points of non
    zero intensity with an allowed mass at some charge; a point is certainly
    kept (it and its neighbour at z +- 1 have each other as neighbours, so both
    keep two neighbours whatever the order of the engine's threads); a point
    may be kept. The engine's float arithmetic is reproduced."""
    f = np.float32
    mz = np.asarray(x, float).astype(f)
    it = (np.round(np.asarray(y, float) * 1e6) / 1e6).astype(f)   # the "%.6f" text the engine reads
    n = len(mz)
    if n < 2:
        return 0, False, False
    zl, ml = _engine_blur(zzsig, msig)
    if f(zzsig) < 0 or f(msig) < 0:
        return 1, True, True   # not predicted
    a, lb, ub = f(adduct), f(masslb), f(massub)
    sig = f(mzsig)
    if int(psfun) == 0:
        sig = f(sig / f(2.35482))
    thr = f(sig * f(2.0))
    zh = (zl - 1) // 2
    zs = (np.arange(numz) + int(startz)).astype(f)
    lo_mz, hi_mz = f(mz[0] - thr), f(mz[-1] + thr)

    def allowed(m):
        return (m < ub) & (m > lb)

    def neighbour(idx, j, j2):
        """index of the neighbour of points idx at charge index j (array) for
        charge index j2 (array), -1 where none"""
        out = np.full(idx.shape, -1, np.int64)
        ok = (j2 >= 0) & (j2 < numz)
        z2 = np.where(ok, zs[np.clip(j2, 0, numz - 1)], f(1))
        ok &= z2 != 0
        m = zs[j] * (mz[idx] - a)
        point = (m + z2 * a) / z2
        ok &= ~((point < lo_mz) | (point > hi_mz))
        ind = _nearfast(mz, point)
        ok &= allowed(z2 * (mz[ind] - a)) & (np.abs(point - mz[ind]) < thr)
        out[ok] = ind[ok]
        return out

    nzp = np.where(it > 0)[0]
    in_range, maybe = 0, False
    for c0 in range(0, len(nzp), 20000):
        idx = nzp[c0:c0 + 20000]
        I = np.repeat(idx, numz)
        J = np.tile(np.arange(numz), len(idx))
        keep = allowed(zs[J] * (mz[I] - a))
        I, J = I[keep], J[keep]
        in_range += len(I)
        if not len(I):
            continue
        if ml >= 3:
            return in_range, True, True
        self_ok = neighbour(I, J, J) >= 0
        num = self_ok.astype(int)
        for dz in range(-zh, zh + 1):
            if dz == 0:
                continue
            K = neighbour(I, J, J + dz)
            has = K >= 0
            num += has
            if np.any(has & self_ok):
                sel = np.where(has & self_ok)[0]
                Kb, Jb = K[sel], J[sel] + dz
                back = neighbour(Kb, Jb, J[sel])
                kself = neighbour(Kb, Jb, Jb) >= 0
                if np.any((back == I[sel]) & kself):
                    return in_range, True, True
        if np.any(num >= 2):
            maybe = True
    return in_range, False, maybe


def check_grid(p, method):
    """Refuses settings that would run for hours or exhaust the memory (for
    example a mass step of 0.0001 Da over a wide mass range), and UniDec
    settings its engine cannot run. For UniDec this is the check before the
    data are known; unidec_deconvolute repeats it with the engine input."""
    lo, hi = float(p["mass_range"][0]), float(p["mass_range"][1])
    step = float(p.get("mass_step", 1.0))
    z0, z1 = map(float, p["z_range"])
    if not np.all(np.isfinite([lo, hi, step, z0, z1])):
        raise ValueError("Mass, charge and step values must be finite numbers")
    if lo <= 0 or z0 < 1 or z1 < z0 or z0 != int(z0) or z1 != int(z1):
        raise ValueError("Use positive masses and whole-number charges in increasing order")
    if step <= 0:
        raise ValueError("The mass step must be larger than 0")
    if step < 0.001:
        raise ValueError("A mass step of %g Da is finer than any mass spectrometer resolves (it would only make "
                         "the calculation extremely slow); use 0.005 Da or more" % step)
    if hi <= lo:
        raise ValueError("The upper mass must be larger than the lower mass")
    nm = (hi - lo) / step
    nz = int(p["z_range"][1]) - int(p["z_range"][0]) + 1
    if method == "unidec":
        if nz < 2:
            raise ValueError("UniDec needs at least two charge states (its engine stops with a single one, e.g. "
                             "1 to 1). Use a charge range such as 1 to 2, or Maximum entropy for one charge state.")
        # the engine stores masses as 32 bit floats
        ulp = float(np.ldexp(1.0, int(np.floor(np.log2(hi))) - 23))
        if step < 2.0 * ulp:
            use = next((v for v in (0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0) if v >= 2.0 * ulp), 1.0)
            raise ValueError("A mass step of %g Da is finer than the UniDec engine can represent at %g Da (it stores "
                             "masses as 32 bit numbers, %.2g Da apart there). Use a mass step of %g Da or more."
                             % (step, hi, ulp, use))
        zl, ml = _engine_blur(p.get("zzsig", 1.0) if p.get("zzsig") not in (None, "") else 1.0,
                              p.get("msig", 0.0) if p.get("msig") not in (None, "") else 0.0)
        if zl * ml < 2:
            raise ValueError("UniDec's engine needs a charge smoothing of 1 or more (or a mass smoothing of 1 Da or "
                             "more): with both below 1 it removes every point and stops (Windows code 0xC0000409). "
                             "Use Charge smoothing 1, the UniDec default.")
        # the mass x charge grid alone (float per cell, int index)
        budget = _memory_budget()
        cells = nm * nz
        need = 4.0 * cells + 16.0 * nm
        if cells > 2147483647.0 or need > budget:
            raise ValueError("Mass range %g to %g Da every %g Da with %d charges gives %.0f x %d grid points and would "
                             "need more than %s of memory%s. Use a larger mass step (isotope peaks need 0.01 to 0.02 "
                             "Da, envelopes 0.5 to 1 Da) or narrower mass and charge ranges."
                             % (lo, hi, step, nz, nm, nz, _fmt_bytes(need),
                                " (more than the UniDec engine can index)" if cells > 2147483647.0 else
                                "; this computer can give it " + _fmt_bytes(budget)))
        return int(nm), nz
    budget = _engine_budget()
    if budget and method == "maxent":
        # the C++ core: limited by the memory of this computer, not by a fixed
        # number of points (its own check repeats this with the data in hand)
        if nm > 2e9:
            raise ValueError("Mass range / mass step gives more than two billion mass points")
        need = 40.0 * nm * nz
        if need > budget:
            raise ValueError("Mass range %g to %g Da every %g Da with %d charges gives %.0f x %d grid points and would "
                             "need about %s of memory; this computer can give it %s. Use a larger mass step (isotope "
                             "peaks need 0.01 to 0.02 Da, envelopes 0.5 to 1 Da) or narrower mass and charge ranges."
                             % (lo, hi, step, nz, nm, nz, _fmt_bytes(need), _fmt_bytes(budget)))
        return int(nm), nz
    if method == "maxent":
        limit = 400000
        if nm > limit:
            raise ValueError("Mass range %g to %g Da every %g Da gives %d mass points; maximum entropy handles at most "
                             "%d. Use a larger mass step (isotope peaks need 0.01 to 0.02 Da) or a narrower mass range."
                             % (lo, hi, step, nm, limit))
    return int(nm), nz


def check_engine_input(n_points, numz, masslb, massub, massbins, endz, mzsig, psfun, zzsig, msig, psig, beta, speedy,
                       window_points, n_spectrum, mz_lo, mz_hi):
    """The whole memory need and the index limits of the engine for its
    input (unidec_deconvolute and the C++ driver, before the engine starts)."""
    zl, ml = _engine_blur(zzsig, msig)
    sig = np.float32(mzsig)
    if int(psfun) == 0:
        sig = np.float32(sig / np.float32(2.35482))
    M = (massub - masslb) / massbins + 2.0 * (6.0 * abs(float(sig)) * endz) / massbins + 2.0
    W = 1.0 if speedy else float(window_points)
    C = zl * ml
    idx = engine_index_problem(n_points, numz, M, C, W, speedy)
    need = engine_memory(n_points, numz, M, C, W, float(psig) >= 1 or float(beta) > 0) + 40.0 * n_spectrum
    budget = _memory_budget()
    if idx or need > budget:
        raise ValueError("UniDec would work on %.0f data points x %d charges (m/z %.2f to %.2f, %s) and %.0f masses: "
                         "its engine would need about %s of memory%s. Use Data reduction: Automatic (bins), a narrower "
                         "m/z or charge range, or a larger mass step."
                         % (n_points, numz, mz_lo, mz_hi, "binned" if speedy else "no data reduction", M,
                            _fmt_bytes(need), (" and " + idx + " would exceed its 32 bit index") if idx else
                            "; this computer can give it " + _fmt_bytes(budget)))


def _memory_budget():
    """Memory one deconvolution may use (bytes): half of the physical memory,
    at least 1 GB (the library's rule; from the operating system when the
    library is not loaded)."""
    b = _engine_budget()
    if b:
        return b
    total = 0.0
    try:
        if os.name == "nt":
            import ctypes

            class _MS(ctypes.Structure):
                _fields_ = [("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
                            ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
                            ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
                            ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
                            ("ullAvailExtendedVirtual", ctypes.c_ulonglong)]
            st = _MS()
            st.dwLength = ctypes.sizeof(_MS)
            if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(st)):
                total = float(st.ullTotalPhys)
        else:
            total = float(os.sysconf("SC_PAGE_SIZE")) * float(os.sysconf("SC_PHYS_PAGES"))
    except Exception:
        total = 0.0
    return max(0.5 * total, 1e9)


def _engine_budget():
    """Memory one deconvolution may use (bytes) when the C++ core is in use,
    else 0."""
    try:
        import msengine_py
        L = msengine_py.lib()
        if L is None or not hasattr(L, "ms_memory_budget"):
            return 0.0
        return float(L.ms_memory_budget())
    except Exception:
        return 0.0


def _fmt_bytes(b):
    return "%.1f GB" % (b / 1e9) if b >= 1e9 else "%.0f MB" % (b / 1e6)


def _engine_unique(data):
    """Merge points indistinguishable to UniDec's float32 reader, summing
    intensities as UniDec intends, without its histogram's last-bin loss."""
    data = np.asarray(data, float)
    if not len(data):
        return data.copy()
    keys = data[:, 0].astype(np.float32)
    if not np.isfinite(keys).all():
        raise ValueError("m/z values exceed the UniDec engine's numeric range")
    unique, inv, count = np.unique(keys, return_inverse=True, return_counts=True)
    if len(unique) == len(data) and np.all(np.diff(keys) > 0):
        return data.copy()
    mz = np.bincount(inv, weights=data[:, 0]) / count
    intensity = np.bincount(inv, weights=data[:, 1])
    return np.column_stack([mz, intensity])


def _export_engine_input(data, path):
    # Six decimals can move two distinct float32 values onto the same side
    # of a rounding boundary. Preserve double precision in the text handoff.
    data = _engine_unique(data)
    # Keep UniDec's established intensity export precision; only the m/z
    # column needs the higher precision to prevent error 104.
    np.savetxt(path, data, fmt=("%.17g", "%.6f"))
    return data


def _data_heights(spec, mass, zs, sign, width):
    """Heights in the spectrum at the m/z of each charge state of a mass."""
    a = _adduct(sign)
    x, y = spec[:, 0], spec[:, 1]
    out = []
    for z in zs:
        mz = (mass + z * a) / z
        lo, hi = np.searchsorted(x, [mz - width, mz + width])
        out.append(float(y[lo:hi].max()) if hi > lo else 0.0)
    return np.array(out)


# Windows cannot create a path longer than 260 characters (MAX_PATH); UniDec
# writes <folder>\\<name>_unidecfiles\\<name>_<suffix>, so a long data file name
# in a deep folder (e.g. OneDrive) fails. Longest suffix UniDec writes: _manualfile.dat
_PATH_LIMIT = 240


def _unidec_place(folder, name):
    """Folder and run name for UniDec whose files stay below the Windows path
    limit: the name as it is when it fits, else shortened with a short hash
    (so different spectra stay apart), else a folder in the temporary folder."""
    def longest(fo, na):
        return len(os.path.join(os.path.abspath(fo), na + "_unidecfiles", na + "_manualfile.dat"))
    if longest(folder, name) <= _PATH_LIMIT:
        return folder, name
    import hashlib
    tag = hashlib.md5(name.encode("utf-8")).hexdigest()[:6]
    for n in range(len(name) - 1, 7, -1):
        short = name[:n].rstrip(" ._-") + "_" + tag
        if longest(folder, short) <= _PATH_LIMIT:
            return folder, short
    import tempfile
    return os.path.join(tempfile.gettempdir(), "MS_Analysis_UniDec"), "run_" + tag


def _autocorr_from_lag2(datatop, config=None, window=10):
    """UniDec's tools.autocorr with two corrections (as the C++ library): the
    lags start at 2 points by index (the test corrx > xdiff let lag 1 through
    when the rounding made it a hair larger than xdiff), and the peaks of the
    autocorrelation must reach 1e-6 of its top (with 0 the round off of the FFT
    between the peaks of a clean spectrum was taken for peaks)."""
    from scipy import signal
    from unidec import tools as ud
    corry = signal.fftconvolve(datatop[:, 1], datatop[:, 1][::-1], mode='same')
    if np.amax(corry) == 0:
        return [[]], [[]]
    corry /= np.amax(corry)
    maxpos1 = np.argmax(datatop[:, 1])
    start = np.amax([maxpos1 - len(datatop) / 10, 0])
    end = np.amin([len(datatop) - 1, maxpos1 + len(datatop) / 10])
    cutdat = datatop[int(start):int(end)]
    if len(cutdat) < 20:
        cutdat = datatop
    xdiff = np.mean(cutdat[1:, 0] - cutdat[:len(cutdat) - 1, 0])
    corrx = np.arange(0.0, len(corry)) * xdiff
    maxpos = np.argmax(corry)
    corrx = corrx - corrx[maxpos]
    autocorr = np.transpose([corrx, corry])
    boo1 = np.arange(len(corry)) - maxpos >= 2
    cpeaks = ud.peakdetect(autocorr[boo1], config, window=window, threshold=1e-6)
    return autocorr, cpeaks


def _auto_peak_width(eng, ud, original_unique):
    """eng.get_auto_peak_width (UniDec's tools.auto_peak_width) measured on
    evenly spaced data: without bins (linflag 2) the data UniDec prepares have
    the zeros between the peaks removed, and the autocorrelation, whose lags
    are counted in points, then mixed gaps with steps (about 40 times too wide
    on clean spectra). Sets the peak width and shape of eng.config."""
    c = eng.config
    data = eng.data.data2
    if getattr(c, "linflag", 2) == 2:
        saved = c.linflag, c.mzbins
        c.linflag, c.mzbins = 0, 0
        ud.removeduplicates = _engine_unique
        try:
            data = ud.dataprep(eng.data.rawdata, c, silent=True)
        finally:
            c.linflag, c.mzbins = saved
            ud.removeduplicates = original_unique
    original_autocorr = ud.autocorr
    ud.autocorr = _autocorr_from_lag2
    try:
        fwhm, psfun, _mid = ud.auto_peak_width(np.asarray(data, float))
    finally:
        ud.autocorr = original_autocorr
    c.automzsig, c.autopsfun = fwhm, psfun
    c.psfun, c.mzsig = psfun, fwhm


def unidec_deconvolute(spec, folder, name, p, progress=None):
    """The UniDec engine; every setting of p that UniDec knows is passed on
    (numit, zzsig, psig, beta, psfun, msig, isotopemode, poolflag, smooth,
    subbuff, subtype, intthresh, linflag, mzbins, peaknorm; see the
    Deconvolute window). Intensities come back in the units of the data."""
    from unidec import engine
    check_grid(p, "unidec")
    spec = _restrict(spec, p.get("mz_range"))
    if len(spec) < 10:
        raise ValueError("Too few data points in the m/z range")
    folder, name = _unidec_place(folder, name)
    os.makedirs(folder, exist_ok=True)
    if progress:
        progress(1, 5, "preparing the data")
    fname = name + ".npz"
    np.savez(os.path.join(folder, fname), data=spec)  # binary: far faster than a text file
    udir = os.path.join(folder, name + "_unidecfiles")
    os.makedirs(udir, exist_ok=True)
    raw_copy = os.path.join(udir, name + "_rawdata.txt")
    if not os.path.isfile(raw_copy):
        open(raw_copy, "w").close()  # UniDec would otherwise write the whole spectrum as text again
    eng = engine.UniDec()
    eng.open_file(fname, folder, silent=True, refresh=True)
    c = eng.config
    c.minmz, c.maxmz = float(spec[0, 0]), float(spec[-1, 0])
    c.startz, c.endz = int(p["z_range"][0]), int(p["z_range"][1])
    c.numz = c.endz - c.startz + 1
    c.masslb, c.massub = float(p["mass_range"][0]), float(p["mass_range"][1])
    c.massbins = float(p.get("mass_step", 1.0))
    # a peak window below two mass steps makes every point of a noisy or
    # isotope resolved spectrum a peak of its own
    c.peakwindow = max(float(p.get("peak_window", 10.0)), 2.0 * float(p.get("mass_step", 1.0)))
    c.peakthresh = float(p.get("peak_thresh", 0.1))
    c.adductmass = carrier_mass(p)
    try:
        c.polarity = "Negative" if (p.get("sign", 1) or 1) < 0 else "Positive"
    except Exception:
        pass
    for key, attr, cast in (("numit", "numit", int), ("zzsig", "zzsig", float), ("psig", "psig", float),
                            ("beta", "beta", float), ("psfun", "psfun", int), ("msig", "msig", float),
                            ("isotopemode", "isotopemode", int), ("poolflag", "poolflag", int),
                            ("smooth", "smooth", float), ("subbuff", "subbuff", float), ("subtype", "subtype", int),
                            ("peaknorm", "peaknorm", int)):
        if p.get(key) is not None and p.get(key) != "":
            setattr(c, attr, cast(p[key]))
    binning = p.get("binning", "auto")
    if binning == "none":
        c.linflag, c.mzbins = 2, 0
    elif binning == "linear":
        c.linflag, c.mzbins = 0, float(p.get("mzbins") or 0.01)
    elif binning == "resolution":
        c.linflag = 1
        c.mzbins = float(p.get("mzbins") or 0) or float(spec[0, 0] / (3.0 * max(estimate_resolution(spec), 1000.0)))
    elif len(spec) > 60000:  # automatic
        if float(p.get("mass_step", 1.0)) >= 0.1:
            # charge state envelopes only: linear m/z bins keep the engine
            # fast and are still far narrower than an envelope
            c.linflag = 0
            c.mzbins = float(max(np.median(np.diff(spec[:, 0])), (spec[-1, 0] - spec[0, 0]) / 60000.0))
        else:
            # isotope resolved work: bins of constant resolving power (TOF),
            # three per isotope peak, so the isotopes stay resolved
            c.linflag = 1
            c.mzbins = float(spec[0, 0] / (3.0 * max(estimate_resolution(spec), 1000.0)))
    width = float(p.get("peak_width", 0) or 0)
    if width > 0:
        c.mzsig = width
    # UniDec's own intensity threshold acts on the normalised data (0 to 1),
    # not on the counts shown; the minimum intensity of the window is applied
    # before (run_method), so it stays off here
    c.intthresh = 0.0
    psfun = int(getattr(c, "psfun", 0) or 0)
    if progress:
        progress(2, 5, "processing the data")
    # The engine process owns this temporary patch; do not modify the
    # installed UniDec library or its deconvolution algorithm.
    from unidec import tools as ud
    original_unique = ud.removeduplicates
    ud.removeduplicates = _engine_unique
    try:
        eng.process_data(silent=True)
    finally:
        ud.removeduplicates = original_unique
    eng.data.data2 = _export_engine_input(eng.data.data2, c.infname)
    d2 = np.asarray(getattr(eng.data, "data2", []), float)
    if d2.ndim != 2 or len(d2) < 3 or not (d2[:, 1] > 0).any():
        raise ValueError("No data left to deconvolute in the m/z range (check the m/z range and the minimum "
                         "intensity)")
    if width <= 0:
        try:
            _auto_peak_width(eng, ud, original_unique)
        except Exception:
            c.mzsig = 0
        if getattr(c, "psfun", None) is None:
            c.psfun = psfun  # the automatic width failed and left no peak shape (crashes the settings export)
        # the automatic width can fail on weak or noisy spectra (it came out
        # negative for a negative mode spectrum, which stops the engine)
        if not np.isfinite(c.mzsig) or c.mzsig <= 0 or c.mzsig > 0.2 * (spec[-1, 0] - spec[0, 0]):
            c.mzsig = estimate_peak_width(spec)
    if getattr(c, "linflag", 2) == 0 and c.mzbins > 0 and c.mzsig < 2.5 * c.mzbins and not p.get("native_bins"):
        c.mzsig = 2.5 * c.mzbins  # a peak narrower than the bins crashes the engine
    width_x = 0.0  # m/z where the width shown in the notes applies (constant resolution bins)
    d2w = np.asarray(eng.data.data2, float)
    if getattr(c, "linflag", 2) == 1 and c.mzbins > 0 and len(d2w) > 1:
        # constant resolution bins: the engine's peak shape is a kernel in points
        # set at the first bin (MakePeakShape1D), so the width it applies grows
        # with m/z; the width measured (or set) at the tallest peak is converted
        # to the first point, at least two bins (2.7 forced two of the widest bins)
        width_x = float(d2w[int(np.argmax(d2w[:, 1])), 0])
        if width_x > 0:
            c.mzsig = c.mzsig * float(d2w[0, 0]) / width_x
        c.mzsig = max(c.mzsig, 2.0 * c.mzbins)
    # what the engine will do with this input: its memory, its index limits
    # and whether any point survives its neighbour rule (else it stops with
    # 0xC0000409); the same as the C++ driver
    plan_note = ""

    def check_engine():
        d2 = np.asarray(eng.data.data2, float)
        speedy = int(getattr(c, "linflag", 2)) != 2
        ps = int(getattr(c, "psfun", 0) or 0)
        wpts = 1 if speedy else engine_window_points(d2[:, 0], c.mzsig, ps)
        check_engine_input(len(d2), int(c.numz), float(c.masslb), float(c.massub), float(c.massbins), int(c.endz),
                           float(c.mzsig), ps, float(c.zzsig), float(c.msig), float(c.psig), float(c.beta), speedy, wpts,
                           len(spec), float(d2[0, 0]), float(d2[-1, 0]))
        return engine_points(d2[:, 0], d2[:, 1], int(c.startz), int(c.numz), float(c.adductmass), float(c.masslb),
                             float(c.massub), float(c.zzsig), float(c.msig), float(c.mzsig), ps)
    in_range, sure, _maybe = check_engine()
    if not sure and int(getattr(c, "linflag", 2)) == 2 and in_range > 0:
        # data with gaps (zero runs removed, or isolated peaks): the masses lack
        # the neighbouring charge state the engine requires and nothing may be
        # left; on bins of the data's own spacing every m/z has a point again
        # (zero intensity ones count as neighbours); the measured width is kept
        bw = float(np.median(np.diff(spec[:, 0])))
        c.linflag, c.mzbins = 0, bw
        ud.removeduplicates = _engine_unique
        try:
            eng.process_data(silent=True)
        finally:
            ud.removeduplicates = original_unique
        eng.data.data2 = _export_engine_input(eng.data.data2, c.infname)
        d2 = np.asarray(eng.data.data2, float)
        if d2.ndim != 2 or len(d2) < 3 or not (d2[:, 1] > 0).any():
            raise ValueError("No data left to deconvolute in the m/z range (check the m/z range and the minimum "
                             "intensity)")
        plan_note = ("data put on linear bins of %.4g m/z (their own spacing): without data reduction their gaps left "
                     "no mass with the neighbouring charge state the UniDec engine needs, and it would have stopped" % bw)
        in_range, sure, _maybe = check_engine()
    if not sure:
        d2 = np.asarray(eng.data.data2, float)
        if in_range == 0:
            raise ValueError("No data point in m/z %.2f to %.2f gives a mass between %g and %g Da at charges %d to %d. "
                             "Check the m/z, mass and charge ranges." % (d2[0, 0], d2[-1, 0], c.masslb, c.massub,
                                                                         c.startz, c.endz))
        raise ValueError("UniDec cannot deconvolute this spectrum with these settings: no peak in m/z %.2f to %.2f has "
                         "the same mass at the next higher or lower charge (masses %g to %g Da, charges %d to %d), so "
                         "the UniDec engine would remove every point and stop (Windows code 0xC0000409). The m/z range "
                         "probably holds no charge state series, only isolated peaks (for example singly charged ions, "
                         "or a minimum intensity above the protein signal). Check the ranges and the minimum intensity."
                         % (d2[0, 0], d2[-1, 0], c.masslb, c.massub, c.startz, c.endz))
    if progress:
        progress(3, 5, "deconvoluting (UniDec engine, all processor cores)")
    # the engine's own output (its peak scores) is kept in <name>_log.txt, as
    # the C++ driver does; UniDec's exe_call would discard it
    log = {"text": ""}
    original_call = ud.exe_call

    def exe_call(call, silent=False):
        import subprocess
        r = subprocess.run(call, shell=False, capture_output=True, text=True)
        log["text"] = (r.stdout or "") + (r.stderr or "")
        return r.returncode
    ud.exe_call = exe_call
    try:
        out = eng.run_unidec(silent=True)
    finally:
        ud.exe_call = original_call
    try:
        with open(c.outfname + "_log.txt", "w") as fh:
            fh.write(log["text"])
    except Exception:
        pass
    if out not in (0, None):
        raise RuntimeError(engine_exit_text(out))
    if progress:
        progress(4, 5, "picking peaks")
    # peak scores and UniScore as the engine computed them (UD_score.c), the
    # same numbers the C++ driver reads; UniDec's Python copy of the scores
    # is not run (slow on large grids, slightly different numbers)
    eng.pick_peaks(calc_dscore=False)
    dscores = engine_dscores(log["text"])
    uniscore = engine_uniscore(c.outfname + "_error.txt")
    mass = np.array(eng.data.massdat, float)
    data2 = np.asarray(eng.data.data2, float)
    fit = np.asarray(getattr(eng.data, "fitdat", []), float)
    zs = np.arange(c.startz, c.endz + 1)
    # intensities in data units: the tallest mass is given the summed heights
    # of its charge state peaks in the spectrum (as vendor software shows it)
    sign = c.adductmass  # carrier mass (charge states of the tallest mass)
    factor = 1.0
    zdist = None
    if len(mass) and mass[:, 1].max() > 0:
        top = float(mass[int(np.argmax(mass[:, 1])), 0])
        heights = _data_heights(spec, top, zs, sign, max(c.mzsig, 1e-4))
        zrel = None
        try:  # charge distribution of the fit: the mass x charge grid summed over the masses
            mg = np.asarray(getattr(eng.data, "massgrid", []), float)
            if mg.size == len(mass) * len(zs):
                zrel = mg.reshape(len(mass), len(zs)).sum(axis=0)
        except Exception:
            zrel = None
        use = heights > 0
        if zrel is not None and zrel.max() > 0:
            use &= zrel >= 0.05 * zrel.max()
        s_h = float(heights[use].sum())
        if s_h > 0:
            factor = s_h / float(mass[:, 1].max())
        if zrel is not None and zrel.max() > 0 and heights.max() > 0:
            # UniDec's charge distribution (all masses), on the intensity scale of the data
            zdist = np.column_stack([zs, zrel / zrel.max() * float(heights[use].max() if use.any() else heights.max())])
        else:
            zdist = np.column_stack([zs, np.where(use, heights, 0.0)])
    mass[:, 1] *= factor
    fitxy = None
    if fit.ndim == 1 and len(fit) == len(data2) and len(data2):
        ymax = float(spec[:, 1].max()) / max(float(data2[:, 1].max()), 1e-30)
        fitxy = np.column_stack([data2[:, 0], fit * ymax])  # back to the data scale
    peaks = []
    for pk in eng.pks.peaks:
        # UniDec gives the grid point; the apex between the points is more
        # precise (with 0.1 Da steps the grid alone is up to 4 ppm off at 12 kDa)
        peaks.append({"mass": refine_peak_mass(mass, float(pk.mass)), "height": float(pk.height) * factor,
                      "area": float(getattr(pk, "area", 0) or 0) * factor,
                      "score": match_dscore(dscores, float(pk.mass), c.massbins)})
    peaks.sort(key=lambda q: -q["height"])
    kept = []
    for q in peaks:  # one entry per apex
        if all(abs(q["mass"] - r["mass"]) > 0.3 * c.massbins for r in kept):
            kept.append(q)
    peaks = sorted(kept, key=lambda q: q["mass"])
    _rel_frac(peaks)
    pm = np.sort([q["mass"] for q in peaks])
    grouped_used = False
    if c.massbins <= 0.25 and len(pm) >= 5 and np.sum(np.abs(np.diff(pm) - 1.00235) < 0.06) >= 4:
        # resolved isotopes: one entry per species (average mass, most
        # abundant isotope), as for maximum entropy, instead of one per isotope
        rng = []
        grouped = group_isotopes(mass, float(p.get("peak_thresh", 0.1)), ranges=rng)
        if grouped:
            for q, (lo, hi) in zip(grouped, rng):
                q["score"] = species_dscore(dscores, lo, hi, c.massbins)
            peaks = grouped
            grouped_used = True
    # the engine of UniDec 8.2.1 has no isotope mode (the setting is read and
    # ignored: identical outputs for 0, 1 and 2), so the notes do not claim one
    iso_txt = ", isotope mode ignored (the UniDec 8.2.1 engine has none)" if int(getattr(c, "isotopemode", 0) or 0) else ""
    shown = c.mzsig * width_x / float(d2w[0, 0]) if width_x > 0 else c.mzsig
    wtxt = " at m/z %.0f (constant resolution)" % width_x if width_x > 0 else ""
    notes = "UniDec: charge %d to %d, mass %g to %g Da every %g Da, peak width %.3g m/z%s%s" % (
        c.startz, c.endz, c.masslb, c.massub, c.massbins, shown, wtxt, iso_txt)
    if plan_note:
        notes += "; " + plan_note
    r2 = float(getattr(c, "error", 0) or 0)
    z_main = int(zdist[int(np.argmax(zdist[:, 1])), 0]) if zdist is not None and len(zdist) and zdist[:, 1].max() > 0 else 0
    # a width set by hand (not the one run_method measured before the minimum intensity)
    width_set = 0.0 if p.get("auto_width_measured") else width
    width_data = estimate_peak_width(spec) if width_set > 0 else 0.0
    quality = unidec_quality(peaks, uniscore, r2, float(shown), float(c.massbins), z_main, (c.masslb, c.massub),
                             (c.startz, c.endz), resolved=grouped_used, width_set=width_set, width_data=width_data)
    if quality["level"] == "poor":
        notes = quality["text"] + "; " + notes
    iso_pk = isotope_peaks(mass) if c.massbins < 0.1 and len(mass) else None
    return {"method": "UniDec", "mass": _thin_mass_axis(mass), "fit": fitxy, "data": spec, "peaks": peaks,
            "r2": r2, "uniscore": uniscore, "quality": quality,
            "zdist": zdist, "notes": notes, "folder": udir,
            "isotope_peaks": iso_pk}


# Reliability of a UniDec result, from UniDec's own scores (Kostelic and
# Marty, Methods Mol. Biol. 2022, 2500, 159-180; UD_score.c of the engine):
# the DScore of a peak is the product of its uniqueness score (how much of
# the data at its charge states it explains alone), its peak shape score (the
# same mass peak at every charge state), its charge state score (a smooth,
# single charge state distribution) and its FWHM score (peaks separated by a
# dip to half height); the UniScore is R squared times the mean DScore
# weighted by the squared peak intensities. The two factors fail differently:
# artifacts (a few peaks read at many charge states) can fit the spectrum well
# (high R squared) while their DScores are low; real masses with unsuited
# settings keep good DScores while R squared is low. So the level is decided
# on the intensity weighted DScore of the listed masses and on R squared
# separately; the score shown is the UniScore. Thresholds from the test
# spectra (README of the tests folder): real proteins and peptides gave mean
# DScores of 0.33 to 0.64, artifact results 0.07 to 0.12.
QUALITY_D_POOR = 0.2       # intensity weighted DScore below this: poor
QUALITY_D_GOOD = 0.35      # ... at least this (and R squared below) for good
QUALITY_R2_POOR = 0.2      # R squared below this: poor
QUALITY_R2_GOOD = 0.6
QUALITY_D_LOW = 0.2        # a listed mass with a DScore below this is "low"


def unidec_quality(peaks, uniscore, r2, mzsig, step, z_main, mass_range, z_range, resolved=False, width_set=0.0,
                   width_data=0.0):
    """res["quality"] of a UniDec result: score 0 to 1 (the UniScore),
    level good, fair or poor, and one or two plain sentences on what is wrong
    and what to try. peaks: the listed masses (height, score = DScore);
    mzsig: the peak width used (m/z); step: mass step; z_main: the charge
    with the largest share of the result; resolved: the masses are species
    of grouped isotopes; width_set: a peak width set by hand (0: measured)
    and width_data: the width measured on the data."""
    u = min(max(float(uniscore), 0.0), 1.0)
    r2 = float(r2)
    if not peaks:
        return {"score": 0.0, "level": "poor",
                "text": "Unreliable result: UniDec found no mass above the peak threshold. Check the m/z, mass and "
                        "charge ranges and the minimum intensity."}
    h2 = [float(q["height"]) ** 2 for q in peaks]
    d = [float(q.get("score", 0) or 0) for q in peaks]
    wd = sum(a * b for a, b in zip(h2, d)) / max(sum(h2), 1e-300)
    top = max(peaks, key=lambda q: q["height"])
    n = len(peaks)
    n_low = sum(1 for q in peaks if float(q.get("score", 0) or 0) < QUALITY_D_LOW and q is not top)
    pct = 100.0 * max(r2, 0.0)
    if wd < QUALITY_D_POOR or r2 < QUALITY_R2_POOR:
        level = "poor"
    elif wd >= QUALITY_D_GOOD and r2 >= QUALITY_R2_GOOD:
        level = "good"
    else:
        level = "fair"
    # isotopes resolved in the data: the peak width is below 0.8 isotope spacings at the main charge
    iso_resolved = z_main > 0 and mzsig > 0 and mzsig < 0.8 * 1.00235 / z_main
    if level == "poor" and wd < QUALITY_D_POOR:
        if resolved and r2 >= QUALITY_R2_POOR and n <= 10:
            text = ("Unreliable result: UniDec gives the masses low scores (mean DScore %.2f) because their isotope "
                    "peaks overlap in the mass spectrum. Use a mass step of 0.5 to 1 Da for average masses, or a "
                    "smaller peak width (now %.3g m/z) if the isotopes are resolved in the data." % (wd, mzsig))
        elif n >= 10:
            text = ("Unreliable result: UniDec gives its %d masses low scores (mean DScore %.2f), the pattern of a few "
                    "peaks read at many charge states (harmonics), not of charge state series. The m/z range "
                    "probably holds no protein of %g to %g Da at charges %d to %d: check where its charge states are "
                    "and set the m/z, mass and charge ranges to match." % (n, wd, mass_range[0], mass_range[1],
                                                                           z_range[0], z_range[1]))
        else:
            text = ("Unreliable result: UniDec gives the masses low scores (mean DScore %.2f): they are not supported "
                    "by consistent charge state series. Check the m/z, mass and charge ranges, the minimum intensity "
                    "and the peak width." % wd)
    elif level == "poor":
        if width_set > 0 and width_data > 0 and width_set > 3.0 * width_data:
            hint = ("The peak width of %.3g m/z is much wider than the peaks of the data (about %.3g m/z): set it to "
                    "0 (measured)." % (width_set, width_data))
        elif width_set > 0 and width_data > 0 and width_set < width_data / 3.0:
            hint = ("The peak width of %.3g m/z is much narrower than the peaks of the data (about %.3g m/z): set it "
                    "to 0 (measured)." % (width_set, width_data))
        else:
            hint = "Most of the signal is not explained by masses in the range: check the mass, charge and m/z ranges."
        text = "Unreliable result: UniDec explains only %.0f %% of the spectrum (R squared %.2f). %s" % (pct, r2, hint)
    elif level == "fair" and r2 < QUALITY_R2_GOOD:
        if iso_resolved and step >= 0.1:
            hint = ("The isotope peaks are resolved in the data but a mass step of %g Da cannot follow them: for "
                    "isotope resolved masses use 0.01 to 0.05 Da over a narrow mass range." % step)
        else:
            hint = ("Part of the signal is not explained: widen the mass or charge range, or check the baseline and "
                    "other compounds in the m/z range.")
        text = ("The masses are supported by charge state series (mean DScore %.2f) but UniDec explains only %.0f %% "
                "of the spectrum (R squared %.2f). %s" % (wd, pct, r2, hint))
    elif level == "fair":
        text = ("UniDec's scores are moderate (mean DScore %.2f, R squared %.2f). Masses with a score below %.1f may be "
                "artifacts; narrower mass and charge ranges usually help." % (wd, r2, QUALITY_D_LOW))
    else:
        text = "UniDec's scores are good (UniScore %.2f, mean DScore %.2f, R squared %.2f)." % (u, wd, r2)
        if n_low:
            text += (" %d minor mass%s ha%s a low score (below %.1f) and may be an artifact%s."
                     % (n_low, "es" if n_low > 1 else "", "ve" if n_low > 1 else "s", QUALITY_D_LOW,
                        "s" if n_low > 1 else ""))
    return {"score": u, "level": level, "text": text}


def engine_dscores(log_text):
    """(mass, intensity, DScore) of the peaks the engine scored, from its lines
    "Peak: Mass: m Int: i DScore: d" (UD_score.c score_from_peaks)."""
    out = []
    pos = 0
    while True:
        pos = log_text.find("Peak: Mass:", pos)
        if pos < 0:
            break
        parts = log_text[pos:pos + 200].split()
        try:
            out.append((float(parts[2]), float(parts[4]), float(parts[6])))
        except (IndexError, ValueError):
            pass
        pos += 11
    return out


def engine_uniscore(error_file):
    """The engine's UniScore ("uniscore = x" in <name>_error.txt), 0 if absent."""
    try:
        for line in open(error_file):
            k, _, v = line.partition("=")
            if k.strip() == "uniscore":
                return float(v)
    except (OSError, ValueError):
        pass
    return 0.0


def match_dscore(dscores, mass, massbins):
    """DScore of the engine's peak at this grid mass (within half a mass step)."""
    for m, _i, d in dscores:
        if abs(m - mass) <= 0.5 * massbins:
            return float(d)
    return 0.0


def species_dscore(dscores, lo, hi, massbins):
    """DScore of a species of resolved isotopes: the mean of the DScores of
    the engine's peaks from its first to its last isotope peak, weighted by
    their intensities (0 when the engine scored none of them)."""
    sw = sd = 0.0
    for m, i, d in dscores:
        if lo - 0.5 * massbins <= m <= hi + 0.5 * massbins:
            sw += i
            sd += i * d
    return float(sd / sw) if sw > 0 else 0.0


def refine_peak_mass(mass, m0):
    """Apex of the peak of a mass spectrum near the grid point m0: Gaussian
    (parabola through the logarithm) of the points above half height, or of
    the three top points."""
    x, y = mass[:, 0], mass[:, 1]
    if len(x) < 3:
        return float(m0)
    E = _lib31()
    if E is not None:
        r = E.refine_peak_mass(mass, m0)
        if r is not None:
            return r
    k = int(np.argmin(np.abs(x - m0)))
    top = y[k]
    if top <= 0 or k == 0 or k == len(y) - 1 or max(y[k - 1], y[k + 1]) > top:
        return float(m0)  # not a maximum of the spectrum: left as it is
    a, b = k, k
    while a > 0 and y[a - 1] >= 0.5 * top and y[a - 1] <= y[a]:
        a -= 1
    while b < len(y) - 1 and y[b + 1] >= 0.5 * top and y[b + 1] <= y[b]:
        b += 1
    try:
        if b - a >= 2:
            xx = x[a:b + 1] - x[k]
            c2, c1, _c0 = np.polyfit(xx, np.log(np.maximum(y[a:b + 1], 1e-300)), 2, w=y[a:b + 1] / top)
            if c2 < 0 and abs(c1 / (2 * c2)) <= (x[b] - x[a]):
                return float(x[k] - c1 / (2 * c2))
        if y[k - 1] > 0 and y[k + 1] > 0:
            la, lb, lc = np.log(y[k - 1]), np.log(top), np.log(y[k + 1])
            den = la - 2 * lb + lc
            if den < 0:
                off = max(-0.5, min(0.5, 0.5 * (la - lc) / den))
                step = x[k + 1] - x[k] if off >= 0 else x[k] - x[k - 1]
                return float(x[k] + off * step)
    except Exception:
        pass
    return float(x[k])


def _rel_frac(peaks):
    top = max([q["height"] for q in peaks] + [1e-30])
    tot = sum(q["height"] for q in peaks) or 1.0
    for q in peaks:
        q["rel"] = 100.0 * q["height"] / top
        q["frac"] = 100.0 * q["height"] / tot


# ============================================================================
# maximum entropy
# ============================================================================
def averagine_sigma(mass):
    """Standard deviation (Da) of the isotope distribution of a protein of
    this mass (averagine, Senko 1995): 13C, 15N, 18O, 34S, 2H."""
    n = mass / 111.1254
    atoms = {"C": 4.9384 * n, "H": 7.7583 * n, "N": 1.3577 * n, "O": 1.4773 * n, "S": 0.0417 * n}
    var = 0.0
    for el, cnt in atoms.items():
        iso = {"C": [(0, 0.9893), (1.00336, 0.0107)], "H": [(0, 0.999885), (1.00628, 0.000115)],
               "N": [(0, 0.99636), (0.99704, 0.00364)], "O": [(0, 0.99757), (1.00422, 0.00038), (2.00425, 0.00205)],
               "S": [(0, 0.9499), (0.99939, 0.0075), (1.99580, 0.0425)]}[el]
        m1 = sum(d * w for d, w in iso)
        m2 = sum(d * d * w for d, w in iso)
        var += cnt * (m2 - m1 * m1)
    return float(np.sqrt(var))


def _raw_isotope_areas(data, masses, p):
    """Area of each isotope peak (masses: its neutral masses) in the m/z data,
    summed over the charge states of the deconvolution: at every charge z the
    points within half an isotope spacing of (M + z a) / z are integrated and
    the lowest point around the peak is taken as its floor."""
    x, y = np.asarray(data[:, 0], float), np.asarray(data[:, 1], float)
    a = carrier_mass(p)
    z0, z1 = int(p["z_range"][0]), int(p["z_range"][1])
    out = np.zeros(len(masses))
    for z in range(max(1, z0), max(z0, z1) + 1):
        hw = 0.45 / z
        for i, m in enumerate(masses):
            c = (m + z * a) / z
            lo, hi = np.searchsorted(x, [c - hw, c + hw])
            if hi - lo < 2:
                continue
            seg, xs = y[lo:hi], x[lo:hi]
            floor = float(seg.min())
            out[i] += float(np.sum(0.5 * (seg[1:] + seg[:-1] - 2 * floor) * np.diff(xs)))
    return out


def refine_species(mass, peaks, spacing=1.00235, data=None, p=None):
    """Most abundant isotope of isotope resolved species, from the m/z data.

    For a protein of about 12 kDa the isotopes +7 and +8 are nearly equal
    (averagine: 100 and 95 %), and the deconvolution can tip the balance by a
    percent or two: Protein POS.d, 3.17 to 3.31 min, UniDec gave 12228.245 at
    100 % and 12229.249 at 99.2 %, while the isotope peaks in the m/z data,
    summed over the 13 charge states, are 100 and 100.1 % (3.21 to 3.29 min:
    100 and 102 %). So the isotopes next to the deconvolution's most abundant
    one are measured in the m/z data (_raw_isotope_areas) and the tallest of
    them there becomes q["apex"] (its mass as deconvoluted); q["apex_dec"]
    keeps the deconvolution's choice. When the next isotope is at least 95 %
    as tall, it is given as q["apex2"] with its height in percent
    (q["apex2_rel"]), so both are named. Without the m/z data the heights of
    the deconvoluted isotope peaks are used."""
    if mass is None or len(mass) < 5 or not peaks:
        return
    x = np.asarray(mass[:, 0], float)
    y = np.asarray(mass[:, 1], float)
    mm = np.column_stack([x, y])
    E = _lib31()
    use = data is not None and p is not None and len(data) > 10 and p.get("z_range")
    if E is not None:
        xc, yc = np.ascontiguousarray(x), np.ascontiguousarray(y)
        dx = np.ascontiguousarray(np.asarray(data[:, 0], float)) if use else None
        dy = np.ascontiguousarray(np.asarray(data[:, 1], float)) if use else None
        a = carrier_mass(p) if use else 0.0
        z0, z1 = (int(p["z_range"][0]), int(p["z_range"][1])) if use else (0, 0)
    for q in peaks:
        try:
            if not q.get("n_iso") or q.get("apex") is None or q.get("apex_dec") is not None:
                continue
            r = E.refine_species_one(xc, yc, dx, dy, a, z0, z1, spacing, float(q["apex"])) if E is not None else None
            if r is not None:
                checked, apex, apex2, apex2_rel, source, best = r
                if checked:
                    q["apex_dec"] = float(q["apex"])
                    q["apex_source"] = source
                    if best != 0:
                        q["apex"] = apex
                    if apex2 is not None:
                        q["apex2"] = apex2
                        q["apex2_rel"] = apex2_rel
                continue
            ap = float(q["apex"])
            ks = list(range(-3, 4))
            pos, hdec = {}, {}
            for k in ks:
                c = ap + k * spacing
                lo, hi = np.searchsorted(x, [c - 0.3, c + 0.3])
                if hi > lo:
                    t = lo + int(np.argmax(y[lo:hi]))
                    hdec[k] = float(y[t])
                    pos[k] = ap if k == 0 else refine_peak_mass(mm, float(x[t]))
            if not hdec.get(0, 0.0) > 0:
                continue
            ks = [k for k in ks if k in pos]
            heights, source = dict(hdec), "deconvolution"
            if data is not None and p is not None and len(data) > 10 and p.get("z_range"):
                areas = _raw_isotope_areas(data, np.array([pos[k] for k in ks]), p)
                if areas.max() > 0 and areas[ks.index(0)] > 0:
                    heights, source = {k: float(v) for k, v in zip(ks, areas)}, "m/z data"
            # the tallest of the deconvolution's apex and its two neighbours
            best = max((k for k in (-1, 0, 1) if k in heights), key=lambda k: heights[k])
            q["apex_dec"] = ap
            q["apex_source"] = source
            if best != 0:
                q["apex"] = float(pos[best])
            hb = heights[best]
            cand = [k for k in (best - 1, best + 1) if k in heights and heights[k] >= 0.95 * hb]
            if cand:
                k2 = max(cand, key=lambda k: heights[k])
                q["apex2"] = float(pos[k2])
                q["apex2_rel"] = 100.0 * heights[k2] / hb
        except Exception as ex:  # never lose a result over this
            print("[species] most abundant isotope not checked: %s" % ex)

def estimate_resolution(spec, n=5):
    """Resolving power m/FWHM: median over the tallest well separated peaks
    (half height crossings interpolated between the data points; taking the
    first point below half height would make the peaks look up to a third
    wider on a coarsely sampled TOF spectrum)."""
    x, y = np.asarray(spec[:, 0], float), np.asarray(spec[:, 1], float)
    if len(x) < 5 or y.max() <= 0:
        return 10000.0
    order = np.argsort(y)[::-1]
    used, vals = [], []
    for i in order[:5000]:
        if y[i] < 0.05 * y[order[0]] or len(vals) >= n:
            break
        if any(abs(x[i] - u) < 0.2 for u in used):
            continue
        used.append(x[i])
        w = _fwhm_at(x, y, int(i))
        if w and w > 0:
            vals.append(x[i] / w)
    return float(np.median(vals)) if vals else float(x[order[0]] / max(estimate_peak_width(spec), 1e-9))


def _fwhm_at(x, y, i):
    half = y[i] / 2.0
    lo, hi = i, i
    while lo > 0 and y[lo] > half:
        lo -= 1
    while hi < len(y) - 1 and y[hi] > half:
        hi += 1
    if y[lo] > half or y[hi] > half or lo == i or hi == i:
        return None
    xl = x[lo] + (half - y[lo]) * (x[lo + 1] - x[lo]) / max(y[lo + 1] - y[lo], 1e-30)
    xr = x[hi - 1] + (y[hi - 1] - half) * (x[hi] - x[hi - 1]) / max(y[hi - 1] - y[hi], 1e-30)
    return float(xr - xl)


def maxent_deconvolute(spec, p, progress=None):
    """p: z_range, mass_range, mass_step, peak_width (m/z, 0 = auto), sign,
    peak_window, peak_thresh, and for high resolution data:
      resolution  resolving power R (0 = measured); the m/z axis is then
                  handled in ln(m/z), so the peak width grows with m/z as on
                  a TOF instrument
      isotopes    "envelope": one peak per mass with the width of its
                  isotope distribution (average masses, as MaxEnt at low
                  resolving power); "resolved": every isotope peak (use a mass
                  step of 0.02 to 0.05 Da and a narrow mass range)
      threads     worker threads (0 or absent: all processor cores)."""
    import msengine_py
    return msengine_py.maxent(spec, p, progress)


def _thin_mass_axis(mass, limit=120000, floor=1e-5):
    """Mass spectrum with fewer points for the window and the files: the
    full sampling is kept around the peaks (above `floor` of the tallest
    point, widened by a margin), the empty stretches between them keep
    every k-th point (and stay at their values, i.e. near zero)."""
    x, y = mass[:, 0], mass[:, 1]
    if len(x) <= limit or y.max() <= 0:
        return mass
    keep = y >= floor * y.max()
    r = max(2, len(x) // 20000)
    cs = np.concatenate([[0], np.cumsum(keep)])  # dilation by r points: any kept point within r
    keep = (cs[np.minimum(np.arange(len(x)) + r + 1, len(x))] - cs[np.maximum(np.arange(len(x)) - r, 0)]) > 0
    k = int(np.ceil(len(x) / float(limit)))
    keep[::k] = True
    keep[0] = keep[-1] = True
    return mass[keep]


def _find_peaks(y, height, distance):
    """Local maxima of y at least `height` high and `distance` points apart
    (the higher one wins), as scipy.signal.find_peaks; numpy only, since
    importing scipy.signal costs seconds in the worker process."""
    y = np.asarray(y, float)
    if len(y) < 3:
        return np.zeros(0, int)
    cand = np.where((y[1:-1] > y[:-2]) & (y[1:-1] >= y[2:]) & (y[1:-1] >= height))[0] + 1
    # a flat top counts once: scipy takes the middle of a plateau
    if len(cand) > 1:
        keep = np.ones(len(cand), bool)
        for k in range(1, len(cand)):
            if cand[k] - cand[k - 1] == 1 and y[cand[k]] == y[cand[k - 1]]:
                keep[k] = False
        cand = cand[keep]
    if distance <= 1 or len(cand) < 2:
        return cand
    order = cand[np.argsort(-y[cand], kind="stable")]
    taken = np.zeros(len(y), bool)
    out = []
    for i in order:
        if not taken[i]:
            out.append(i)
            taken[max(0, i - distance + 1):i + distance] = True
    return np.array(sorted(out), int)


def isotope_peaks(mass, threshold=0.01, max_n=400):
    """Apex masses (centroid of the top half) and heights of the isotope
    peaks of a resolved mass spectrum, for labels."""
    x, y = mass[:, 0], mass[:, 1]
    if len(x) < 5 or y.max() <= 0:
        return []
    E = _lib31()
    if E is not None:
        r = E.isotope_peaks(mass, threshold, max_n)
        if r is not None:
            return r
    dx = float(np.median(np.diff(x)))
    idx = _find_peaks(y, threshold * y.max(), max(1, int(0.4 / max(dx, 1e-9))))
    out = []
    for i in idx:
        lo, hi = i, i
        while lo > 0 and y[lo - 1] >= 0.5 * y[i] and y[lo - 1] <= y[lo]:
            lo -= 1
        while hi < len(y) - 1 and y[hi + 1] >= 0.5 * y[i] and y[hi + 1] <= y[hi]:
            hi += 1
        seg = slice(lo, hi + 1)
        cen = float(np.sum(x[seg] * y[seg]) / np.sum(y[seg])) if hi > lo else float(x[i])
        out.append((cen, float(y[i])))
    out.sort(key=lambda t: -t[1])
    return sorted(out[:max_n])


def label_maxima_index(x, y, n=6, min_sep=None, rel=0.05):
    """Indices of the n highest local maxima of y at least rel x the tallest of them, further
    apart than min_sep (default: a fortieth of the x span), tallest first; the selection of
    unilcms.label_maxima (on the points in view)."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if len(x) < 3:
        return []
    E = _lib31()
    if E is not None:
        r = E.label_maxima(x, y, n, min_sep, rel)
        if r is not None:
            return r
    lm = np.where((y[1:-1] >= y[:-2]) & (y[1:-1] > y[2:]))[0] + 1
    if not len(lm):
        return []
    order = lm[np.argsort(y[lm], kind="stable")[::-1]]
    top = y[order[0]]
    sep = min_sep if min_sep is not None else (x[-1] - x[0]) / 40.0
    chosen = []
    for i in order:
        if y[i] < rel * top or len(chosen) >= n:
            break
        if all(abs(x[i] - x[j]) > sep for j in chosen):
            chosen.append(int(i))
    return chosen


def group_isotopes(mass, threshold=0.05, spacing=1.00235, ranges=None):
    """Isotope resolved mass spectrum -> one entry per species: isotope peaks
    about 1 Da apart are grouped; mass = average mass (intensity weighted
    over the group), apex = most abundant isotope peak. ranges (a list, if
    given) receives the mass of the first and last isotope peak of each
    species returned."""
    E = _lib31()
    res = E.group_isotopes(mass, threshold, spacing) if E is not None else None
    if res is not None:
        out = [{"mass": a, "apex": ap, "height": h, "area": ar, "n_iso": ni, "z": "most abundant %.3f" % ap}
               for a, ap, h, ar, ni, _f, _l in res]
        if ranges is not None:
            ranges.extend((f0, f1) for _a, _ap, _h, _ar, _ni, f0, f1 in res)
        if not out:
            return []
        tallest = max(q["height"] for q in out)
        tot = sum(q["height"] for q in out) or 1.0
        for q in out:
            q["rel"] = 100.0 * q["height"] / tallest
            q["frac"] = 100.0 * q["height"] / tot
        return out
    iso = pick_mass_peaks(mass, 0.35, threshold * 0.2)
    if not iso:
        return []
    iso.sort(key=lambda q: q["mass"])

    def chains(pk):
        groups, cur = [], [pk[0]]
        for q in pk[1:]:
            if abs(q["mass"] - cur[-1]["mass"] - spacing) <= 0.12 or q["mass"] - cur[-1]["mass"] < 0.6:
                cur.append(q)
            else:
                groups.append(cur)
                cur = [q]
        groups.append(cur)
        return groups

    def dips(groups):
        # an isotope distribution has one maximum: a deep dip in the heights
        # inside a chain means two species next to each other (e.g. +16 Da)
        split = []
        for g in groups:
            while len(g) >= 5:
                h = np.array([q["height"] for q in g])
                cut = None
                for j in range(1, len(g) - 1):
                    # a local minimum with at least twice its height on both sides
                    if h[j] <= h[j - 1] and h[j] <= h[j + 1] and h[j] <= 0.5 * min(h[:j].max(), h[j + 1:].max()):
                        cut = j
                        break
                if cut is None:
                    break
                split.append(g[:cut])
                g = g[cut + 1:] if h[cut + 1] > h[cut] else g[cut:]
            if g:
                split.append(g)
        return split

    # each species is limited to its isotope envelope: the most abundant
    # isotope +- 2.5 averagine sigma (+- 4 Da at 12 kDa); what lies outside is
    # grouped again (another species when it has its own maximum, dropped as
    # the tail of this one when it only falls away from it). 2.7 let a weak
    # species swallow a run of small neighbours (12273 Da in the 12 kDa test
    # protein: 55 isotope peaks, average mass 4.8 Da above its envelope)
    groups = []

    def limit(g, side):
        hs = [q["height"] for q in g]
        t = int(np.argmax(hs))
        if (side < 0 and t == len(g) - 1) or (side > 0 and t == 0):
            return
        sig = averagine_sigma(g[t]["mass"])
        lo, hi = g[t]["mass"] - 2.5 * sig - 0.5, g[t]["mass"] + 2.5 * sig + 0.5
        groups.append([q for q in g if lo <= q["mass"] <= hi])
        left = [q for q in g if q["mass"] < lo]
        right = [q for q in g if q["mass"] > hi]
        if left:
            for part in dips(chains(left)):
                limit(part, -1)
        if right:
            for part in dips(chains(right)):
                limit(part, 1)
    for g in dips(chains(iso)):
        limit(g, 0)
    groups.sort(key=lambda g: g[0]["mass"])
    x, y = mass[:, 0], mass[:, 1]
    out = []
    for g in groups:
        lo, hi = g[0]["mass"] - 0.5, g[-1]["mass"] + 0.5
        m = (x >= lo) & (x <= hi)
        avg = float(np.sum(x[m] * y[m]) / np.sum(y[m])) if np.any(m) and np.sum(y[m]) > 0 else g[0]["mass"]
        top = max(g, key=lambda q: q["height"])
        dx = float(np.median(np.diff(x))) if len(x) > 1 else 1.0
        out.append({"mass": avg, "apex": top["mass"], "height": top["height"],
                    "area": float(np.sum(y[m]) * dx) if np.any(m) else 0.0, "n_iso": len(g),
                    "z": "most abundant %.3f" % top["mass"]})
    tallest = max(q["height"] for q in out)
    keep = [k for k, q in enumerate(out) if q["height"] >= threshold * tallest]
    if ranges is not None:
        ranges.extend((groups[k][0]["mass"], groups[k][-1]["mass"]) for k in keep)
    out = [out[k] for k in keep]
    tot = sum(q["height"] for q in out) or 1.0
    for q in out:
        q["rel"] = 100.0 * q["height"] / tallest
        q["frac"] = 100.0 * q["height"] / tot
    return out


def estimate_peak_width(spec):
    """FWHM (m/z) of the tallest peak, half height crossings interpolated."""
    x, y = np.asarray(spec[:, 0], float), np.asarray(spec[:, 1], float)
    i = int(np.argmax(y))
    w = _fwhm_at(x, y, i)
    if w is None or not np.isfinite(w) or w <= 0:
        w = 5 * float(np.median(np.diff(x)))
    return float(w)


# ============================================================================
# IsoDec
# ============================================================================
def isodec_deconvolute(spec, p, log=None):
    from unidec.IsoDec.runtime import IsoDecRuntime
    spec = _restrict(spec, p.get("mz_range"))
    spec = spec[spec[:, 1] > 0]
    if len(spec) < 5:
        raise ValueError("Too few data points in the m/z range")
    rt = IsoDecRuntime(phaseres=int(p.get("phaseres", 8) or 8))
    try:
        rt.config.adductmass = carrier_mass(p)
    except Exception:
        pass
    for key in ("css_thresh", "matchtol", "minpeaks", "maxshift", "knockdown_rounds", "min_score_diff",
                "minareacovered", "isotopethreshold", "datathreshold", "zscore_threshold",
                "background_subtraction"):
        v = p.get("iso_" + key)
        if v is not None and v != "" and hasattr(rt.config, key):
            cur = getattr(rt.config, key)
            setattr(rt.config, key, type(cur)(v) if isinstance(cur, (int, float)) else v)
    zlo, zhi = int(p["z_range"][0]), int(p["z_range"][1])
    pks = rt.batch_process_spectrum(spec, window=int(p.get("iso_window", 5)),
                                    threshold=float(p.get("iso_thresh", 0.0001)), refresh=True)
    lo_m, hi_m = float(p["mass_range"][0]), float(p["mass_range"][1])
    peaks = []
    masses = list(getattr(pks, "masses", []) or [])
    if masses:
        for mm in masses:
            zs = sorted(set(int(z) for z in np.atleast_1d(getattr(mm, "zs", [])) if zlo <= int(z) <= zhi))
            if not zs or not (lo_m <= mm.monoiso <= hi_m):
                continue
            inten = float(getattr(mm, "totalintensity", 0) or getattr(mm, "apexintensity", 0) or 0)
            peaks.append({"mass": float(mm.monoiso), "height": inten, "area": inten,
                          "z": ", ".join(str(z) for z in zs),
                          "mz": float(np.atleast_1d(getattr(mm, "mzs", [0]))[0]),
                          "avg": float(getattr(mm, "avgmass", 0) or 0)})
    else:
        for mp in getattr(pks, "peaks", []):
            if not (zlo <= int(mp.z) <= zhi) or not (lo_m <= mp.monoiso <= hi_m):
                continue
            inten = float(getattr(mp, "matchedintensity", 0) or getattr(mp, "peakint", 0) or 0)
            peaks.append({"mass": float(mp.monoiso), "height": inten, "area": inten, "z": str(int(mp.z)),
                          "mz": float(mp.mz), "avg": float(getattr(mp, "avgmass", 0) or 0)})
    peaks.sort(key=lambda q: q["mass"])
    top = max([q["height"] for q in peaks] + [1e-12])
    peaks = [q for q in peaks if q["height"] >= float(p.get("peak_thresh", 0) or 0) * top]
    tot = sum(q["height"] for q in peaks) or 1.0
    for q in peaks:
        q["rel"] = 100.0 * q["height"] / top
        q["frac"] = 100.0 * q["height"] / tot
    step = float(p.get("mass_step", 0.01) or 0.01)
    step = min(step, 0.05)
    if peaks:
        # drawn as narrow Gaussians, each evaluated only around its mass; the
        # grid is thinned when the masses spread over a wide range (a 100 to
        # 50000 Da range at 0.01 Da would be 5 million points)
        span = peaks[-1]["mass"] - peaks[0]["mass"] + 15
        step = max(step, span / 400000.0)
        grid = np.arange(peaks[0]["mass"] - 5, peaks[-1]["mass"] + 10, step)
        y = np.zeros(len(grid))
        sig = max(3 * step, 0.02)
        for q in peaks:
            a, b = np.searchsorted(grid, [q["mass"] - 6 * sig, q["mass"] + 6 * sig])
            if b > a:
                y[a:b] = np.maximum(y[a:b], q["height"] * np.exp(-0.5 * ((grid[a:b] - q["mass"]) / sig) ** 2))
        mass = np.column_stack([grid, y])
    else:
        mass = np.zeros((0, 2))
    return {"method": "IsoDec", "mass": mass, "fit": None, "data": spec, "peaks": peaks, "r2": None,
            "zdist": None, "notes": "IsoDec: %d monoisotopic mass%s from isotope distributions, charge %d to %d" % (
                len(peaks), "es" if len(peaks) != 1 else "", zlo, zhi)}


# ============================================================================
# one entry point for the window and the worker process
# ============================================================================
def min_intensity(spec, p):
    """Minimum intensity in data units (counts, as on the y axis) from
    p["min_int"] and p["min_int_unit"] ("abs" or "pct" of the base peak, the
    tallest point in the m/z range used). Returns (threshold, base peak)."""
    y = np.asarray(spec, float)[:, 1] if len(spec) else np.zeros(0)
    base = float(y.max()) if len(y) else 0.0
    v = float(p.get("min_int", 0) or 0)
    if v <= 0:
        return 0.0, base
    return (v / 100.0 * base if p.get("min_int_unit") == "pct" else v), base


def peak_mask(spec, thr):
    """Points that belong to a peak whose top reaches thr: every such peak is
    kept whole, down to the valleys on both sides (searched within about
    five peak widths), so peak shapes and widths are not cut. Returns
    (mask, number of peaks)."""
    x, y = np.asarray(spec[:, 0], float), np.asarray(spec[:, 1], float)
    keep = np.zeros(len(y), bool)
    if len(y) < 3:
        return keep | (y >= thr), int((y >= thr).sum())
    E = _lib31()
    if E is not None:
        r = E.peak_mask(spec, thr)
        if r is not None:
            return r
    from scipy.signal import find_peaks, peak_prominences
    idx, _ = find_peaks(y, height=thr, plateau_size=(1, None))
    if not len(idx):
        return keep, 0
    fw = estimate_peak_width(spec)
    i0 = int(np.argmax(y))
    lo, hi = max(0, i0 - 50), min(len(x) - 1, i0 + 50)
    dx = float(np.median(np.diff(x[lo:hi + 1]))) if hi > lo else float(np.median(np.diff(x)))
    n_fw = fw / dx if dx > 0 else 5.0
    wlen = int(max(7, min(10 * n_fw, len(y)))) | 1
    _, left, right = peak_prominences(y, idx, wlen=wlen)
    for a, b in zip(left, right):
        keep[int(a):int(b) + 1] = True
    return keep, int(len(idx))


def data_summary(spec, p, with_mask=False):
    """What the deconvolution uses: m/z range, base peak, minimum intensity
    and the peaks and points kept (for the notes and the report)."""
    spec = _restrict(spec, p.get("mz_range"))
    thr, base = min_intensity(spec, p)
    y = spec[:, 1] if len(spec) else np.zeros(0)
    if thr > 0 and len(y):
        mask, npk = peak_mask(spec, thr)
    else:
        mask = y > 0
        npk = int(((y[1:-1] > y[:-2]) & (y[1:-1] >= y[2:]) & (y[1:-1] > 0)).sum()) if len(y) > 2 else 0
    out = {"mz": (float(spec[0, 0]), float(spec[-1, 0])) if len(spec) else None, "base": base,
           "base_mz": float(spec[int(np.argmax(y)), 0]) if len(y) else None, "threshold": thr,
           "unit": p.get("min_int_unit", "abs"), "value": float(p.get("min_int", 0) or 0),
           "points": len(y), "points_used": int(mask.sum()), "peaks_used": npk}
    if with_mask:
        out["_mask"] = mask
    return out


def _fmt_int(v):
    """Intensity for the notes: 850, 42000, 1.25e6."""
    a = abs(float(v))
    if a >= 1e5:
        return ("%.2e" % v).replace("e+0", "e").replace("e+", "e")
    return "%.0f" % v if a >= 100 else "%.3g" % v


def summary_text(s):
    """One line for the notes, e.g. m/z 400 to 2000, minimum intensity 5.00e3
    (2 % of the base peak 2.50e5)."""
    if not s or not s.get("mz"):
        return ""
    t = "m/z %.2f to %.2f" % s["mz"]
    if s["threshold"] > 0:
        pct = 100.0 * s["threshold"] / s["base"] if s["base"] > 0 else 0.0
        t += (", peaks from %s up (%.3g %% of the base peak %s): %d peaks kept whole, %d of %d points" % (
            _fmt_int(s["threshold"]), pct, _fmt_int(s["base"]), s["peaks_used"], s["points_used"], s["points"]))
    else:
        t += ", no minimum intensity (base peak %s)" % _fmt_int(s["base"])
    return t


def _engine_method(method, spec, folder, name, p, progress):
    """UniDec or IsoDec through the C++ core; None when the core is not
    available or failed (the Python code is used then). A ValueError (bad
    settings) is raised as from the Python code."""
    try:
        import msengine_py
        if method == "unidec":
            return msengine_py.unidec(spec, folder, name, p, progress=progress)
        return msengine_py.isodec(spec, p)
    except ValueError:
        raise
    except Exception as ex:
        if str(ex).startswith("The UniDec engine"):
            # the engine itself failed: the Python code would run the same
            # engine again on the same input (twice the wait, same failure)
            raise
        import sys
        mp = sys.modules.get("msengine_py")
        if method == "unidec" and mp is not None and getattr(mp, "unidec_in_library", lambda: False)():
            # the library runs UniDec's engine itself (3.1): its failure is the
            # result; the Python code would need unidec.exe, which MS Analysis
            # no longer uses
            raise
        print("[engine] Python %s used: %s" % (method, ex))
        return None


def run_method(method, spec, folder, name, p, progress=None):
    """Minimum intensity and baseline subtraction (if asked), then the
    method; progress(i, n, text). Returns the result dict of the method.

    The minimum intensity is in the units of the spectrum (or % of its base
    peak in the m/z range) and is judged on the data as shown: peaks whose
    top is below it are removed (set to zero, after any baseline
    subtraction); the peaks above it are kept whole, for every method."""
    if method == "maxent":
        # All processing/calculations happen in the native pipeline. Never silently
        # switch to a different Python solver when a native component is missing.
        import msengine_py
        callback = None if progress is None else lambda i, n, chi: progress(
            i, n, "maximum entropy (C++): round %d of up to %d" % (i, n))
        return msengine_py.maxent(spec, dict(p), progress=callback)
    spec = _restrict(spec, p.get("mz_range"))
    p = dict(p)
    summ = data_summary(spec, p, with_mask=True)
    keep = summ.pop("_mask")
    thr = summ["threshold"]
    if method in ("unidec", "isodec"):
        # the C++ core when available: it does the whole of this function
        # (minimum intensity, baseline, binning, the engine, peak picking)
        res = _engine_method(method, spec, folder, name, p, progress)
        if res is not None:
            inp = subtract_baseline(spec, factor=float(p.get("baseline_width", 15.0))) if p.get("baseline") else spec
            if thr > 0:
                inp = np.array(inp, float)
                inp[~keep, 1] = 0.0
            res["input_data"] = inp
            res["data_info"] = summ
            refine_species(res.get("mass"), res.get("peaks"), data=inp, p=p)
            return res
    if p.pop("baseline", False):
        spec = subtract_baseline(spec, factor=float(p.get("baseline_width", 15.0)))
    if thr > 0 and method == "unidec" and len(spec) > 10:
        # UniDec with gaps of zeros in the data: (1) its automatic peak width
        # is measured on the data before the weak peaks are removed (on the
        # gaps it came out far too wide or narrow); (2) the data are put on
        # bins of their own spacing: its engine stopped (0xC0000409) on data
        # with zero gaps when they were not binned (uniform LC-MS data keep
        # their exact grid)
        if not float(p.get("peak_width", 0) or 0) > 0:
            try:
                from unidec import tools as ud
                original_autocorr = ud.autocorr
                ud.autocorr = _autocorr_from_lag2  # (as the C++ library)
                try:
                    fwhm, psfun, _mid = ud.auto_peak_width(np.array(spec, float))
                finally:
                    ud.autocorr = original_autocorr
                if np.isfinite(fwhm) and 0 < fwhm < 0.2 * (spec[-1, 0] - spec[0, 0]):
                    p["peak_width"] = float(fwhm)
                    p["auto_width_measured"] = True
                    if psfun is not None:
                        p["psfun"] = int(psfun)
            except Exception:
                pass
        if p.get("binning", "auto") in ("auto", "none") and len(spec) <= 60000:
            p["binning"], p["mzbins"] = "linear", float(np.median(np.diff(spec[:, 0])))
            p["native_bins"] = True  # bins = the data spacing: the peak width is not widened
    if thr > 0:
        spec = np.array(spec, float)
        spec[~keep, 1] = 0.0
        if summ["peaks_used"] < 1 or summ["points_used"] < 5 or not (spec[:, 1] > 0).any():
            raise ValueError("No data above the minimum intensity (%s; the base peak in the m/z range is %s). "
                             "Lower the minimum intensity." % (_fmt_int(thr), _fmt_int(summ["base"])))
    if method == "unidec":
        res = unidec_deconvolute(spec, folder, name, p, progress=progress)
    elif method == "isodec":
        if progress:
            progress(1, 2, "IsoDec: isotope distributions")
        res = isodec_deconvolute(spec, p)
    else:
        raise ValueError("unknown method %r" % method)
    res["input_data"] = spec
    res["data_info"] = summ
    extra = summary_text(summ)
    if extra and extra not in (res.get("notes") or ""):  # the C++ core writes it itself
        res["notes"] = (res.get("notes") or "") + "; " + extra
    refine_species(res.get("mass"), res.get("peaks"), data=spec, p=p)
    return res

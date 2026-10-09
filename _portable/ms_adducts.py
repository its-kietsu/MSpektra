# -*- coding: utf-8 -*-
"""Suggested neutral mass of a compound from its averaged ESI+ and ESI- spectra of a unit resolution
quadrupole LC-MS (Shimadzu LCMS-2020 style: profile spectra binned at 0.05 to 0.1 m/z, masses good to
about 0.2 to 0.3). No wx here. The calculation runs in the C++ library (msengine, src/adducts.cpp,
ms_neutral_masses2) when it is there; the code below is the fallback and the reference (same steps, same
arithmetic in the same order: tests/check_adducts.py compares the two bit for bit).

    neutral_masses(pos, neg, tol=0.3, min_rel=0.05, max_ions=12, top=3)
    neutral_masses_info(pos, neg, tol=0.3, min_rel=0.05, max_ions=12, top=3)

pos, neg: (m/z array, intensity array) or an N x 2 array of the positive and the negative spectrum, or
None. neutral_masses returns a list (best first, at most top) of dicts {"mass", "score", "both" (ions of
both polarities support it), "assumed" (a single ion taken as [M+H]+ or [M-H]-), "additive" (None, or the
mobile phase additive the mass matches: "TFA", "formic acid", "acetic acid", "DFA"), "m2" (None, or the
mass of its M+2 partner: one Br or Cl), "m2_ratio" (None, or the height of the partner's strongest ion /
the height of its counterpart ion here: about 1 for one Br, about 0.3 for one Cl), "ions": [{"polarity":
"+" or "-", "adduct": "[M+H]+", "mz", "rel" (height relative to the base peak of its spectrum), "mass"
(the neutral mass this ion gives)}], "m2_ions": the ions of the M+2 partner (same form)}; [] when there is
nothing. neutral_masses_info returns {"masses": that list, "base": {"+": apex m/z of the strongest peak of
the positive spectrum or None, "-": the same of the negative one}, "explained": {"+": True when the best
candidate explains that peak (one of its ions or of its partner's ions), "-": ...}}. Bad arguments raise
ValueError.

The algorithm:
 1. Points: those with a NaN or infinite m/z or intensity are left out, negative intensities count as
    zero, the rest is sorted by m/z (stable sort). Base peak: the highest intensity; none above zero: no
    peaks.
 2. Maxima: runs of equal intensity whose neighbours on both sides are lower or missing; a neighbour
    counts only when it lies within tol in m/z (a wider gap breaks the profile, so that centroid lists
    and spectra without their zero points work as well). The maximum is the middle point of the run
    ((first + last) // 2), its height the intensity of the run. Kept: height > 0 and height >= min_rel x
    base.
 3. Apex: the intensity weighted mean m/z of the points of the peak, from the run outwards as long as the
    next point is connected (within tol), not higher than the point before it and at least half the
    height (sums from left to right).
 4. Ions: the maxima ranked by height (descending), then position. A maximum is dropped when one of
    better rank lies at d = (its apex - the other's apex) with |d| <= tol (the same peak: a jagged top or
    a shoulder), |d - ISO/2| <= tol (isotope of a doubly charged ion) or |d - ISO| <= tol (isotope;
    ISO = 1.00335483507, 13C - 12C). The first max_ions of the others are the ions of the spectrum,
    rel = height / base.
 5. Candidates: each ion with each adduct of its polarity gives M = (z x mz - d) / k (the ion's m/z is
    (k M + d) / z, d from the exact masses: proton 1.00727646688, Na+, K+, NH4+, Cl-, HCOO-, CH3COO- from
    the atomic masses minus or plus an electron). Only M > 0.
 6. Groups: the candidates sorted by M (then ion: positive ions first, each polarity in rank order; then
    adduct code); a candidate joins the group of the one before it when it is at most f x tol heavier,
    f = the larger z / k of the two adducts (the error of M is z / k x the error of the m/z: 2 for doubly
    charged ions, 1, 0.5 for dimers). In a group each peak keeps only its most common adduct
    (lowest code), then each adduct only its strongest ion (highest rel, then the first ion). A group
    needs two ions or more after that. Its ions in the order: positive first, rel descending, m/z, code.
    Mass = sum(rel x M) / sum(rel) over the ions in that order.
 7. Score = sum(rel x w) in that order, w = 1 for [M+H]+, [M-H]- and [M+Na]+, 0.8 for the other singly
    charged ions of one molecule, 0.5 for doubly charged and dimer ions; then x 1.5 when ions of both
    polarities support the group, x 0.5 when it rests on doubly charged and dimer ions only. Not
    normalized: a well supported compound scores about 1 or more, a coincidence of two weak peaks a few
    hundredths.
 8. Mobile phase additives: a candidate whose mass lies within max(tol, 0.5) of TFA (113.99286), formic
    acid (46.00548), acetic acid (60.02113) or DFA (difluoroacetic acid, 96.00229, the LC-MS alternative
    to TFA) is marked with that additive (the first in this order). At least 0.5: unit resolution
    quadrupoles are often 0.3 to 0.4 off at low m/z (TFA [M+H]+ measured at 115.35).
 9. M+2 partners (one Br or Cl; M2 = 1.997): the groups in the order of mass (stable); for each A not
    removed, a heavier group B not removed is its partner when |mass B - mass A - M2| <= tol and every
    ion of B at least half as high as B's strongest has an ion of A of the same adduct at m/z
    k x M2 / z lower (within tol; the counterparts of weaker ions can be missing: below min_rel, beyond
    the max_ions strongest, or lost in the isotope lump of a doubly charged ion); of several, the best
    scored (then the lighter). A keeps its mass and ions, gets m2 = mass of B, the ions of B, m2_ratio =
    rel of B's strongest ion (the first of the highest) / rel of its counterpart in A, and score A +
    score B (also when B scored higher: the lighter, monoisotopic one stays); B is removed.
10. Reinterpretations: the groups by score (descending; then mass, stable); one whose peaks (its ions
    and its partner's: the picked maxima) are all peaks of better groups kept is removed (it only
    rereads ions already explained: 2M and M/2 readings, [M+K]+ / [M+Cl]- readings of an [M+H]+ /
    [M-H]- pair). The peaks of the additive groups are background: a group other than an additive with
    fewer than two peaks outside them is removed too (it rests on one peak and an ion of the mobile
    phase).
11. No group left other than additives: every ion of no group kept alone, taken as [M+H]+ (positive)
    or [M-H]- (negative), the adducts "assumed [M+H]+" / "assumed [M-H]-", w = 0.1 (score 0.1 x rel;
    marked by step 8 too); step 9 among these (a single ion and its M+2 partner).
12. Additives: score x 0.1, and they come after all other candidates.
13. Sorted: additives last, then by score (descending), then mass (stable); the first top are returned.
14. Base peaks: per polarity the apex m/z of the strongest picked peak (always the first ion: nothing
    ranks above it) and whether the best candidate explains it (it is one of its ions or of its
    partner's ions).
"""
import bisect

import numpy as np

PROTON = 1.00727646688
ISO = 1.00335483507
M2 = 1.997   # 81Br - 79Br 1.99795, 37Cl - 35Cl 1.99705

# code: (name, polarity, k, z, d, weight); m/z = (k M + d) / z; a lower code is the more common adduct
ADDUCTS = [
    ("[M+H]+", 1, 1, 1, PROTON, 1.0),
    ("[M-H]-", -1, 1, 1, -PROTON, 1.0),
    ("[M+Na]+", 1, 1, 1, 22.989220702091, 1.0),
    ("[M+NH4]+", 1, 1, 1, 18.033825553441, 0.8),
    ("[M+HCOO]-", -1, 1, 1, 44.998202851279, 0.8),
    ("[M+K]+", 1, 1, 1, 38.963157906491, 0.8),
    ("[M+Cl]-", -1, 1, 1, 34.969401261909, 0.8),
    ("[M+CH3COO]-", -1, 1, 1, 59.013852915739, 0.8),
    ("[M+2H]2+", 1, 1, 2, 2.0 * PROTON, 0.5),
    ("[M-2H]2-", -1, 1, 2, -2.0 * PROTON, 0.5),
    ("[2M+H]+", 1, 2, 1, PROTON, 0.5),
    ("[2M+Na]+", 1, 2, 1, 22.989220702091, 0.5),
    ("[2M-H]-", -1, 2, 1, -PROTON, 0.5),
    ("assumed [M+H]+", 1, 1, 1, PROTON, 0.1),
    ("assumed [M-H]-", -1, 1, 1, -PROTON, 0.1),
]
ADDUCT_NAMES = [a[0] for a in ADDUCTS]
N_REGULAR = 13
ASSUMED_POS, ASSUMED_NEG = 13, 14
# mobile phase additives (neutral monoisotopic masses): CF3COOH, HCOOH, CH3COOH, CHF2COOH
ADDITIVES = [("TFA", 113.99286375956), ("formic acid", 46.0054793036), ("acetic acid", 60.02112936806),
             ("DFA", 96.00228562906)]


def _lib():
    """msengine_py when the C++ library is in use, else None (the Python code runs)."""
    try:
        import msengine_py
        return msengine_py if msengine_py.available() else None
    except Exception:
        return None


def check_args(tol, min_rel, max_ions, top):
    """The arguments as numbers; ValueError for values the calculation does not accept (the library
    returns MS_BAD_ARG for the same ones)."""
    tol, min_rel = float(tol), float(min_rel)
    if not (0 < tol <= 5):
        raise ValueError("the tolerance must be above 0 and at most 5")
    if not (0 <= min_rel <= 1):
        raise ValueError("min_rel must be 0 to 1")
    max_ions, top = int(max_ions), int(top)
    if not (1 <= max_ions <= 1000):
        raise ValueError("max_ions must be 1 to 1000")
    if not (1 <= top <= 1000):
        raise ValueError("top must be 1 to 1000")
    return tol, min_rel, max_ions, top


def spectrum_xy(s):
    """(m/z, intensity) float64 arrays of a spectrum given as (mz, it), an N x 2 array or None."""
    if s is None:
        return np.zeros(0), np.zeros(0)
    if isinstance(s, np.ndarray) and s.ndim == 2:
        if s.shape[1] == 2:
            x, y = s[:, 0], s[:, 1]
        elif s.shape[0] == 2:
            x, y = s[0], s[1]
        else:
            raise ValueError("a spectrum is (mz, intensity) or an N x 2 array")
    else:
        if len(s) != 2:
            raise ValueError("a spectrum is (mz, intensity) or an N x 2 array")
        x, y = s
    x = np.ascontiguousarray(x, dtype=np.float64).ravel()
    y = np.ascontiguousarray(y, dtype=np.float64).ravel()
    if len(x) != len(y):
        raise ValueError("m/z and intensity of a spectrum differ in length")
    return x, y


def _spectrum_ions(x, y, pol, tol, min_rel, max_ions):
    """Steps 1 to 4: the ions (polarity, apex m/z, rel) of one spectrum, in rank order."""
    ok = np.isfinite(x) & np.isfinite(y)
    x = x[ok]
    y = np.where(y[ok] > 0, y[ok], 0.0)
    m = len(x)
    if m == 0:
        return []
    o = np.argsort(x, kind="stable")
    x, y = x[o], y[o]
    base = float(np.max(y))
    if not base > 0:
        return []
    thr = min_rel * base
    # 2. maxima (runs of equal connected points)
    conn = (x[1:] - x[:-1]) <= tol
    same = conn & (y[1:] == y[:-1])
    starts = np.flatnonzero(np.concatenate(([True], ~same)))
    ends = np.flatnonzero(np.concatenate((~same, [True])))
    h = y[starts]
    left = ~np.concatenate(([False], conn))[starts] | (y[np.maximum(starts - 1, 0)] < h)
    right = ~np.concatenate((conn, [False]))[ends] | (y[np.minimum(ends + 1, m - 1)] < h)
    keep = np.flatnonzero(left & right & (h > 0) & (h >= thr))
    xs, ys = x.tolist(), y.tolist()
    mx = []   # (pos, height, apex)
    for q in keep.tolist():
        a, b = int(starts[q]), int(ends[q])
        hq = ys[a]
        # 3. apex
        lo, hi = a, b
        half = 0.5 * hq
        while lo > 0 and xs[lo] - xs[lo - 1] <= tol and ys[lo - 1] <= ys[lo] and ys[lo - 1] >= half:
            lo -= 1
        while hi < m - 1 and xs[hi + 1] - xs[hi] <= tol and ys[hi + 1] <= ys[hi] and ys[hi + 1] >= half:
            hi += 1
        sx = sw = 0.0
        for i in range(lo, hi + 1):
            sx += xs[i] * ys[i]
            sw += ys[i]
        mx.append(((a + b) // 2, hq, sx / sw))
    # 4. rank, isotopes and shoulders
    k = len(mx)
    order = sorted(range(k), key=lambda p: (-mx[p][1], mx[p][0]))
    rank = [0] * k
    for r, p in enumerate(order):
        rank[p] = r
    by_apex = sorted(range(k), key=lambda p: (mx[p][2], p))
    apx = [mx[p][2] for p in by_apex]
    out = []
    for r in range(k):
        if len(out) >= max_ions:
            break
        pos, hq, aq = mx[order[r]]
        lo = bisect.bisect_left(apx, aq - (ISO + 2.0 * tol))
        hi = bisect.bisect_right(apx, aq + 2.0 * tol)
        drop = False
        for j in range(lo, hi):
            p = by_apex[j]
            if rank[p] >= r:
                continue
            d = aq - mx[p][2]
            if abs(d) <= tol or abs(d - 0.5 * ISO) <= tol or abs(d - ISO) <= tol:
                drop = True
                break
        if not drop:
            out.append((pol, aq, hq / base))
    return out


def _finish(members, ions, assumed, tol):
    """A candidate (steps 6 to 8): its ions in their order, mass, score, both, additive."""
    members = sorted(members, key=lambda c: (-ions[c[1]][0], -ions[c[1]][2], ions[c[1]][1], c[2]))
    sm = sr = s = 0.0
    pos = neg = single = False
    for M, i, code in members:
        pol, mz, rel = ions[i]
        sm += rel * M
        sr += rel
        s += rel * ADDUCTS[code][5]
        if pol > 0:
            pos = True
        else:
            neg = True
        if ADDUCTS[code][2] == 1 and ADDUCTS[code][3] == 1:
            single = True
    both = pos and neg
    if both:
        s = s * 1.5
    if not single:
        s = s * 0.5
    mass = sm / sr
    additive = -1
    for a, (name, am) in enumerate(ADDITIVES):
        if abs(mass - am) <= max(tol, 0.5):
            additive = a
            break
    return {"m": members, "m2m": [], "mass": mass, "score": s, "m2": None, "m2_ratio": None, "both": both,
            "assumed": assumed, "additive": additive, "removed": False}


def _m2_partner(A, B, ions, tol):
    """Step 9: the index in A's ions of the counterpart of each ion of B (None: an ion weaker than half
    of B's strongest without one) when B is A's M+2 partner, else None."""
    if not (B["mass"] > A["mass"]) or not (abs(B["mass"] - A["mass"] - M2) <= tol):
        return None
    top = max(ions[ib][2] for Mb, ib, code in B["m"])
    cp = []
    for Mb, ib, code in B["m"]:
        name, pol, k, z, d, w = ADDUCTS[code]
        found = None
        for j, (Ma, ia, ca) in enumerate(A["m"]):
            if ca == code and abs(ions[ib][1] - ions[ia][1] - k * M2 / z) <= tol:
                found = j
                break
        if found is None and not ions[ib][2] < 0.5 * top:
            return None
        cp.append(found)
    return cp


def _merge_m2(G, idx, ions, tol):
    """Step 9 among the candidates G[i], i in idx: in the order of mass, B joins A."""
    by_mass = sorted(idx, key=lambda q: G[q]["mass"])
    for ia in by_mass:
        A = G[ia]
        if A["removed"]:
            continue
        best = best_cp = None
        for ib in by_mass:
            if ib == ia or G[ib]["removed"]:
                continue
            cp = _m2_partner(A, G[ib], ions, tol)
            if cp is None:
                continue
            if best is None or G[ib]["score"] > G[best]["score"] or (G[ib]["score"] == G[best]["score"]
                                                                      and G[ib]["mass"] < G[best]["mass"]):
                best, best_cp = ib, cp
        if best is None:
            continue
        B = G[best]
        k = 0   # B's strongest ion (the first of the highest)
        for i in range(1, len(B["m"])):
            if ions[B["m"][i][1]][2] > ions[B["m"][k][1]][2]:
                k = i
        A["m2"] = B["mass"]
        A["m2_ratio"] = ions[B["m"][k][1]][2] / ions[A["m"][best_cp[k]][1]][2]
        A["m2m"] = list(B["m"])
        A["score"] = A["score"] + B["score"]
        B["removed"] = True


def _ion_dicts(members, ions):
    return [{"polarity": "+" if ions[i][0] > 0 else "-", "adduct": ADDUCTS[code][0], "mz": ions[i][1],
             "rel": ions[i][2], "mass": M} for M, i, code in members]


def neutral_masses_info_py(pos, neg, tol=0.3, min_rel=0.05, max_ions=12, top=3):
    """The Python reference of neutral_masses_info (the algorithm in the module's description)."""
    tol, min_rel, max_ions, top = check_args(tol, min_rel, max_ions, top)
    px, py = spectrum_xy(pos)
    nx, ny = spectrum_xy(neg)
    ions = _spectrum_ions(px, py, 1, tol, min_rel, max_ions)
    npos = len(ions)
    ions += _spectrum_ions(nx, ny, -1, tol, min_rel, max_ions)
    # 5. candidates (M, ion, code)
    cands = []
    for i, (pol, mz, rel) in enumerate(ions):
        for code in range(N_REGULAR):
            name, apol, k, z, d, w = ADDUCTS[code]
            if apol != pol:
                continue
            M = (z * mz - d) / k
            if M > 0:
                cands.append((M, i, code))
    cands.sort()
    # 6. groups
    G = []
    s = 0
    while s < len(cands):
        e = s + 1
        while e < len(cands):
            P, Q = ADDUCTS[cands[e - 1][2]], ADDUCTS[cands[e][2]]
            if not (cands[e][0] - cands[e - 1][0] <= tol * max(P[3] / P[2], Q[3] / Q[2])):
                break
            e += 1
        g = cands[s:e]
        a = [c for c in g if not any(o[1] == c[1] and o[2] < c[2] for o in g)]
        b = []
        for i, c in enumerate(a):
            rc = ions[c[1]][2]
            if not any(j != i and o[2] == c[2] and (ions[o[1]][2] > rc or (ions[o[1]][2] == rc and o[1] < c[1]))
                       for j, o in enumerate(a)):
                b.append(c)
        if len(b) >= 2:
            G.append(_finish(b, ions, False, tol))
        s = e
    # 9. M+2 partners and 10. reinterpretations among the groups
    ng = len(G)
    _merge_m2(G, range(ng), ions, tol)
    order = sorted([q for q in range(ng) if not G[q]["removed"]], key=lambda q: (-G[q]["score"], G[q]["mass"]))
    background = set(i for g in G if not g["removed"] and g["additive"] >= 0 for M, i, code in g["m"] + g["m2m"])
    covered = set()
    for q in order:
        g = G[q]
        peaks = set(i for M, i, code in g["m"] + g["m2m"])
        if peaks <= covered or (g["additive"] < 0 and len(peaks - background) < 2):
            g["removed"] = True
        else:
            covered |= peaks
    # 11. no group left other than additives: each ion of none alone, assumed [M+H]+ / [M-H]-, and their
    # M+2 partners
    if not any(not g["removed"] and g["additive"] < 0 for g in G):
        for i, (pol, mz, rel) in enumerate(ions):
            if i in covered:
                continue
            code = ASSUMED_POS if pol > 0 else ASSUMED_NEG
            name, apol, k, z, d, w = ADDUCTS[code]
            M = (z * mz - d) / k
            if M > 0:
                G.append(_finish([(M, i, code)], ions, True, tol))
        _merge_m2(G, range(ng, len(G)), ions, tol)
    # 12. additives x 0.1 and after the others; 13. best first
    out = []
    for g in G:
        if not g["removed"]:
            if g["additive"] >= 0:
                g["score"] = g["score"] * 0.1
            out.append(g)
    out.sort(key=lambda g: (g["additive"] >= 0, -g["score"], g["mass"]))
    out = out[:top]
    # 14. base peaks and whether the best candidate explains them
    base, explained = {}, {}
    for p, first in (("+", 0 if npos > 0 else None), ("-", npos if len(ions) > npos else None)):
        base[p] = ions[first][1] if first is not None else None
        explained[p] = bool(first is not None and out and any(i == first for M, i, code in out[0]["m"] + out[0]["m2m"]))
    masses = [{"mass": g["mass"], "score": g["score"], "both": g["both"], "assumed": g["assumed"],
               "additive": ADDITIVES[g["additive"]][0] if g["additive"] >= 0 else None, "m2": g["m2"],
               "m2_ratio": g["m2_ratio"], "ions": _ion_dicts(g["m"], ions), "m2_ions": _ion_dicts(g["m2m"], ions)}
              for g in out]
    return {"masses": masses, "base": base, "explained": explained}


def neutral_masses_py(pos, neg, tol=0.3, min_rel=0.05, max_ions=12, top=3):
    """The Python reference of neutral_masses."""
    return neutral_masses_info_py(pos, neg, tol, min_rel, max_ions, top)["masses"]


def neutral_masses_info(pos, neg, tol=0.3, min_rel=0.05, max_ions=12, top=3):
    """Suggested neutral masses of a compound from its positive and negative spectrum and the base peaks
    (see the module's description): {"masses": the list of neutral_masses, "base": {"+": apex m/z of the
    strongest peak of the positive spectrum or None, "-": ...}, "explained": {"+": True when the best
    candidate explains it, "-": ...}}. The library's ms_neutral_masses2 when it is there, else the Python
    reference."""
    tol, min_rel, max_ions, top = check_args(tol, min_rel, max_ions, top)
    px, py = spectrum_xy(pos)
    nx, ny = spectrum_xy(neg)
    E = _lib()
    if E is not None:
        try:
            r = E.neutral_masses_info(px, py, nx, ny, tol, min_rel, max_ions, top)
        except AttributeError:   # an msengine_py without the function
            r = None
        if r is not None:
            return r
    return neutral_masses_info_py((px, py), (nx, ny), tol, min_rel, max_ions, top)


def neutral_masses(pos, neg, tol=0.3, min_rel=0.05, max_ions=12, top=3):
    """Suggested neutral masses of a compound from its positive and negative spectrum (see the module's
    description): the list of candidate dicts, best first."""
    return neutral_masses_info(pos, neg, tol, min_rel, max_ions, top)["masses"]

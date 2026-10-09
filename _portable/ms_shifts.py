# -*- coding: utf-8 -*-
"""Mass shift finder on the masses of a deconvolution result (MS Analysis 3.35): which masses differ by a
known modification or adduct (Na, oxidation, water, carbamylation, ...), by k tags (a dye or another label
attached to the protein) plus at most two other shifts, which differences match nothing (unknown adducts),
and the degree of conjugation (the share of the species with 0, 1, 2, ... tags). No wx here. The
calculation runs in the C++ library (msengine, src/shifts.cpp, ms_mass_shifts) when it is there; the code
below is the fallback and the reference (same steps, same arithmetic in the same order:
tests/check_shifts.py compares the two).

    find_shifts(mass, height, area, shifts, tags, tol_da, tol_ppm, ref=-1, kmax=4, max_extra=2, flags=0,
                min_rel=0.0)

mass, height, area: the masses of the result (Da) with their heights and areas (area None: zeros); shifts:
the mass differences of the shifts in use (Da, signed: a loss is negative), tags: the tag masses (Da,
positive); tol_da + tol_ppm x 1e-6 x (the heavier mass of a pair) is the tolerance of a pair; ref: index
of the reference ("unmodified") mass, -1 = automatic (the tallest, or the lightest of the tag series: step
11); kmax: at most this many tags per molecule (1 to
20); max_extra: other shifts with tags (0 to 2); flags: ISOTOPES (isotope resolved masses: an explanation
may differ by one isotope, 13C - 12C, in either direction) and COLLAPSE (the list holds the isotope peaks
of each species: a peak 1 or 2 isotopes above another one is its isotope peak); min_rel (0 to 1): masses
below this fraction of the tallest one are left out. Returns a dict:
    "ref": index of the reference used (-1 for an empty list),
    "pairs": [{"from", "to" (indexes: the parent and the child, see step 8),
               "delta" (mass[to] - mass[from]), "value" (the explanation, NaN: unknown), "error" (delta -
               value), "kind" (0 a shift, 1 tags with up to max_extra shifts, 2 unknown, 3 isotope peak),
               "sign" (+1 as written, -1 reversed: the parent carries it), "tag" (-1: none), "k", "s1", "s2"
               (shift indexes, -1: none; s1 <= s2), "iso" (-1, 0, +1), "n_alt" (other explanations of
               another value within the tolerance), "alt" (the next best of them: dict with tag, k, s1, s2,
               sign, iso, error; or None),
               "ref" (a pair with the reference), "link" (the pair drawn by default)}],
    "species": per input mass {"status" (0 reference, 1 explained by its pair with the reference, 2
               explained through another mass, 3 unknown, 4 isotope peak of another mass, 5 below min_rel),
               "parent" (the mass its composition comes through; status 4: the first mass of its chain; -1
               none), "link" (index into pairs, -1 none), "value" (its composition in Da), "error", "comp" ((tag counts),
               (shift counts), isotope offset), "height", "area" (with COLLAPSE: of its isotope peaks too)},
    "conj": per tag {"height": [% for k = 0 .. kmax], "area": [...], "avg_height", "avg_area" (tags per
               molecule; NaN without data), "n": species counted}.
Bad arguments raise ValueError (the library returns MS_BAD_ARG for the same ones).

The algorithm (src/shifts.cpp is the same):
 1. Masses sorted (stable: by mass, then index). Heights and areas: not finite or negative count as 0.
    The tolerance of two masses a, b: tol_da + tol_ppm x 1e-6 x max(a, b).
 2. Isotope peaks (COLLAPSE only): in the order of mass, a mass j is the isotope peak of a lighter mass i
    when |m_j - m_i - n ISO| <= tol (n = 1 or 2, ISO = 1.00335483507); of several, the smallest n, then
    the smallest error, then the nearest i. It belongs to the species of i (the first mass of its chain):
    the species' height is the largest of its chain, its area the sum (in the order of mass). Pairs of
    kind 3 (from the first mass of the chain) record it; isotope peaks take no further part.
 3. Species below min_rel x the height of the tallest species take no further part (status 5). Reference:
    the given one (its species, also below min_rel), else the tallest species (the first in input order of
    equal heights).
 4. Explanations: sign x (one shift), or sign x (k x tag + up to max_extra shifts, a shift may repeat, two
    shifts that cancel (|s1 + s2| <= 0.001 Da, e.g. water and loss of water) are left out), k = 1 ..
    kmax, one tag kind at a time; with ISOTOPES each also + iso x ISO, iso = -1 or +1. Value: base = k x
    tag (0 without a tag), base += s1, base += s2 (s1 <= s2), value = sign x base + iso x ISO. Terms: 1
    for the tags (any k), 1 per shift, 1 per isotope offset.
 5. A difference d (tolerance tol) is explained by the candidates with |d - value| <= tol. The best has
    the lowest cost = terms + |d - value| / (0.5 tol) (each half tolerance of error counts as one more
    term), then the smaller error, then sign +, tag, k, s1, s2, iso (ascending). The alternative ("alt")
    is the next one in this order whose value differs from the best's by more than 1e-6 Da (one of the
    same value is the same explanation: water and the reversed loss of water, a tag + double oxidation and
    a tag + 2 oxidations); n_alt counts those of another value.
 6. Explained pairs: every two species a, b (a lighter or equal) with |m_b - m_a| > tol whose difference
    is explained (it then is in both directions; each direction has its own best explanation: a loss as
    written going down, an addition going up).
 7. Compositions, from the reference outwards (its composition is empty): a species j next to a placed
    species i (an explained pair) can have the composition of i plus the explanation from i to j:
    signed counts of each tag and shift and the isotope offset; its value is the sum counts x masses,
    tags first, then the shifts, each in index order, then iso x ISO; its error (m_j - m_ref) - value
    must be within the tolerance of (ref, j), every tag count within -kmax .. kmax, the shift counts at
    most max_extra + 1 together, the isotope offset -1 to 1; its cost is terms (1 per tag kind used, 1 per
    shift, 1 per isotope) + |error| / (0.5 tol). Of the compositions a species can have, the cheapest is kept (then
    the smaller error, then the parent placed first); repeatedly the species with the cheapest one (then
    the smaller error, then its sorted position) is placed, until none is left: the others are unknown.
 8. Pairs: the explained pairs and the unknown pairs (the reference and an unknown species whose difference
    nothing explains), oriented from the parent: of two placed species the
    one placed first, a placed species before an unknown one, of two unknown ones the one nearer in mass
    to the reference (equal: the lighter). Order: by the sorted position of the lighter, then of the
    heavier; the isotope pairs follow.
 9. Links (what is drawn by default): a placed species links to its nearest consistent parent (placed
    before it, with an explained pair, its composition plus the explanation is the species'
    composition; nearest in mass, then sorted position); an unknown species to the parent of its
    explained pair with the fewest terms (then nearest, smaller error, sorted position), else to the
    reference with its unknown pair.
10. Degree of conjugation, per tag: the species (reference and placed ones, in input order) with k of
    this tag and no other tag, 0 <= k <= kmax, summed by height and by area: % = 100 x sum_k / total,
    tags per molecule = (sum of k x sum_k) / total.
11. Reference for tags (only with an automatic reference, ref < 0, and at least one tag): when the tallest
    species is a conjugate (common at high labelling), the unmodified species lies below it and its tag
    counts are negative. Among the reference and the placed species (status 0 to 2), those with every
    tag count <= 0 and at least one < 0: if there are any, the one with the lowest total tag count, then
    the fewest other terms (|shift counts| + |isotope offset|), then the tallest (species height), then
    input order, becomes the reference and steps 3 to 10 run again with it ("ref" reports it).

Masses of the shifts: monoisotopic from the atomic masses of AME 2020 (Wang et al., Chin. Phys. C 2021,
45, 030003), average from the standard atomic weights of IUPAC 2005 (Wieser, Pure Appl. Chem. 2006, 78,
2051), the values Unimod (unimod.org) and ExPASy use, each summed over the elemental composition of the
shift (tests/check_shifts.py compares them with the Unimod values).
"""
import bisect
import math

import numpy as np

ISO = 1.00335483507  # 13C - 12C (as ms_adducts)
ISOTOPES = 1  # flag: isotope resolved masses (an explanation may differ by one isotope)
COLLAPSE = 2  # flag: the list holds the isotope peaks of each species
MAX_N, MAX_SHIFTS, MAX_TAGS, MAX_K = 2000, 64, 3, 20
SETTINGS_REV = [0]  # revision of the settings "mass_shifts": every save adds 1 (deconv_shifts, method presets)

# monoisotopic atomic masses (u), AME 2020
MONO = {"H": 1.00782503190, "C": 12.0, "N": 14.00307400425, "O": 15.99491461926, "F": 18.99840316207,
        "Na": 22.98976928195, "P": 30.97376199768, "S": 31.97207117354, "Cl": 34.968852694, "K": 38.96370648482}
# standard atomic weights, IUPAC 2005 (the average masses of Unimod and ExPASy)
AVERAGE = {"H": 1.00794, "C": 12.0107, "N": 14.0067, "O": 15.9994, "F": 18.9984032, "Na": 22.98976928,
           "P": 30.973762, "S": 32.065, "Cl": 35.453, "K": 39.0983}
ELEMENT_ORDER = ["H", "C", "N", "O", "F", "Na", "P", "S", "Cl", "K"]

# the default list: name and elemental composition of the shift (Unimod names in brackets)
DEFAULT_SHIFTS = [
    ("Na", {"Na": 1, "H": -1}),                                   # Na adduct, Na - H (Cation:Na)
    ("K", {"K": 1, "H": -1}),                                     # K adduct, K - H (Cation:K)
    ("oxidation", {"O": 1}),                                      # (Oxidation)
    ("double oxidation", {"O": 2}),                               # (Dioxidation)
    ("water", {"H": 2, "O": 1}),                                  # hydrolysis, e.g. of a succinimide
    ("loss of water", {"H": -2, "O": -1}),                        # (Dehydrated)
    ("disulfide (−2 H)", {"H": -2}),                         # (Dehydro x 2)
    ("carbamylation", {"C": 1, "H": 1, "N": 1, "O": 1}),          # urea (Carbamyl)
    ("acetylation", {"C": 2, "H": 2, "O": 1}),                    # (Acetyl)
    ("methylation", {"C": 1, "H": 2}),                            # (Methyl)
    ("formylation", {"C": 1, "O": 1}),                            # (Formyl)
    ("phosphorylation", {"H": 1, "P": 1, "O": 3}),                # (Phospho)
    ("phosphate adduct", {"H": 3, "P": 1, "O": 4}),               # H3PO4
    ("sulfate adduct", {"H": 2, "S": 1, "O": 4}),                 # H2SO4
    ("TFA adduct", {"C": 2, "H": 1, "F": 3, "O": 2}),             # CF3COOH
    ("hexose", {"C": 6, "H": 10, "O": 5}),                        # glycation (Hex)
    ("β-mercaptoethanol", {"C": 2, "H": 4, "O": 1, "S": 1}),  # mixed disulfide (BME)
]


def composition_mass(comp):
    """(monoisotopic, average) mass of an elemental composition {element: count} (summed in the order of
    ELEMENT_ORDER)."""
    mono = avg = 0.0
    for el in ELEMENT_ORDER:
        c = comp.get(el, 0)
        if c:
            mono += c * MONO[el]
            avg += c * AVERAGE[el]
    unknown = set(comp) - set(ELEMENT_ORDER)
    if unknown:
        raise ValueError("unknown element: %s" % ", ".join(sorted(unknown)))
    return mono, avg


def default_shifts():
    """The default list: [{"name", "avg", "mono", "on"}]."""
    out = []
    for name, comp in DEFAULT_SHIFTS:
        mono, avg = composition_mass(comp)
        out.append({"name": name, "avg": avg, "mono": mono, "on": True})
    return out


def default_settings():
    """The settings of the finder (key "mass_shifts" of the program's settings)."""
    return {"shifts": default_shifts(), "tags": [], "kmax": 4, "max_extra": 2, "show_all": False,
            "tol_envelope": None, "tol_resolved": None, "min_rel": 0.01}


def _num(v):
    try:
        v = float(str(v).replace(",", "."))
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def clean_settings(s):
    """Settings read from the file (hand edited, older) made valid: missing keys get their default, wrong
    entries are dropped."""
    d = default_settings()
    if not isinstance(s, dict):
        return d
    if isinstance(s.get("shifts"), list):
        out = []
        for e in s["shifts"]:
            if not isinstance(e, dict):
                continue
            avg, mono = _num(e.get("avg")), _num(e.get("mono"))
            if avg is None and mono is None:
                continue
            avg = mono if avg is None else avg
            mono = avg if mono is None else mono
            if avg == 0 or mono == 0 or abs(avg) > 1e5 or abs(mono) > 1e5:
                continue
            out.append({"name": str(e.get("name") or "shift")[:60], "avg": avg, "mono": mono,
                        "on": bool(e.get("on", True))})
        d["shifts"] = out[:MAX_SHIFTS]
    if isinstance(s.get("tags"), list):
        out = []
        for e in s["tags"]:
            if not isinstance(e, dict):
                continue
            avg, mono = _num(e.get("avg")), _num(e.get("mono"))
            avg = mono if avg is None else avg
            mono = avg if mono is None else mono
            if avg is None or avg <= 0 or mono <= 0 or avg > 1e6 or mono > 1e6:
                continue
            out.append({"name": str(e.get("name") or "tag")[:40], "avg": avg, "mono": mono})
        d["tags"] = out[:MAX_TAGS]
    for key, lo, hi in (("kmax", 1, MAX_K), ("max_extra", 0, 2)):
        v = _num(s.get(key))
        if v is not None:
            d[key] = int(min(hi, max(lo, int(v))))
    d["show_all"] = bool(s.get("show_all", False))
    v = _num(s.get("min_rel"))
    if v is not None and 0.0 <= v <= 1.0:
        d["min_rel"] = v
    for key in ("tol_envelope", "tol_resolved"):
        t = s.get(key)
        if isinstance(t, (list, tuple)) and len(t) == 2 and t[1] in ("Da", "ppm"):
            v = _num(t[0])
            if v is not None and v > 0 and v <= (100.0 if t[1] == "Da" else 1e5):
                d[key] = [v, t[1]]
    return d


# ------------------------------------------------------------------ the masses of a result
def _grouped(res):
    """As deconv_tab._grouped: an isotope resolved result with one entry per species."""
    return any(q.get("n_iso") not in (None, "") for q in (res or {}).get("peaks") or [] if not q.get("added"))


SPACING = 1.00235  # mean spacing of the isotope peaks of a protein (Da), as ms_deconv.refine_species
# averagine (Senko 1995): atoms per 111.1254 Da and the abundances of the isotopes 0, +1, +2, ... (nominal)
AVERAGINE_MASS = 111.1254
AVERAGINE = (("C", 4.9384, (0.9893, 0.0107)), ("H", 7.7583, (0.999885, 0.000115)),
             ("N", 1.3577, (0.99636, 0.00364)), ("O", 1.4773, (0.99757, 0.00038, 0.00205)),
             ("S", 0.0417, (0.9499, 0.0075, 0.0425, 0.0, 0.0001)))
_TOP_CACHE = {}


def averagine_top(mass):
    """Index of the most abundant isotope (0: the monoisotopic one) of an averagine protein of this mass
    (atom counts rounded; the isotope distribution by FFT). 7 at 12 kDa, 31 at 50 kDa."""
    key = int(round(float(mass)))
    if key in _TOP_CACHE:
        return _TOP_CACHE[key]
    u = max(0.0, float(mass)) / AVERAGINE_MASS
    mean = var = 0.0
    counts = []
    for _, per, ab in AVERAGINE:
        n = int(round(per * u))
        p = np.array(ab, float)
        k = np.arange(len(p))
        m1 = float((k * p).sum() / p.sum())
        mean += n * m1
        var += n * (float((k * k * p).sum() / p.sum()) - m1 * m1)
        counts.append((n, p))
    size = 1 << int(math.ceil(math.log2(mean + 12.0 * math.sqrt(var) + 16.0)))
    f = np.ones(size // 2 + 1, complex)
    for n, p in counts:
        if n:
            f = f * np.fft.rfft(p, size) ** n
    top = int(np.argmax(np.fft.irfft(f, size)))
    if len(_TOP_CACHE) < 100000:
        _TOP_CACHE[key] = top
    return top


def _iso_centroid(res, ap):
    """The centroid of the labelled isotope peak at ap (within 0.3 Da), else ap."""
    iso = res.get("isotope_peaks") or []
    if len(iso):
        k = min(range(len(iso)), key=lambda j: abs(iso[j][0] - ap))
        if abs(iso[k][0] - ap) <= 0.3:
            return float(iso[k][0])
    return float(ap)


def _apex_mass(res, q):
    """As deconv_tab._apex_mass: the isotope of a species shown in the mass table, the plot label and
    _peaks.csv (a clicked one, else the deconvolution's most abundant one), at its isotope peak."""
    return _iso_centroid(res, float(q.get("selected_apex", q.get("apex_dec", q.get("apex", q["mass"])))))


def _mono_estimate(res, q):
    """Monoisotopic mass of an isotope resolved species estimated from its most abundant isotope (as the
    deconvolution gives it, a clicked isotope does not change it): that isotope minus the averagine
    offset of its most abundant isotope x SPACING. A mass added by hand counts as a most abundant isotope."""
    top = _iso_centroid(res, float(q.get("apex_dec", q.get("apex", q["mass"]))))
    return top - averagine_top(top) * SPACING


def result_basis(res):
    """How the masses of a deconvolution result are compared: {"basis": "average" or "monoisotopic",
    "kind": "envelope", "species" (isotope resolved species), "monoisotopic" (IsoDec) or "isotope peaks"
    (an isotope resolved spectrum whose peaks are its isotope peaks), "flags", "masses" (compared),
    "shown" (as shown in the mass table), "heights", "areas", "real_areas" (False: the method gives no peak
    areas, so the degree of conjugation by area is not shown), "tol" ([value, unit] suited to it)}.
    Isotope envelopes (proteins) are compared by their average masses; isotope resolved results by
    monoisotopic differences: IsoDec gives monoisotopic masses; resolved species are compared by the
    monoisotopic masses estimated from their most abundant isotope (_mono_estimate: a tag of a few kDa moves
    the most abundant isotope by several isotopes; one isotope more or less is allowed, as the averagine
    estimate can be one off) and shown by the isotope of the mass table; the isotope peaks of a resolved
    spectrum are grouped first."""
    peaks = res.get("peaks") or []
    p = res.get("params") or {}
    try:
        step = float(p.get("mass_step", 1.0) or 1.0)
    except (TypeError, ValueError):
        step = 1.0
    heights = [float(q.get("height", 0.0) or 0.0) for q in peaks]
    areas = [float(q.get("area", 0.0) or 0.0) for q in peaks]
    shown = None
    if res.get("tag") == "isodec":
        kind, flags = "monoisotopic", ISOTOPES
        masses = [float(q["mass"]) for q in peaks]
    elif _grouped(res):
        kind, flags = "species", ISOTOPES
        masses = [_mono_estimate(res, q) if q.get("apex") is not None else float(q["mass"]) for q in peaks]
        shown = [_apex_mass(res, q) if "apex" in q else float(q["mass"]) for q in peaks]
    elif res.get("isotope_peaks") and step < 0.1:
        kind, flags = "isotope peaks", ISOTOPES | COLLAPSE
        masses = [float(q["mass"]) for q in peaks]
    else:
        kind, flags = "envelope", 0
        masses = [float(q["mass"]) for q in peaks]
    if kind == "envelope":
        tol = [max(1.5, 1.5 * step), "Da"]
    else:
        tol = [10.0, "ppm"]
    # peak areas of the method: UniDec gives none for isotope envelopes (0), IsoDec gives the height again
    own = [(h, a) for q, h, a in zip(peaks, heights, areas) if not q.get("added")]
    real_areas = res.get("tag") != "isodec" and any(a > 0 for _, a in own) and \
        not all(a == h for h, a in own)
    return {"basis": "average" if kind == "envelope" else "monoisotopic", "kind": kind, "flags": flags,
            "masses": masses, "shown": list(masses) if shown is None else shown, "heights": heights,
            "areas": areas, "real_areas": real_areas, "tol": tol}


# ------------------------------------------------------------------ the reference calculation
def check_args(mass, height, area, shifts, tags, tol_da, tol_ppm, ref, kmax, max_extra, flags, min_rel=0.0):
    """The arguments as float64 arrays and numbers; ValueError for values the calculation does not accept
    (the library returns MS_BAD_ARG for the same ones)."""
    m = np.ascontiguousarray(mass, dtype=np.float64).ravel()
    n = len(m)
    if n > MAX_N:
        raise ValueError("at most %d masses" % MAX_N)
    h = np.ascontiguousarray(height, dtype=np.float64).ravel()
    a = np.zeros(n) if area is None else np.ascontiguousarray(area, dtype=np.float64).ravel()
    if len(h) != n or len(a) != n:
        raise ValueError("mass, height and area differ in length")
    if n and not (np.all(np.isfinite(m)) and np.all(m > 0)):
        raise ValueError("every mass must be a positive number")
    s = np.ascontiguousarray([] if shifts is None else shifts, dtype=np.float64).ravel()
    t = np.ascontiguousarray([] if tags is None else tags, dtype=np.float64).ravel()
    if len(s) > MAX_SHIFTS:
        raise ValueError("at most %d shifts" % MAX_SHIFTS)
    if len(t) > MAX_TAGS:
        raise ValueError("at most %d tags" % MAX_TAGS)
    if len(s) and not (np.all(np.isfinite(s)) and np.all(s != 0) and np.all(np.abs(s) <= 1e5)):
        raise ValueError("a shift must be a number other than 0, at most 100000 Da")
    if len(t) and not (np.all(np.isfinite(t)) and np.all(t > 0) and np.all(t <= 1e6)):
        raise ValueError("a tag mass must be above 0, at most 1000000 Da")
    tol_da, tol_ppm = float(tol_da), float(tol_ppm)
    if not (math.isfinite(tol_da) and math.isfinite(tol_ppm)) or tol_da < 0 or tol_ppm < 0 or \
            tol_da > 1000 or tol_ppm > 1e5 or not (tol_da > 0 or tol_ppm > 0):
        raise ValueError("the tolerance must be above 0 (at most 1000 Da and 100000 ppm)")
    ref, kmax, max_extra, flags = int(ref), int(kmax), int(max_extra), int(flags)
    if ref < -1 or ref >= n:
        raise ValueError("ref must be -1 or the index of a mass")
    if not (1 <= kmax <= MAX_K):
        raise ValueError("kmax must be 1 to %d" % MAX_K)
    if not (0 <= max_extra <= 2):
        raise ValueError("max_extra must be 0 to 2")
    if flags < 0 or flags > 3:
        raise ValueError("unknown flags")
    min_rel = float(min_rel)
    if not (0.0 <= min_rel <= 1.0):
        raise ValueError("min_rel must be 0 to 1")
    return m, h, a, s, t, tol_da, tol_ppm, ref, kmax, max_extra, flags, min_rel


def _candidates(s, t, kmax, max_extra, isos):
    """Step 4: (value, terms, sign, tag, k, s1, s2, iso) of every explanation."""
    out = []
    ns, nt = len(s), len(t)
    for sign in (1, -1):
        for iso in isos:
            for a in range(ns):
                base = float(s[a])
                v = base if sign > 0 else -base
                if iso:
                    v = v + iso * ISO
                out.append((v, 1 + abs(iso), sign, -1, 0, a, -1, iso))
            for ti in range(nt):
                for k in range(1, kmax + 1):
                    extras = [(-1, -1)]
                    if max_extra >= 1:
                        extras += [(a, -1) for a in range(ns)]
                    if max_extra >= 2:
                        # (two shifts that cancel, e.g. water and loss of water, are no explanation)
                        extras += [(a, b) for a in range(ns) for b in range(a, ns)
                                   if abs(float(s[a]) + float(s[b])) > 1e-3]
                    for a, b in extras:
                        base = k * float(t[ti])
                        if a >= 0:
                            base = base + float(s[a])
                        if b >= 0:
                            base = base + float(s[b])
                        v = base if sign > 0 else -base
                        if iso:
                            v = v + iso * ISO
                        out.append((v, 1 + (a >= 0) + (b >= 0) + abs(iso), sign, ti, k, a, b, iso))
    out.sort(key=lambda c: c[0])
    return out


def _key(c, e, tol):
    """Order of explanations (step 5): cost = terms + error / (0.5 tol), then the error, sign +, tag, k, s1,
    s2, iso."""
    return (c[1] + e / (0.5 * tol), e, 0 if c[2] > 0 else 1, c[3], c[4], c[5], c[6], c[7])


def _explain(cands, values, d, tol):
    """Step 5: (best candidate, its error, n_alt, (alt candidate, its error) or None) or None."""
    margin = 1e-9 * (abs(d) + 1.0)
    i0 = bisect.bisect_left(values, d - tol - margin)
    i1 = bisect.bisect_right(values, d + tol + margin)
    found = []
    for c in cands[i0:i1]:
        e = abs(d - c[0])
        if e <= tol:
            found.append((_key(c, e, tol), c))
    if not found:
        return None
    found.sort(key=lambda kc: kc[0])
    best = found[0][1]
    others = [c for _, c in found[1:] if abs(c[0] - best[0]) > 1e-6]  # of the same value: the same explanation
    alt = others[0] if others else None
    return best, d - best[0], len(others), (alt, d - alt[0]) if alt else None


def _canon(tc, sc, iso, s, t):
    """Value of a composition (step 7): tags, then shifts, in index order, then the isotope offset."""
    v = 0.0
    for i, c in enumerate(tc):
        if c:
            v = v + c * float(t[i])
    for i, c in enumerate(sc):
        if c:
            v = v + c * float(s[i])
    if iso:
        v = v + iso * ISO
    return v


def _terms(tc, sc, iso):
    return sum(1 for c in tc if c) + sum(abs(c) for c in sc) + abs(iso)


def find_shifts_py(mass, height, area, shifts, tags, tol_da, tol_ppm, ref=-1, kmax=4, max_extra=2, flags=0,
                   min_rel=0.0):
    """The Python reference of find_shifts (see the module's description)."""
    m, h, a, s, t, tol_da, tol_ppm, ref, kmax, max_extra, flags, min_rel = check_args(
        mass, height, area, shifts, tags, tol_da, tol_ppm, ref, kmax, max_extra, flags, min_rel)
    n, ns, nt = len(m), len(s), len(t)
    ppm = tol_ppm * 1e-6

    def tol_of(x, y):
        return tol_da + ppm * (x if x > y else y)

    # 1. order, heights and areas
    order = sorted(range(n), key=lambda i: (float(m[i]), i))
    pos = [0] * n
    for p_, i in enumerate(order):
        pos[i] = p_
    hh = [float(v) if math.isfinite(v) and v > 0 else 0.0 for v in h]
    aa = [float(v) if math.isfinite(v) and v > 0 else 0.0 for v in a]
    sh, sa = list(hh), list(aa)  # species heights and areas
    root = list(range(n))
    iso_pairs = []
    # 2. isotope peaks
    if flags & COLLAPSE:
        for pj in range(n):
            j = order[pj]
            mj = float(m[j])
            best = None
            for pi in range(pj - 1, -1, -1):
                i = order[pi]
                d = mj - float(m[i])
                tl = tol_of(float(m[i]), mj)
                if d > 2.0 * ISO + tl:
                    break
                for nn in (1, 2):
                    e = abs(d - nn * ISO)
                    if e <= tl:
                        key = (nn, e, pj - pi)
                        if best is None or key < best[0]:
                            best = (key, i)
                        break
            if best is not None:
                r0 = root[best[1]]
                root[j] = r0
                if hh[j] > sh[r0]:
                    sh[r0] = hh[j]
                sa[r0] = sa[r0] + aa[j]
                d = mj - float(m[r0])
                nn = int(round(d / ISO))
                iso_pairs.append({"from": r0, "to": j, "delta": d, "value": nn * ISO, "error": d - nn * ISO,
                                  "kind": 3, "sign": 1, "tag": -1, "k": 0, "s1": -1, "s2": -1, "iso": nn,
                                  "n_alt": 0, "alt": None, "ref": False, "link": False})
    # 3. masses below min_rel x the tallest species are left out; the reference
    top = max([sh[i] for i in order if root[i] == i] + [0.0])
    low = set(i for i in order if root[i] == i and sh[i] < min_rel * top)
    r = -1
    if ref >= 0:
        r = root[ref]
        low.discard(r)
    else:
        species = [i for i in order if root[i] == i]
        if species:
            r = min(species, key=lambda i: (-sh[i], i))
    out_species = []
    for i in range(n):
        out_species.append({"status": 4 if root[i] != i else (5 if i in low else 3),
                            "parent": root[i] if root[i] != i else -1,
                            "link": -1, "value": float("nan"), "error": float("nan"),
                            "comp": (tuple([0] * nt), tuple([0] * ns), 0), "height": sh[i] if root[i] == i else hh[i],
                            "area": sa[i] if root[i] == i else aa[i]})
    if r < 0:
        return {"ref": -1, "pairs": iso_pairs, "species": out_species,
                "conj": [{"height": [float("nan")] * (kmax + 1), "area": [float("nan")] * (kmax + 1),
                          "avg_height": float("nan"), "avg_area": float("nan"), "n": 0} for _ in range(nt)]}
    mr = float(m[r])
    # 4. explanations
    cands = _candidates(s, t, kmax, max_extra, (0, -1, 1) if flags & ISOTOPES else (0,))
    values = [c[0] for c in cands]
    # 6. explained pairs (both directions)
    sp = [i for i in order if root[i] == i and i not in low]
    expl = {}  # (from, to) -> explanation of mass[to] - mass[from]
    nbrs = {i: [] for i in sp}  # the species each species has an explained pair with, in sorted order
    pair_list = []  # (lighter, heavier, explained) in the order of the output
    for x in range(len(sp)):
        ia = sp[x]
        ma = float(m[ia])
        for y in range(x + 1, len(sp)):
            ib = sp[y]
            mb = float(m[ib])
            tl = tol_of(ma, mb)
            if abs(mb - ma) <= tl:
                continue
            e1 = _explain(cands, values, mb - ma, tl)
            e2 = _explain(cands, values, ma - mb, tl) if e1 is not None else None
            if e1 is None or e2 is None:
                if ia == r or ib == r:
                    pair_list.append((ia, ib, False))
                continue
            expl[(ia, ib)] = e1
            expl[(ib, ia)] = e2
            nbrs[ia].append(ib)
            nbrs[ib].append(ia)
            pair_list.append((ia, ib, True))

    def counts(ex):
        cb = ex[0]
        tc, sc = [0] * nt, [0] * ns
        if cb[3] >= 0:
            tc[cb[3]] += cb[2] * cb[4]
        if cb[5] >= 0:
            sc[cb[5]] += cb[2]
        if cb[6] >= 0:
            sc[cb[6]] += cb[2]
        return tc, sc, cb[7]

    # 7. compositions, from the reference outwards
    rank = {r: 0}
    comp = {r: ([0] * nt, [0] * ns, 0)}
    status = {r: 0}
    out_species[r].update(status=0, value=0.0, error=0.0)
    best = {}  # species -> (key, comp, value, error, origin)

    def relax(i):
        tc0, sc0, io0 = comp[i]
        for j in nbrs[i]:
            if j in rank:
                continue
            tc1, sc1, io1 = counts(expl[(i, j)])
            tc = [x + y for x, y in zip(tc0, tc1)]
            sc = [x + y for x, y in zip(sc0, sc1)]
            io = io0 + io1
            if any(abs(x) > kmax for x in tc) or sum(abs(x) for x in sc) > max_extra + 1 or abs(io) > 1:
                continue  # more complex than a direct explanation can be
            mj = float(m[j])
            v = _canon(tc, sc, io, s, t)
            e = (mj - mr) - v
            tl = tol_of(mr, mj)
            if abs(e) > tl:
                continue
            key = (_terms(tc, sc, io) + abs(e) / (0.5 * tl), abs(e), rank[i])
            if j not in best or key < best[j][0]:
                best[j] = (key, (tc, sc, io), v, e, i)
    relax(r)
    while best:
        j = min(best, key=lambda q: (best[q][0][0], best[q][0][1], pos[q]))
        key, cj, v, e, i = best.pop(j)
        rank[j] = len(rank)
        comp[j] = cj
        status[j] = 1 if i == r else 2
        out_species[j].update(status=status[j], parent=i, value=v, error=e,
                              comp=(tuple(cj[0]), tuple(cj[1]), cj[2]))
        relax(j)
    for j in sp:
        if j not in status:
            status[j] = 3
    out_species[r]["comp"] = (tuple([0] * nt), tuple([0] * ns), 0)
    # 8. pairs, oriented from the parent
    pairs = []
    pair_at = {}  # (from, to) -> index into pairs
    for ia, ib, ok_ in pair_list:
        ra, rb = rank.get(ia), rank.get(ib)
        if ra is not None and rb is not None:
            p, c = (ia, ib) if ra < rb else (ib, ia)
        elif ra is not None:
            p, c = ia, ib
        elif rb is not None:
            p, c = ib, ia
        else:
            da, db = abs(float(m[ia]) - mr), abs(float(m[ib]) - mr)
            p, c = (ib, ia) if db < da else (ia, ib)
        d = float(m[c]) - float(m[p])
        isref = p == r or c == r
        if not ok_:
            if status[c] != 3:
                continue  # the species is explained through another mass
            pairs.append({"from": p, "to": c, "delta": d, "value": float("nan"), "error": float("nan"),
                          "kind": 2, "sign": 1, "tag": -1, "k": 0, "s1": -1, "s2": -1, "iso": 0,
                          "n_alt": 0, "alt": None, "ref": True, "link": False})
            continue
        cb, err, nalt, alt = expl[(p, c)]
        pr = {"from": p, "to": c, "delta": d, "value": cb[0], "error": err, "kind": 1 if cb[3] >= 0 else 0,
              "sign": cb[2], "tag": cb[3], "k": cb[4], "s1": cb[5], "s2": cb[6], "iso": cb[7], "n_alt": nalt,
              "alt": None, "ref": isref, "link": False}
        if alt is not None:
            ca, ea = alt
            pr["alt"] = {"tag": ca[3], "k": ca[4], "s1": ca[5], "s2": ca[6], "sign": ca[2], "iso": ca[7],
                         "error": ea}
        pair_at[(p, c)] = len(pairs)
        pairs.append(pr)
    # 9. links
    for j in sp:
        if j == r:
            continue
        mj = float(m[j])
        bl = None
        if status[j] in (1, 2):
            for i in nbrs[j]:
                if rank.get(i) is None or rank[i] >= rank[j]:
                    continue
                tc0, sc0, io0 = comp[i]
                tc1, sc1, io1 = counts(expl[(i, j)])
                if [x + y for x, y in zip(tc0, tc1)] == list(comp[j][0]) and \
                        [x + y for x, y in zip(sc0, sc1)] == list(comp[j][1]) and io0 + io1 == comp[j][2]:
                    key = (abs(mj - float(m[i])), pos[i])
                    if bl is None or key < bl[0]:
                        bl = (key, pair_at[(i, j)])
        else:
            for i in nbrs[j]:
                k_ = pair_at.get((i, j))
                if k_ is None:
                    continue
                pr = pairs[k_]
                tc1, sc1, io1 = counts(expl[(i, j)])
                key = (_terms(tc1, sc1, io1), abs(mj - float(m[i])), abs(pr["error"]), pos[i])
                if bl is None or key < bl[0]:
                    bl = (key, k_)
            if bl is None:
                for k_, pr in enumerate(pairs):
                    if pr["kind"] == 2 and pr["to"] == j:
                        bl = ((), k_)
                        break
        if bl is not None:
            pairs[bl[1]]["link"] = True
            out_species[j]["link"] = bl[1]
    pairs += iso_pairs
    # 10. degree of conjugation
    conj = []
    for ti in range(nt):
        H, A = [0.0] * (kmax + 1), [0.0] * (kmax + 1)
        cnt = 0
        for j in range(n):
            if status.get(j) not in (0, 1, 2):
                continue
            tc = comp[j][0]
            kk = tc[ti]
            if kk < 0 or kk > kmax or any(tc[x] for x in range(nt) if x != ti):
                continue
            H[kk] = H[kk] + sh[j]
            A[kk] = A[kk] + sa[j]
            cnt += 1
        res = {"n": cnt}
        for key, S in (("height", H), ("area", A)):
            tot = 0.0
            wk = 0.0
            for k in range(kmax + 1):
                tot = tot + S[k]
                wk = wk + k * S[k]
            if tot > 0:
                res[key] = [100.0 * S[k] / tot for k in range(kmax + 1)]
                res["avg_" + key] = wk / tot
            else:
                res[key] = [float("nan")] * (kmax + 1)
                res["avg_" + key] = float("nan")
        conj.append(res)
    # 11. reference for tags: the lightest of a tag series below the tallest species
    if ref < 0 and nt > 0:
        cand = []
        for j in range(n):
            if status.get(j) not in (0, 1, 2):
                continue
            tc, sc, io = comp[j]
            if all(x <= 0 for x in tc) and any(x < 0 for x in tc):
                cand.append((sum(tc), sum(abs(x) for x in sc) + abs(io), -sh[j], j))
        if cand:
            return find_shifts_py(m, h, a, s, t, tol_da, tol_ppm, min(cand)[3], kmax, max_extra, flags, min_rel)
    return {"ref": r, "pairs": pairs, "species": out_species, "conj": conj}


def _lib():
    """msengine_py when the C++ library is in use, else None (the Python code runs)."""
    try:
        import msengine_py
        return msengine_py if msengine_py.available() else None
    except Exception:
        return None


def find_shifts(mass, height, area, shifts, tags, tol_da, tol_ppm, ref=-1, kmax=4, max_extra=2, flags=0,
                min_rel=0.0):
    """The mass shifts of a list of masses (see the module's description): the library's ms_mass_shifts
    when it is there, else the Python reference."""
    args = check_args(mass, height, area, shifts, tags, tol_da, tol_ppm, ref, kmax, max_extra, flags, min_rel)
    E = _lib()
    if E is not None:
        try:
            r = E.mass_shifts(*args)
        except AttributeError:  # an msengine_py without the function
            r = None
        if r is not None:
            return r
    return find_shifts_py(*args)


# ------------------------------------------------------------------ text
def _signed(c, name):
    a = abs(c)
    return ("%d %s" % (a, name)) if a > 1 else name


def comp_text(comp, shift_names, tag_names, empty="reference"):
    """A composition as text, e.g. "PDI + water", "2 PDI + Na", "−Na", "PDI (isotope +1)"."""
    tc, sc, iso = comp
    parts = []
    for i, c in enumerate(tc):
        if c:
            parts.append((c, _signed(c, tag_names[i] if i < len(tag_names) else "tag %d" % (i + 1))))
    for i, c in enumerate(sc):
        if c:
            parts.append((c, _signed(c, shift_names[i] if i < len(shift_names) else "shift %d" % (i + 1))))
    if not parts:
        txt = empty if not iso else ""
    else:
        txt = ""
        for n_, (c, p) in enumerate(parts):
            if n_ == 0:
                txt = ("−" + p) if c < 0 else p
            else:
                txt += (" − " if c < 0 else " + ") + p
    if iso:
        txt = (txt + " " if txt else "") + "(isotope %+d)" % iso
    return txt


def pair_comp(pr, nt, ns):
    """The explanation of a pair (or its alternative dict) as a composition."""
    tc, sc = [0] * nt, [0] * ns
    sg = pr["sign"]
    if pr.get("tag", -1) >= 0:
        tc[pr["tag"]] += sg * pr["k"]
    if pr.get("s1", -1) >= 0:
        sc[pr["s1"]] += sg
    if pr.get("s2", -1) >= 0:
        sc[pr["s2"]] += sg
    return tuple(tc), tuple(sc), pr.get("iso", 0)

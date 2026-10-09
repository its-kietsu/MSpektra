"""
Exact mass tools for HRMS Analysis (MSpektra): formula parser, ion m/z,
isotope patterns and a formula finder.

Monoisotopic masses and abundances: IUPAC / NIST (AME 2012, CIAAW 2013).
"""
import re
import itertools

import numpy as np

ELECTRON = 0.00054857990946
PROTON = 1.00727646688

ISOTOPES = {
    "H": [(1.00782503207, 0.999885), (2.0141017778, 0.000115)],
    "Li": [(6.015122795, 0.0759), (7.01600455, 0.9241)],
    "B": [(10.0129370, 0.199), (11.0093054, 0.801)],
    "C": [(12.0, 0.9893), (13.0033548378, 0.0107)],
    "N": [(14.0030740048, 0.99636), (15.0001088982, 0.00364)],
    "O": [(15.99491461956, 0.99757), (16.99913170, 0.00038), (17.9991610, 0.00205)],
    "F": [(18.99840322, 1.0)],
    "Na": [(22.9897692809, 1.0)],
    "Mg": [(23.985041700, 0.7899), (24.98583692, 0.1000), (25.982592929, 0.1101)],
    "Al": [(26.98153863, 1.0)],
    "Si": [(27.9769265325, 0.92223), (28.976494700, 0.04685), (29.97377017, 0.03092)],
    "P": [(30.97376163, 1.0)],
    "S": [(31.97207100, 0.9499), (32.97145876, 0.0075), (33.96786690, 0.0425), (35.96708076, 0.0001)],
    "Cl": [(34.96885268, 0.7576), (36.96590259, 0.2424)],
    "K": [(38.96370668, 0.932581), (39.96399848, 0.000117), (40.96182576, 0.067302)],
    "Ca": [(39.96259098, 0.96941), (41.95861801, 0.00647), (42.9587666, 0.00135), (43.9554818, 0.02086)],
    "Mn": [(54.9380451, 1.0)],
    "Fe": [(53.9396105, 0.05845), (55.9349375, 0.91754), (56.9353940, 0.02119), (57.9332756, 0.00282)],
    "Co": [(58.9331950, 1.0)],
    "Ni": [(57.9353429, 0.680769), (59.9307864, 0.262231), (60.9310560, 0.011399), (61.9283451, 0.036345),
           (63.9279660, 0.009256)],
    "Cu": [(62.9295975, 0.6915), (64.9277895, 0.3085)],
    "Zn": [(63.9291422, 0.4917), (65.9260334, 0.2773), (66.9271273, 0.0404), (67.9248442, 0.1845),
           (69.9253193, 0.0061)],
    "Se": [(73.9224764, 0.0089), (75.9192136, 0.0937), (76.9199140, 0.0763), (77.9173091, 0.2377),
           (79.9165213, 0.4961), (81.9166994, 0.0873)],
    "Br": [(78.9183371, 0.5069), (80.9162906, 0.4931)],
    "Ru": [(95.907598, 0.0554), (97.905287, 0.0187), (98.9059393, 0.1276), (99.9042195, 0.1260),
           (100.9055821, 0.1706), (101.9043493, 0.3155), (103.905433, 0.1862)],
    "Rh": [(102.905504, 1.0)],
    "Pd": [(101.905609, 0.0102), (103.904036, 0.1114), (104.905085, 0.2233), (105.903486, 0.2733),
           (107.903892, 0.2646), (109.905153, 0.1172)],
    "Ag": [(106.905097, 0.51839), (108.904752, 0.48161)],
    "Sn": [(111.904818, 0.0097), (113.902779, 0.0066), (114.903342, 0.0034), (115.901741, 0.1454),
           (116.902952, 0.0768), (117.901603, 0.2422), (118.903308, 0.0859), (119.9021947, 0.3258),
           (121.9034390, 0.0463), (123.9052739, 0.0579)],
    "I": [(126.904473, 1.0)],
    "Cs": [(132.905451933, 1.0)],
    "Ir": [(190.9605940, 0.373), (192.9629264, 0.627)],
    "Pt": [(189.959932, 0.00012), (191.9610380, 0.00782), (193.9626803, 0.3286), (194.9647911, 0.3378),
           (195.9649515, 0.2521), (197.9678930, 0.07356)],
    "Au": [(196.9665687, 1.0)],
    "Hg": [(195.965833, 0.0015), (197.9667690, 0.0997), (198.9682799, 0.1687), (199.9683260, 0.2310),
           (200.9703023, 0.1318), (201.9706430, 0.2986), (203.9734939, 0.0687)],
}
for _el, _iso in ISOTOPES.items():
    _iso.sort(key=lambda t: -t[1])  # most abundant first (monoisotopic for light elements)
# monoisotopic mass: the most abundant isotope of each element (IUPAC)
MONO = {el: max(iso, key=lambda t: t[1])[0] for el, iso in ISOTOPES.items()}

ADDUCTS = [
    # label, formula to add, formula to remove, electrons removed (+) or added (-), charge
    ("[M+H]+", "H", "", 1, 1),
    ("[M+Na]+", "Na", "", 1, 1),
    ("[M+K]+", "K", "", 1, 1),
    ("[M+NH4]+", "NH4", "", 1, 1),
    ("[M+2H]2+", "H2", "", 2, 2),
    ("[M+H+Na]2+", "HNa", "", 2, 2),
    ("[M]+ (cation or radical)", "", "", 1, 1),
    ("[M-H]-", "", "H", -1, -1),
    ("[M+Cl]-", "Cl", "", -1, -1),
    ("[M+HCOO]-", "CHO2", "", -1, -1),
    ("[M+CH3COO]-", "C2H3O2", "", -1, -1),
    ("[M-2H]2-", "", "H2", -2, -2),
    ("[M]- (anion or radical)", "", "", -1, -1),
    ("M (neutral)", "", "", 0, 0),
]
ADDUCT_LABELS = [a[0] for a in ADDUCTS]


def _lib():
    """msengine_py when the C++ library is in use: the patterns, the formula finder and the
    pattern matching run there and the Python code is the fallback (the reference)."""
    try:
        import msengine_py
        return msengine_py if msengine_py.available() else None
    except Exception:
        return None

_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)|(\()|(\))(\d*)")


def parse_formula(text):
    """'C6H12O6', 'Ca(OH)2', 'C2H5OH', 'CuSO4·5H2O' (a part after · , . or *
    may start with a count) -> {'C': 6, ...}. Raises ValueError."""
    text = (text or "").replace(" ", "")
    if not text:
        raise ValueError("Empty formula")
    parts = re.split(r"[·•⋅*.]", text)
    if len(parts) > 1:
        out = {}
        for part in parts:
            m = re.match(r"(\d+)(.*)$", part)
            k, body = (int(m.group(1)), m.group(2)) if m else (1, part)
            if not body:
                raise ValueError("Cannot read the formula '%s'" % text)
            for e, c in _parse_one(body).items():
                out[e] = out.get(e, 0) + c * k
        return out
    return _parse_one(text)


def _parse_one(text):
    stack = [{}]
    pos = 0
    for m in _TOKEN.finditer(text):
        if m.start() != pos:
            raise ValueError("Cannot read the formula at '%s'" % text[pos:])
        pos = m.end()
        el, n, op, cl, cn = m.groups()
        if el:
            if el not in ISOTOPES:
                raise ValueError("Unknown element %s" % el)
            stack[-1][el] = stack[-1].get(el, 0) + (int(n) if n else 1)
        elif op:
            stack.append({})
        elif cl:
            if len(stack) < 2:
                raise ValueError("Unbalanced parentheses")
            grp = stack.pop()
            k = int(cn) if cn else 1
            for e, c in grp.items():
                stack[-1][e] = stack[-1].get(e, 0) + c * k
    if pos != len(text) or len(stack) != 1:
        raise ValueError("Cannot read the formula '%s'" % text)
    return stack[0]


def format_formula(f):
    order = ["C", "H"] + sorted(e for e in f if e not in ("C", "H"))
    if "C" not in f:
        order = sorted(f)
    out = ""
    for e in order:
        n = f.get(e, 0)
        if n > 0:
            out += e + (str(n) if n > 1 else "")
    return out


def mono_mass(f):
    return float(sum(MONO[e] * n for e, n in f.items()))


def ion_formula(f, adduct):
    lab, add, rem, el, z = ADDUCTS[adduct] if isinstance(adduct, int) else next(a for a in ADDUCTS if a[0] == adduct)
    g = dict(f)
    if add:
        for e, n in parse_formula(add).items():
            g[e] = g.get(e, 0) + n
    if rem:
        for e, n in parse_formula(rem).items():
            g[e] = g.get(e, 0) - n
            if g[e] < 0:
                raise ValueError("The formula has too few %s atoms for %s" % (e, lab))
    return g, el, z


def ion_mz(f, adduct):
    g, el, z = ion_formula(f, adduct)
    m = mono_mass(g) - el * ELECTRON
    return m / abs(z) if z else m


def isotope_pattern(f, adduct=None, min_rel=1e-4, merge=None):
    """Isotope pattern of the ion (m/z, relative intensity %, max = 100),
    grouped into the nominal isotope peaks (M, M+1, ...) as seen at the
    resolution of a QTOF. merge=None groups by nominal mass; a value in m/z
    keeps the fine structure that is further apart than this."""
    if adduct is not None:
        g, el, z = ion_formula(f, adduct)
    else:
        g, el, z = dict(f), 0, 0
    E = _lib()
    if E is not None and (merge is None or merge > 0):
        # the elements in the order of the dict: the order of the sums below
        text = "".join("%s%d" % (e, int(n)) for e, n in g.items() if n > 0)
        r = E.isotope_pattern(text, -1.0 if merge is None else merge, min_rel) if text else None
        if r is not None and len(r[0]):
            arr = np.column_stack(r)
            arr[:, 0] = arr[:, 0] - el * ELECTRON
            if z:
                arr[:, 0] /= abs(z)
            return arr
    peaks = {0: (0.0, 1.0)}  # key: nominal offset -> (mass, prob) coarse
    dist = [(0.0, 1.0)]
    for e, n in g.items():
        if n <= 0:
            continue
        iso = ISOTOPES[e]
        for _ in range(int(n)):
            new = {}
            for m, p in dist:
                for im, ip in iso:
                    mm = m + im
                    key = round(mm * 1e4)  # 0.1 mDa fine structure
                    q = p * ip
                    if key in new:
                        om, op = new[key]
                        new[key] = ((om * op + mm * q) / (op + q), op + q)
                    else:
                        new[key] = (mm, q)
            top = max(v[1] for v in new.values())
            dist = [v for v in new.values() if v[1] >= top * min_rel * 0.01]
    dist.sort()
    mono = dist[0][0] if dist else 0.0
    mono_ref = mono_mass(g)
    groups = {}
    for m, p in dist:
        k = int(round((m - mono_ref) / 1.00235)) if merge is None else round(m / merge)
        a = groups.get(k)
        groups[k] = (a[0] + m * p, a[1] + p) if a else (m * p, p)
    out = sorted((v[0] / v[1], v[1]) for v in groups.values())
    arr = np.array(out)
    arr[:, 0] = arr[:, 0] - el * ELECTRON
    if z:
        arr[:, 0] /= abs(z)
    arr[:, 1] = arr[:, 1] / arr[:, 1].max() * 100.0
    return arr[arr[:, 1] >= min_rel * 100.0]


def rdb(f):
    """Ring and double bond equivalents (C, Si: 4-valent; N, P: 3; H,
    halogens, alkali: 1)."""
    tetra = f.get("C", 0) + f.get("Si", 0)
    tri = f.get("N", 0) + f.get("P", 0) + f.get("B", 0)
    mono = sum(f.get(e, 0) for e in ("H", "F", "Cl", "Br", "I", "Na", "K", "Li"))
    return tetra - mono / 2.0 + tri / 2.0 + 1.0


def find_formulas(mz, adduct=0, ppm=5.0, limits=None, observed=None, max_results=30):
    """Candidate formulas for an ion. limits: {'C': (0, 60), 'H': (0, 120),
    'N': (0, 10), 'O': (0, 20), 'S': (0, 3), 'P': (0, 2), 'F': ..., ...}.
    observed: optional isotope peaks [(m/z, rel %)] of the measured ion used
    to rank candidates. Returns a list of dicts sorted by score."""
    lab, add, rem, el, z = ADDUCTS[adduct]
    zabs = abs(z) or 1
    # [M]+ and [M]-: M is the ion itself, an even electron ion (pyridinium,
    # quaternary ammonium, carboxylate of a zwitterion: half integer RDB) or a
    # radical ion (integer RDB); the other ions are made from a neutral,
    # even electron molecule M (integer RDB)
    ion_itself = z != 0 and not add and not rem
    ion_mass = mz * zabs + el * ELECTRON
    addf = parse_formula(add) if add else {}
    remf = parse_formula(rem) if rem else {}
    neutral = ion_mass - mono_mass(addf) + mono_mass(remf)
    tol = ion_mass * ppm * 1e-6
    lim = {"C": (0, 80), "H": (0, 160), "N": (0, 10), "O": (0, 20), "S": (0, 3), "P": (0, 2)}
    if limits:
        lim.update(limits)
    hetero = [e for e in lim if e not in ("C", "H") and lim[e][1] > 0]
    ranges = [range(lim[e][0], lim[e][1] + 1) for e in hetero]
    mC, mH = MONO["C"], MONO["H"]
    out = []
    E = _lib()
    res = None
    if E is not None:
        try:
            res = E.find_formulas(mz, neutral, tol, zabs, el, list(lim), [int(lim[e][0]) for e in lim],
                                  [int(lim[e][1]) for e in lim], add, rem, ion_itself)
        except (TypeError, ValueError):
            res = None
    if res is not None:
        counts, theo, errs, rdbs = res
        for k in range(len(theo)):
            f = {"C": int(counts[k, 0]), "H": int(counts[k, 1])}
            f.update({e: int(n) for e, n in zip(hetero, counts[k, 2:]) if n})
            out.append({"formula": format_formula(f), "f": f, "mz": float(theo[k]), "ppm": float(errs[k]),
                        "rdb": float(rdbs[k])})
    for combo in (itertools.product(*ranges) if res is None else ()):  # the Python fallback
        rest = sum(MONO[e] * n for e, n in zip(hetero, combo))
        if rest > neutral + tol:
            continue
        cmax = min(lim["C"][1], int((neutral - rest) / mC) + 1)
        for c in range(lim["C"][0], cmax + 1):
            h = int(round((neutral - rest - c * mC) / mH))
            if h < lim["H"][0] or h > lim["H"][1]:
                continue
            m = rest + c * mC + h * mH
            if abs(m - neutral) > tol:
                continue
            f = {"C": c, "H": h}
            f.update({e: n for e, n in zip(hetero, combo) if n})
            r = rdb(f)
            if r < -0.5 or abs(r * 2 - round(r * 2)) > 1e-6:
                continue
            if c > 0:  # plausibility (Kind and Fiehn, seven golden rules, loose)
                # the lower bound counts the halogens with H (perfluorinated and perchlorinated compounds)
                hx = h + f.get("F", 0) + f.get("Cl", 0) + f.get("Br", 0) + f.get("I", 0)
                if not (0.1 <= hx / c and h / c <= 6.0) or f.get("N", 0) / c > 4 or f.get("O", 0) / c > 3:
                    continue
            elif h > 8 and not f.get("N"):
                continue
            # integer RDB for the neutral molecule (even electron)
            if abs(r - round(r)) > 1e-6 and not ion_itself:
                continue
            ion = dict(f)
            for e, n in addf.items():
                ion[e] = ion.get(e, 0) + n
            ok = True
            for e, n in remf.items():
                ion[e] = ion.get(e, 0) - n
                ok = ok and ion[e] >= 0
            if not ok:
                continue
            theo = (mono_mass(ion) - el * ELECTRON) / zabs
            err = (mz - theo) / theo * 1e6
            out.append({"formula": format_formula(f), "f": f, "mz": theo, "ppm": err, "rdb": r})
    if observed is not None and len(observed) >= 2:
        obs = np.asarray(observed, float)
        for c in out:
            pat = isotope_pattern(c["f"], adduct)
            c["iso"] = isotope_match(pat, obs)
    for c in out:
        s = abs(c["ppm"]) / max(ppm, 1e-9)
        if "iso" in c:
            s = s * 0.5 + (100.0 - c["iso"]) / 100.0 * 1.5
        c["score"] = s
    out.sort(key=lambda c: c["score"])
    return out[:max_results]


def isotope_match(pattern, observed):
    """0 to 100: agreement of the relative intensities of the isotope peaks
    M, M+1, M+2, M+3 of a theoretical pattern and measured peaks. The
    measured peaks start at the monoisotopic peak; the pattern is aligned on
    its peak nearest to it (for B, Fe, Li, Se, Pt, Hg, Sn, Ru or Pd the
    pattern has lighter peaks before the monoisotopic one)."""
    pat = np.asarray(pattern, float)
    obs = np.asarray(observed, float)
    if not len(pat) or not len(obs):
        return 0.0
    s = int(np.argmin(np.abs(pat[:, 0] - obs[0, 0])))
    pat = pat[s:s + 4]
    E = _lib()
    if E is not None:
        r = E.isotope_match(pat, obs)
        if r is not None:
            return r
    mono = pat[0, 0]
    o = []
    for m, a in pat:
        d = np.abs(obs[:, 0] - m)
        k = int(np.argmin(d))
        o.append(obs[k, 1] if d[k] <= max(0.02, m * 10e-6) else 0.0)
    o = np.array(o, float)
    if o[0] <= 0:
        return 0.0
    o = o / o[0] * 100.0
    t = pat[:, 1] / pat[0, 1] * 100.0
    w = np.maximum(t, 5.0)
    diff = np.sum(np.abs(o - t) / w * np.minimum(t, 100)) / np.sum(np.minimum(t, 100))
    return float(max(0.0, 100.0 * (1.0 - diff)))


def measured_pattern(spec, mono_mz, z=1, n=4, ppm=10.0):
    """Heights of the isotope peaks M, M+1/z, ... in a spectrum."""
    spec = np.asarray(spec, float)
    E = _lib()
    if E is not None and spec.ndim == 2 and len(spec):
        r = E.measured_pattern(spec, mono_mz, z, n, ppm)
        if r is not None:
            return r
    out = []
    for k in range(n):
        m = mono_mz + k * 1.00335 / max(abs(z), 1)
        tol = max(m * ppm * 1e-6, 0.004)
        w = (spec[:, 0] >= m - tol * 3) & (spec[:, 0] <= m + tol * 3)
        if not np.any(w):
            out.append((m, 0.0))
            continue
        i = np.where(w)[0][np.argmax(spec[w, 1])]
        out.append((float(spec[i, 0]), float(spec[i, 1])))
    return out

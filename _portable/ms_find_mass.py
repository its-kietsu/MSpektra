"""Find mass (LCMS and HRMS Analysis): where in the run does the charge state pattern of a neutral mass appear?

1. the charge states of M with m/z >= mzmin (500) inside the m/z range of the scan event, each a mass chromatogram of
   the window (M +- w + z*carrier)/z (negative ions: the carrier is a lost proton, (M +- w - z*1.007276)/z);
2. each state against the same state of decoy masses (M +- 50..170 Da, moved out by w - 10 for windows wider than
   10 Da, so that they stay outside the window): a state is present where it is 3x the median of its decoys;
   candidate scans have >= 4 neighbouring states present;
3. in each candidate scan the summed signal of the states (each +- min(w, 5) Da) over the masses M-150-w..M+150+w
   must peak inside the window, that peak the best within +-30 Da of itself and > 3x the median; every such peak
   inside the window (>= 20 % of the best) is a mass found there (a wide window can hold several species);
4. found regions: runs of >= 8 found scans (gaps of 2 allowed);
5. trace: the sum over the states of (state - decoy median), clipped at 0;
6. harmonic note: M x 1/2, 2, 3/2, 2/3 whose charge states share most of M's m/z while M's own states are weak.

HRMS Analysis (step 0.1): the same search on the 1 Da grid; the mass of a found scan is then the highest point of
the profile on a 0.1 Da grid within 1 Da of the 1 Da peak (close to the average mass of a protein whose isotopes are
resolved), and the mass of a region their median to 0.1 Da. With a calibration in use the windows and spectra are in
calibrated m/z, as everywhere in the window.

The mass chromatograms and the scans are read by msengine (one pass for all windows: ms_xic_multi), and the profiles
of step 3 are computed there for a block of candidate scans at a time (ms_mass_profile, the values of _profile); the
rest is a few array operations on them."""
import numpy as np

PROTON = 1.007276
CARRIERS = (("H+", 1.007276), ("Na+", 22.989218), ("NH4+", 18.033823), ("K+", 38.963158))
MZMIN = 500.0
DECOYS = (-170, -130, -110, -70, -50, 50, 70, 110, 130, 170)
GRID = np.arange(-150, 151, 1.0)
NEED = 4  # neighbouring states present in a candidate scan
BLOCK = 128  # candidate scans per ms_mass_profile call (cancel is checked between the calls)
HARMONICS = ((0.5, "half"), (2.0, "double"), (1.5, "one and a half times"), (2 / 3, "two thirds"))


class Cancelled(Exception):
    pass


def carrier_mass(name, polarity="+"):
    """Mass added per charge: the carrier for positive ions, minus a proton for negative ions."""
    if polarity == "-":
        return -PROTON
    return dict(CARRIERS).get(name, PROTON)


def charges(M, lo, hi, c=PROTON, mzmin=MZMIN):
    lo = max(lo or 0.0, mzmin)
    hi = hi if hi else float("inf")
    return [z for z in range(1, 400) if lo <= (M + z * c) / z <= hi]


def _windows(M, w, zs, c):
    out = []
    for z in zs:
        a, b = (M - w + z * c) / z, (M + w + z * c) / z
        out.append((0.5 * (a + b), 0.5 * (b - a)))
    return out


def _xics(f, e, wins, check):
    """Mass chromatograms of the windows [(centre, half width)] as rows (the values of f.chromatogram)."""
    try:
        import msengine_py as E
        Y = E.xic_multi(f, e, [m for m, h in wins], [h for m, h in wins]) if wins else None
    except Exception:
        Y = None
    if Y is not None:
        return np.asarray(Y, float)
    out = []
    for m, h in wins:
        check()
        out.append(np.asarray(f.chromatogram(e, "xic", m, h)[1], float))
    return np.array(out).reshape(len(wins), -1)


def decoys(w):
    """Decoy mass offsets: DECOYS, moved out by w - 10 Da for windows wider than 10 Da (outside the window)."""
    sh = max(0.0, float(w) - 10.0)
    return tuple(d + (sh if d > 0 else -sh) for d in DECOYS)


def grid(w):
    """Mass offsets of the profile: GRID, wider by the window beyond 0 (whole Da)."""
    r = 150 + int(np.ceil(max(0.0, float(w))))
    return GRID if r <= 150 + 5 else np.arange(-r, r + 1, 1.0)


def fine_grid(w, step):
    """Mass offsets in steps of step (< 1 Da) over every offset a found peak of the 1 Da grid can lie within 1 Da
    of (|offset| <= max(w, 3) + 1)."""
    n = int(round((np.ceil(max(float(w), 3.0)) + 1.0) / step))
    return np.round(np.arange(-n, n + 1) * step, 6)


def _finer(G, Gf, Pf, k):
    """The offset of the highest point of the fine profile Pf (offsets Gf) within 1 Da of G[k]."""
    sel = np.abs(Gf - G[k]) <= 1.0 + 1e-9
    j = int(np.argmax(np.where(sel, Pf, -np.inf)))
    return float(Gf[j])


def _profiles(f, e, times, M, w, zs, c, G):
    """_profile of the scans nearest to times (the spectra of f.scan_spectrum) in one msengine call: (rows, number
    of points of each spectrum), or None (then _profile on each spectrum)."""
    try:
        import msengine_py as E
        return E.mass_profile(f, e, times, M, min(float(w), 5.0), zs, c, G)
    except Exception:
        return None


def _profile(mz, it, M, w, zs, c, G=GRID):
    """Summed signal of the charge states over the masses M + G (each state +- min(w, 5) Da: a wide window must not
    blur two species into one)."""
    w = min(float(w), 5.0)
    cs = np.concatenate([[0.0], np.cumsum(it, dtype=float)])
    zz = np.asarray(zs, float)[:, None]
    cen = (M + G[None, :] + zz * c) / zz
    hi = cs[np.searchsorted(mz, cen + w / zz)]
    lo = cs[np.searchsorted(mz, cen - w / zz)]
    P = np.zeros(len(G))
    for k in range(len(zs)):
        P += hi[k] - lo[k]
    return P


def _species(lists, main, tol=3.0, share=0.3):
    """Masses seen in at least share of the found scans of a region (each within tol Da of their median), the main
    mass first: a wide window can hold several species (b-lactoglobulin A and B in 18300 +- 100)."""
    flat = sorted(m for ms in lists for m in ms)
    groups = []
    for m in flat:
        if groups and m - groups[-1][-1] <= tol:
            groups[-1].append(m)
        else:
            groups.append([m])
    out = [float(np.median(gr)) for gr in groups if len(gr) >= share * max(1, len(lists))]
    out.sort(key=lambda m: abs(m - main))
    return out


def _runs(mask, min_len=8, gap=2):
    idx = np.where(mask)[0]
    if not len(idx):
        return []
    segs, a, b = [], idx[0], idx[0]
    for i in idx[1:]:
        if i - b <= gap + 1:
            b = i
        else:
            segs.append((a, b))
            a = b = i
    segs.append((a, b))
    return [(int(a), int(b)) for a, b in segs if b - a + 1 >= min_len]


def search(f, e, M, w, carrier="H+", mzmin=MZMIN, harmonics=True, cancel=None, progress=None, step=1.0):
    """Search scan event e of the file f (EngineFile) for the neutral mass M (window +- w Da). cancel() -> True
    stops it (Cancelled); progress(fraction). step < 1 (HRMS): the masses found to that step (see above).
    Returns a dict of plain values and arrays (see the keys at the end)."""
    def check():
        if cancel is not None and cancel():
            raise Cancelled()

    def tell(x):
        if progress is not None:
            try:
                progress(float(x))
            except Exception:
                pass
    M, w = float(M), float(w)
    ev = f.events[e]
    pol = ev.get("polarity") or "+"
    c = carrier_mass(carrier, pol)
    zs = charges(M, ev.get("mz_low"), ev.get("mz_high"), c, mzmin)
    t = np.asarray(f.chromatogram(e, "tic")[0], float)
    res = dict(M=M, w=w, e=int(e), carrier=carrier if pol != "-" else "H+", polarity=pol, c=float(c),
               mzmin=float(mzmin), t=t, trace=np.zeros(len(t)), ok=np.zeros(len(t), bool), zs=list(zs),
               regions=[], harmonic=None, few=len(zs) < NEED)
    if len(zs) < NEED:
        tell(1.0)
        return res
    check()
    wins = _windows(M, w, zs, c)
    dec = decoys(w)
    for d in dec:
        wins += _windows(M + d, w, zs, c)
    Y = _xics(f, e, wins, check)
    nz = len(zs)
    X = Y[:nz]
    Dz = Y[nz:].reshape(len(dec), nz, -1)
    tell(0.3)
    check()
    ref = np.median(Dz, axis=0)
    present = X > 3 * ref
    run = np.zeros(X.shape[1], int)
    cur = np.zeros(X.shape[1], int)
    for i in range(nz):
        cur = np.where(present[i], cur + 1, 0)
        run = np.maximum(run, cur)
    trace = (X - ref).clip(0).sum(axis=0)
    ok = np.zeros(len(t), bool)
    best = np.full(len(t), np.nan)
    peaks = {}  # scan -> masses found in it (several when the window holds several species)
    cand = np.where(run >= NEED)[0]
    G = grid(w)
    near = np.abs(G) <= max(w, 3)
    fine = float(step) < 1.0
    Gf = fine_grid(w, float(step)) if fine else np.zeros(0)
    GG = np.concatenate([G, Gf]) if fine else G  # (one profile over both grids)
    ng = len(G)
    binw = getattr(f, "_binw", None)  # (the bin width in use: the search does not change it)
    rows = None
    for j, i in enumerate(cand):
        if j % 16 == 0:
            check()
            tell(0.3 + 0.6 * j / max(1, len(cand)))
        if j % BLOCK == 0:
            rows = _profiles(f, e, t[cand[j:j + BLOCK]], M, w, zs, c, GG)
        if rows is not None:
            if not rows[1][j % BLOCK]:
                continue
            P = rows[0][j % BLOCK]
        else:
            s = f.scan_spectrum(e, float(t[i]), binw=binw) if binw is not None else f.scan_spectrum(e, float(t[i]))
            s = np.asarray(s[0] if isinstance(s, tuple) else s)
            if not len(s):
                continue
            P = _profile(s[:, 0], s[:, 1], M, w, zs, c, GG)
        if fine:
            P, Pf = P[:ng], P[ng:]
        med = np.median(P)
        k = int(np.argmax(np.where(near, P, -1)))
        if P[k] >= P[np.abs(G - G[k]) <= 30].max() - 1e-9 and P[k] > 3 * med:
            ok[i] = True
            best[i] = M + (_finer(G, Gf, Pf, k) if fine else G[k])
            if w > 10:
                js = [j for j in np.where(near & (P >= 0.2 * P[k]) & (P > 3 * med))[0]
                      if P[j] >= P[np.abs(G - G[j]) <= 30].max() - 1e-9]
                peaks[int(i)] = [float(M + (_finer(G, Gf, Pf, j) if fine else G[j])) for j in js]
    regions = []
    for a, b in _runs(ok):
        sel = slice(a, b + 1)
        states = [int(z) for z, p in zip(zs, present[:, sel].mean(axis=1)) if p >= 0.5]
        g = dict(t0=float(t[a]), t1=float(t[b]), mass=float(np.nanmedian(best[sel])),
                 height=float(trace[sel].max()), states=states)
        if peaks:
            g["masses"] = _species([peaks[i] for i in range(a, b + 1) if i in peaks], g["mass"])
        if fine:  # (to the step: a median of two values ends in 5)
            nd = max(0, int(round(-np.log10(float(step)))))
            g["mass"] = round(g["mass"], nd)
            if "masses" in g:
                g["masses"] = [round(m, nd) for m in g["masses"]]
        regions.append(g)
    res.update(trace=trace, ok=ok, regions=regions)
    tell(0.9)
    if harmonics and regions:
        check()
        res["harmonic"] = _harmonic(f, e, M, w, zs, c, X, t, regions, ev, mzmin, check)
    tell(1.0)
    return res


def _harmonic(f, e, M, w, zs, c, X, t, regions, ev, mzmin, check):
    """A related mass whose charge states include most of M's (M's own states, not shared, must be weak and its
    own states strong): (mass, name, states shared) or None."""
    sel = np.zeros(len(t), bool)
    for g in regions:
        sel |= (t >= g["t0"]) & (t <= g["t1"])
    mz_m = np.array([(M + z * c) / z for z in zs])
    for r, name in HARMONICS:
        M2 = M * r
        z2 = charges(M2, ev.get("mz_low"), ev.get("mz_high"), c, mzmin)
        if len(z2) < 3:
            continue
        mz_2 = np.array([(M2 + z * c) / z for z in z2])
        sh_m = np.array([np.min(np.abs(mz_2 - m)) < max(0.3, w / z) for m, z in zip(mz_m, zs)])
        sh_2 = np.array([np.min(np.abs(mz_m - m)) < max(0.3, w * r / z) for m, z in zip(mz_2, z2)])
        if sh_m.sum() < 3 or (~sh_2).sum() < 2:
            continue
        check()
        X2 = _xics(f, e, _windows(M2, w * r, z2, c), check)
        sig, sig2 = X[:, sel].mean(axis=1), X2[:, sel].mean(axis=1)
        own_m = sig[~sh_m].mean() if (~sh_m).any() else 0.0  # M's states that M2 does not explain
        own_2 = sig2[~sh_2].mean()  # M2's states that M does not explain
        if own_m < 0.3 * sig[sh_m].mean() and own_2 > 0.3 * sig[sh_m].mean():
            return (float(M2), name, int(sh_m.sum()))
    return None


# ------------------------------------------------------------------ display
def merged(regions, gap=0.15):
    """Found regions closer than gap (min) as one."""
    out = []
    for g in sorted(regions, key=lambda g: g["t0"]):
        if out and g["t0"] - out[-1]["t1"] <= gap:
            o = out[-1]
            o["t1"] = max(o["t1"], g["t1"])
            if g["height"] > o["height"]:
                o.update(height=g["height"], mass=g["mass"], states=list(g["states"]))
        else:
            out.append(dict(g, states=list(g["states"])))
    return out


def core(res, g, frac):
    """The main part of a region: where the trace is >= frac of the region's top (t0, t1)."""
    t, y = res["t"], res["trace"]
    sel = (t >= g["t0"]) & (t <= g["t1"])
    if not np.any(sel):
        return g["t0"], g["t1"]
    on = sel & (y >= frac * y[sel].max())
    tt = t[on]
    return (float(tt.min()), float(tt.max())) if len(tt) else (g["t0"], g["t1"])


def shaded(res):
    """The parts shaded: [(t0, t1, region)] of the regions >= 5 % of the trace's top, where the trace is >= 10 %
    of the region's top."""
    if not res.get("regions"):
        return []
    top = float(np.max(res["trace"])) if len(res["trace"]) else 0.0
    out = []
    for g in merged(res["regions"]):
        if g["height"] >= 0.05 * top:
            a, b = core(res, g, 0.1)
            out.append((a, b, g))
    return out


def best_region(res):
    regs = merged(res.get("regions") or [])
    return max(regs, key=lambda g: g["height"]) if regs else None


def mz_of(res, z):
    return (res["M"] + z * res["c"]) / z

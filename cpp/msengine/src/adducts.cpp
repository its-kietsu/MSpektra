// Suggested neutral mass of a compound from its averaged ESI+ and ESI- spectra of a unit resolution
// quadrupole LC-MS (Shimadzu LCMS-2020 style: profile spectra binned at 0.05 to 0.1 m/z, masses good to
// about 0.2 to 0.3). Port of the reference in ms_adducts.py: the same steps, the same arithmetic in the
// same order, so that both give the same numbers bit for bit (tests/check_adducts.py compares them).
//
//  1. Points: those with a NaN or infinite m/z or intensity are left out, negative intensities count as
//     zero, the rest is sorted by m/z (stable sort). Base peak: the highest intensity; none above zero: no
//     peaks.
//  2. Maxima: runs of equal intensity whose neighbours on both sides are lower or missing; a neighbour
//     counts only when it lies within tol in m/z (a wider gap breaks the profile, so that centroid lists
//     and spectra without their zero points work as well). The maximum is the middle point of the run
//     ((first + last) / 2, rounded down), its height the intensity of the run. Kept: height > 0 and
//     height >= min_rel x base.
//  3. Apex: the intensity weighted mean m/z of the points of the peak, from the run outwards as long as
//     the next point is connected (within tol), not higher than the point before it and at least half the
//     height (sums from left to right).
//  4. Ions: the maxima ranked by height (descending), then position. A maximum is dropped when one of
//     better rank lies at d = (its apex - the other's apex) with |d| <= tol (the same peak: a jagged top
//     or a shoulder), |d - ISO/2| <= tol (isotope of a doubly charged ion) or |d - ISO| <= tol (isotope;
//     ISO = 1.00335483507, 13C - 12C). The first max_ions of the others are the ions of the spectrum, rel =
//     height / base.
//  5. Candidates: each ion with each adduct of its polarity gives M = (z x mz - d) / k (the ion's m/z is
//     (k M + d) / z; d with the exact masses of the ions: proton 1.00727646688, Na+, K+, NH4+, Cl-,
//     HCOO-, CH3COO- from the atomic masses minus or plus an electron). Only M > 0.
//  6. Groups: the candidates sorted by M (then ion: positive ions first, each polarity in rank order;
//     then adduct code); a candidate joins the group of the one before it when it is at most f x tol
//     heavier, f = the larger z / k of the two adducts (the error of M is z / k x the error of the m/z:
//     2 for doubly charged ions, 1, 0.5 for dimers). In a group each peak keeps only its most common
//     adduct (lowest code), then each adduct only its strongest ion (highest rel, then the first ion).
//     A group needs two ions or more after that. Its ions in the order: positive first, rel descending,
//     m/z, code. Mass = sum(rel x M) / sum(rel) over the ions in that order.
//  7. Score = sum(rel x w) in that order, w = 1 for [M+H]+, [M-H]- and [M+Na]+, 0.8 for the other singly
//     charged ions of one molecule, 0.5 for doubly charged and dimer ions; then x 1.5 when ions of both
//     polarities support the group, x 0.5 when it rests on doubly charged and dimer ions only. Not
//     normalized: a well supported compound scores about 1 or more, a coincidence of two weak peaks a few
//     hundredths.
//  8. Mobile phase additives: a candidate whose mass lies within max(tol, 0.5) of TFA (113.99286), formic
//     acid (46.00548), acetic acid (60.02113) or DFA (difluoroacetic acid, 96.00229, the LC-MS alternative
//     to TFA) is marked with that additive (the first in this order). At least 0.5: unit resolution
//     quadrupoles are often 0.3 to 0.4 off at low m/z (TFA [M+H]+ measured at 115.35).
//  9. M+2 partners (one Br or Cl; M2 = 1.997): the groups in the order of mass (stable); for each A not
//     removed, a heavier group B not removed is its partner when |mass B - mass A - M2| <= tol and every
//     ion of B at least half as high as B's strongest has an ion of A of the same adduct at m/z
//     k x M2 / z lower (within tol; the counterparts of weaker ions can be missing: below min_rel, beyond
//     the max_ions strongest, or lost in the isotope lump of a doubly charged ion); of several, the best
//     scored (then the lighter). A keeps its mass and ions, gets m2 = mass of B, the ions of B, m2_ratio =
//     rel of B's strongest ion (the first of the highest) / rel of its counterpart in A, and score A +
//     score B (also when B scored higher: the lighter, monoisotopic one stays); B is removed.
// 10. Reinterpretations: the groups by score (descending; then mass, stable); one whose peaks (its ions
//     and its partner's: the picked maxima) are all peaks of better groups kept is removed (it only
//     rereads ions already explained: 2M and M/2 readings, [M+K]+ / [M+Cl]- readings of an [M+H]+ /
//     [M-H]- pair). The peaks of the additive groups are background: a group other than an additive with
//     fewer than two peaks outside them is removed too (it rests on one peak and an ion of the mobile
//     phase).
// 11. No group left other than additives: every ion of no group kept alone, taken as [M+H]+ (positive)
//     or [M-H]- (negative), the adducts "assumed [M+H]+" / "assumed [M-H]-", w = 0.1 (score 0.1 x rel;
//     marked by step 8 too); step 9 among these (a single ion and its M+2 partner).
// 12. Additives: score x 0.1, and they come after all other candidates.
// 13. Sorted: additives last, then by score (descending), then mass (stable); the first top are returned.
// 14. Base peaks: per polarity the apex m/z of the strongest picked peak (always the first ion: nothing
//     ranks above it) and whether the best candidate explains it (it is one of its ions or of its
//     partner's ions).
#include "common.h"
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <vector>

namespace ms {
namespace {

const double PROTON = 1.00727646688;
const double ISO = 1.00335483507;
const double M2 = 1.997;   // 81Br - 79Br 1.99795, 37Cl - 35Cl 1.99705
const double NaN = std::numeric_limits<double>::quiet_NaN();

struct Adduct {
    const char* name;
    int pol, k, z;   // polarity, molecules, charges: m/z = (k M + d) / z
    double d, w;     // mass added (signed), weight in the score
};

// the code is the index; a lower code is the more common adduct
const Adduct ADD[MS_ADDUCT_COUNT] = {
    {"[M+H]+", 1, 1, 1, PROTON, 1.0},
    {"[M-H]-", -1, 1, 1, -PROTON, 1.0},
    {"[M+Na]+", 1, 1, 1, 22.989220702091, 1.0},
    {"[M+NH4]+", 1, 1, 1, 18.033825553441, 0.8},
    {"[M+HCOO]-", -1, 1, 1, 44.998202851279, 0.8},
    {"[M+K]+", 1, 1, 1, 38.963157906491, 0.8},
    {"[M+Cl]-", -1, 1, 1, 34.969401261909, 0.8},
    {"[M+CH3COO]-", -1, 1, 1, 59.013852915739, 0.8},
    {"[M+2H]2+", 1, 1, 2, 2.0 * PROTON, 0.5},
    {"[M-2H]2-", -1, 1, 2, -2.0 * PROTON, 0.5},
    {"[2M+H]+", 1, 2, 1, PROTON, 0.5},
    {"[2M+Na]+", 1, 2, 1, 22.989220702091, 0.5},
    {"[2M-H]-", -1, 2, 1, -PROTON, 0.5},
    {"assumed [M+H]+", 1, 1, 1, PROTON, 0.1},
    {"assumed [M-H]-", -1, 1, 1, -PROTON, 0.1},
};
const int N_REGULAR = 13;
const int ASSUMED_POS = 13, ASSUMED_NEG = 14;

// mobile phase additives (neutral monoisotopic masses)
struct Additive {
    const char* name;
    double mass;
};
const Additive ADDITIVES[MS_ADDITIVE_COUNT] = {
    {"TFA", 113.99286375956},          // CF3COOH
    {"formic acid", 46.0054793036},    // HCOOH
    {"acetic acid", 60.02112936806},   // CH3COOH
    {"DFA", 96.00228562906},           // CHF2COOH
};

struct Ion {
    int pol;
    double mz, rel;
};

struct Max {
    long pos, a, b;
    double h, apex;
};

// steps 1 to 4: the ions of one spectrum, in rank order
std::vector<Ion> spectrum_ions(const double* mz0, const double* it0, long n, int pol, double tol, double min_rel,
                               int max_ions) {
    std::vector<Ion> out;
    std::vector<double> x, y;
    x.reserve(n);
    y.reserve(n);
    for (long i = 0; i < n; i++) {
        const double a = mz0[i], b = it0[i];
        if (!std::isfinite(a) || !std::isfinite(b)) continue;
        x.push_back(a);
        y.push_back(b > 0 ? b : 0.0);
    }
    const long m = (long)x.size();
    if (m == 0) return out;
    bool sorted = true;
    for (long i = 1; i < m && sorted; i++)
        if (x[i] < x[i - 1]) sorted = false;
    if (!sorted) {
        std::vector<long> o(m);
        for (long i = 0; i < m; i++) o[i] = i;
        std::stable_sort(o.begin(), o.end(), [&x](long p, long q) { return x[p] < x[q]; });
        std::vector<double> xs(m), ys(m);
        for (long i = 0; i < m; i++) { xs[i] = x[o[i]]; ys[i] = y[o[i]]; }
        x.swap(xs);
        y.swap(ys);
    }
    double base = 0.0;
    for (double v : y)
        if (v > base) base = v;
    if (!(base > 0)) return out;
    const double thr = min_rel * base;

    // 2. maxima (runs of equal connected points)
    std::vector<Max> mx;
    for (long a = 0; a < m;) {
        long b = a;
        while (b + 1 < m && x[b + 1] - x[b] <= tol && y[b + 1] == y[b]) b++;
        const double h = y[a];
        const bool left = a == 0 || !(x[a] - x[a - 1] <= tol) || y[a - 1] < h;
        const bool right = b == m - 1 || !(x[b + 1] - x[b] <= tol) || y[b + 1] < h;
        if (left && right && h > 0 && h >= thr) mx.push_back({(a + b) / 2, a, b, h, 0.0});
        a = b + 1;
    }
    // 3. apex
    for (Max& q : mx) {
        long l = q.a, r = q.b;
        const double half = 0.5 * q.h;
        while (l > 0 && x[l] - x[l - 1] <= tol && y[l - 1] <= y[l] && y[l - 1] >= half) l--;
        while (r < m - 1 && x[r + 1] - x[r] <= tol && y[r + 1] <= y[r] && y[r + 1] >= half) r++;
        double sx = 0.0, sw = 0.0;
        for (long i = l; i <= r; i++) {
            sx += x[i] * y[i];
            sw += y[i];
        }
        q.apex = sx / sw;
    }
    // 4. rank, isotopes and shoulders
    const long k = (long)mx.size();
    std::vector<long> order(k), rank(k), by_apex(k);
    for (long i = 0; i < k; i++) order[i] = by_apex[i] = i;
    std::sort(order.begin(), order.end(), [&mx](long p, long q) {
        return mx[p].h > mx[q].h || (mx[p].h == mx[q].h && mx[p].pos < mx[q].pos);
    });
    for (long r = 0; r < k; r++) rank[order[r]] = r;
    std::sort(by_apex.begin(), by_apex.end(), [&mx](long p, long q) {
        return mx[p].apex < mx[q].apex || (mx[p].apex == mx[q].apex && p < q);
    });
    std::vector<double> apx(k);
    for (long i = 0; i < k; i++) apx[i] = mx[by_apex[i]].apex;
    for (long r = 0; r < k && (long)out.size() < max_ions; r++) {
        const Max& q = mx[order[r]];
        // the maxima that can lie in one of the windows (a wider range; the tests below decide)
        const long lo = (long)(std::lower_bound(apx.begin(), apx.end(), q.apex - (ISO + 2.0 * tol)) - apx.begin());
        const long hi = (long)(std::upper_bound(apx.begin(), apx.end(), q.apex + 2.0 * tol) - apx.begin());
        bool drop = false;
        for (long j = lo; j < hi && !drop; j++) {
            const long p = by_apex[j];
            if (rank[p] >= r) continue;
            const double d = q.apex - mx[p].apex;
            drop = std::fabs(d) <= tol || std::fabs(d - 0.5 * ISO) <= tol || std::fabs(d - ISO) <= tol;
        }
        if (!drop) out.push_back({pol, q.apex, q.h / base});
    }
    return out;
}

struct Cand {
    double M;
    int ion, code;
};

struct Group {
    std::vector<Cand> m, m2m;   // its ions; the ions of its M+2 partner
    double mass = 0.0, score = 0.0, m2 = NaN, m2_ratio = NaN;
    int both = 0, assumed = 0, additive = -1;
    bool removed = false;
};

// ions in their order (positive first, rel descending, m/z, code), mass, score, both, additive (step 8)
void finish(Group& g, const std::vector<Ion>& ions, double tol) {
    std::sort(g.m.begin(), g.m.end(), [&ions](const Cand& a, const Cand& b) {
        const Ion &p = ions[a.ion], &q = ions[b.ion];
        if (p.pol != q.pol) return p.pol > q.pol;
        if (p.rel != q.rel) return p.rel > q.rel;
        if (p.mz != q.mz) return p.mz < q.mz;
        return a.code < b.code;
    });
    double sm = 0.0, sr = 0.0, s = 0.0;
    bool pos = false, neg = false, single = false;
    for (const Cand& c : g.m) {
        const Ion& io = ions[c.ion];
        sm += io.rel * c.M;
        sr += io.rel;
        s += io.rel * ADD[c.code].w;
        (io.pol > 0 ? pos : neg) = true;
        if (ADD[c.code].k == 1 && ADD[c.code].z == 1) single = true;
    }
    g.mass = sm / sr;
    g.both = pos && neg;
    if (g.both) s = s * 1.5;
    if (!single) s = s * 0.5;
    g.score = s;
    g.additive = -1;
    for (int a = 0; a < MS_ADDITIVE_COUNT && g.additive < 0; a++)
        if (std::fabs(g.mass - ADDITIVES[a].mass) <= std::max(tol, 0.5)) g.additive = a;
}

// step 10: B (heavier) is the M+2 partner of A when every ion of B at least half as high as B's strongest
// has an ion of A of the same adduct k x M2 / z lower (within tol); the index in A.m of the counterpart of
// each ion of B (A.m.size(): none, a weaker ion)
bool m2_partner(const Group& A, const Group& B, const std::vector<Ion>& ions, double tol, std::vector<size_t>& cp) {
    if (!(B.mass > A.mass) || !(std::fabs(B.mass - A.mass - M2) <= tol)) return false;
    double top = ions[B.m[0].ion].rel;
    for (const Cand& m : B.m) top = std::max(top, ions[m.ion].rel);
    cp.assign(B.m.size(), A.m.size());
    for (size_t i = 0; i < B.m.size(); i++) {
        const Adduct& D = ADD[B.m[i].code];
        bool found = false;
        for (size_t j = 0; j < A.m.size() && !found; j++)
            if (A.m[j].code == B.m[i].code &&
                std::fabs(ions[B.m[i].ion].mz - ions[A.m[j].ion].mz - D.k * M2 / D.z) <= tol) {
                cp[i] = j;
                found = true;
            }
        if (!found && !(ions[B.m[i].ion].rel < 0.5 * top)) return false;
    }
    return true;
}

// step 9 among the candidates G[lo .. hi): in the order of mass, B joins A
void merge_m2(std::vector<Group>& G, long lo, long hi, const std::vector<Ion>& io, double tol) {
    std::vector<long> by_mass;
    for (long i = lo; i < hi; i++) by_mass.push_back(i);
    std::stable_sort(by_mass.begin(), by_mass.end(), [&G](long p, long q) { return G[p].mass < G[q].mass; });
    std::vector<size_t> cp, best_cp;
    for (long ia : by_mass) {
        Group& A = G[ia];
        if (A.removed) continue;
        long best = -1;
        for (long ib : by_mass) {
            if (ib == ia || G[ib].removed || !m2_partner(A, G[ib], io, tol, cp)) continue;
            if (best < 0 || G[ib].score > G[best].score || (G[ib].score == G[best].score && G[ib].mass < G[best].mass)) {
                best = ib;
                best_cp = cp;
            }
        }
        if (best < 0) continue;
        Group& B = G[best];
        size_t k = 0;   // B's strongest ion (the first of the highest)
        for (size_t i = 1; i < B.m.size(); i++)
            if (io[B.m[i].ion].rel > io[B.m[k].ion].rel) k = i;
        A.m2 = B.mass;
        A.m2_ratio = io[B.m[k].ion].rel / io[A.m[best_cp[k]].ion].rel;
        A.m2m = B.m;
        A.score = A.score + B.score;
        B.removed = true;
    }
}

struct Result {
    std::vector<Ion> ions;
    std::vector<Group> cands;   // best first
    double base_mz[2];
    int explained[2];
};

void neutral_masses(const double* pmz, const double* pit, long np, const double* nmz, const double* nit, long nn,
                    double tol, double min_rel, int max_ions, int top, Result& R) {
    std::vector<Ion>& io = R.ions;
    io = spectrum_ions(pmz, pit, np, 1, tol, min_rel, max_ions);
    const long npos = (long)io.size();
    const std::vector<Ion> neg = spectrum_ions(nmz, nit, nn, -1, tol, min_rel, max_ions);
    io.insert(io.end(), neg.begin(), neg.end());
    // 5. candidates
    std::vector<Cand> c;
    for (int i = 0; i < (int)io.size(); i++)
        for (int a = 0; a < N_REGULAR; a++) {
            const Adduct& A = ADD[a];
            if (A.pol != io[i].pol) continue;
            const double M = (A.z * io[i].mz - A.d) / A.k;
            if (M > 0) c.push_back({M, i, a});
        }
    std::sort(c.begin(), c.end(), [](const Cand& a, const Cand& b) {
        if (a.M != b.M) return a.M < b.M;
        if (a.ion != b.ion) return a.ion < b.ion;
        return a.code < b.code;
    });
    // 6. groups
    std::vector<Group> G;
    size_t s = 0;
    while (s < c.size()) {
        size_t e = s + 1;
        while (e < c.size()) {
            const Adduct &P = ADD[c[e - 1].code], &Q = ADD[c[e].code];
            if (!(c[e].M - c[e - 1].M <= tol * std::max((double)P.z / P.k, (double)Q.z / Q.k))) break;
            e++;
        }
        // each peak its most common adduct, then each adduct its strongest ion
        std::vector<Cand> a;
        for (size_t i = s; i < e; i++) {
            bool best = true;
            for (size_t j = s; j < e && best; j++)
                if (c[j].ion == c[i].ion && c[j].code < c[i].code) best = false;
            if (best) a.push_back(c[i]);
        }
        Group g;
        for (size_t i = 0; i < a.size(); i++) {
            bool best = true;
            const double ri = io[a[i].ion].rel;
            for (size_t j = 0; j < a.size() && best; j++) {
                if (j == i || a[j].code != a[i].code) continue;
                const double rj = io[a[j].ion].rel;
                if (rj > ri || (rj == ri && a[j].ion < a[i].ion)) best = false;
            }
            if (best) g.m.push_back(a[i]);
        }
        if (g.m.size() >= 2) {
            finish(g, io, tol);
            G.push_back(std::move(g));
        }
        s = e;
    }
    // 9. M+2 partners and 10. reinterpretations among the groups
    const long ng = (long)G.size();
    merge_m2(G, 0, ng, io, tol);
    std::vector<long> order;
    for (long i = 0; i < ng; i++)
        if (!G[i].removed) order.push_back(i);
    std::stable_sort(order.begin(), order.end(), [&G](long p, long q) {
        return G[p].score > G[q].score || (G[p].score == G[q].score && G[p].mass < G[q].mass);
    });
    std::vector<char> covered(io.size(), 0), background(io.size(), 0);
    for (const Group& g : G)
        if (!g.removed && g.additive >= 0) {
            for (const Cand& m : g.m) background[m.ion] = 1;
            for (const Cand& m : g.m2m) background[m.ion] = 1;
        }
    for (long i : order) {
        Group& g = G[i];
        std::vector<int> peaks;
        for (const Cand& m : g.m) peaks.push_back(m.ion);
        for (const Cand& m : g.m2m) peaks.push_back(m.ion);
        std::sort(peaks.begin(), peaks.end());
        peaks.erase(std::unique(peaks.begin(), peaks.end()), peaks.end());
        bool all = true;
        int own = 0;
        for (int p : peaks) {
            all = all && covered[p];
            own += !background[p];
        }
        if (all || (g.additive < 0 && own < 2)) {
            g.removed = true;
            continue;
        }
        for (int p : peaks) covered[p] = 1;
    }
    // 11. no group left other than additives: each ion of none alone, assumed [M+H]+ / [M-H]-, and
    // their M+2 partners
    bool analyte = false;
    for (const Group& g : G)
        if (!g.removed && g.additive < 0) analyte = true;
    if (!analyte) {
        for (int i = 0; i < (int)io.size(); i++) {
            if (covered[i]) continue;
            const int code = io[i].pol > 0 ? ASSUMED_POS : ASSUMED_NEG;
            const Adduct& A = ADD[code];
            const double M = (A.z * io[i].mz - A.d) / A.k;
            if (!(M > 0)) continue;
            Group g;
            g.m.push_back({M, i, code});
            g.assumed = 1;
            finish(g, io, tol);
            G.push_back(std::move(g));
        }
        merge_m2(G, ng, (long)G.size(), io, tol);
    }
    // 12. additives x 0.1 and after the others; 13. best first
    std::vector<Group>& out = R.cands;
    out.clear();
    for (Group& g : G)
        if (!g.removed) {
            if (g.additive >= 0) g.score = g.score * 0.1;
            out.push_back(std::move(g));
        }
    std::stable_sort(out.begin(), out.end(), [](const Group& a, const Group& b) {
        const bool pa = a.additive >= 0, pb = b.additive >= 0;
        if (pa != pb) return pb;
        return a.score > b.score || (a.score == b.score && a.mass < b.mass);
    });
    if ((long)out.size() > top) out.resize(top);
    // 14. base peaks (the first ion of each polarity) and whether the best candidate explains them
    const long first[2] = {npos > 0 ? 0 : -1, (long)io.size() > npos ? npos : -1};
    for (int p = 0; p < 2; p++) {
        R.base_mz[p] = first[p] >= 0 ? io[first[p]].mz : NaN;
        int ex = 0;
        if (first[p] >= 0 && !out.empty()) {
            for (const Cand& m : out[0].m) ex = ex || m.ion == first[p];
            for (const Cand& m : out[0].m2m) ex = ex || m.ion == first[p];
        }
        R.explained[p] = ex;
    }
}

void check_args(long np, const double* pmz, const double* pit, long nn, const double* nmz, const double* nit,
                double tol, double min_rel, int max_ions, int top, bool outs) {
    if (!outs || np < 0 || nn < 0 || (np > 0 && (!pmz || !pit)) || (nn > 0 && (!nmz || !nit)))
        throw std::invalid_argument("ms_neutral_masses: bad arguments");
    if (!(tol > 0 && tol <= 5)) throw std::invalid_argument("ms_neutral_masses: the tolerance must be above 0 and at most 5");
    if (!(min_rel >= 0 && min_rel <= 1)) throw std::invalid_argument("ms_neutral_masses: min_rel must be 0 to 1");
    if (max_ions < 1 || max_ions > 1000) throw std::invalid_argument("ms_neutral_masses: max_ions must be 1 to 1000");
    if (top < 1 || top > 1000) throw std::invalid_argument("ms_neutral_masses: top must be 1 to 1000");
}

thread_local Result g_res;
thread_local std::vector<ms_mass_ion> g_ions;
thread_local std::vector<ms_mass_cand2> g_cands2;
thread_local std::vector<ms_mass_cand> g_cands;

// the result in the flat arrays of the C interface
void flatten(const Result& R) {
    g_ions.clear();
    g_cands2.clear();
    for (const Group& g : R.cands) {
        ms_mass_cand2 mc;
        mc.mass = g.mass;
        mc.score = g.score;
        mc.both = g.both;
        mc.assumed = g.assumed;
        mc.additive = g.additive;
        mc.m2 = g.m2;
        mc.m2_ratio = g.m2_ratio;
        mc.first = (long)g_ions.size();
        mc.n = (long)g.m.size();
        for (const Cand& m : g.m) {
            const Ion& q = R.ions[m.ion];
            g_ions.push_back({q.pol, m.code, q.mz, q.rel, m.M});
        }
        mc.m2_first = (long)g_ions.size();
        mc.m2_n = (long)g.m2m.size();
        for (const Cand& m : g.m2m) {
            const Ion& q = R.ions[m.ion];
            g_ions.push_back({q.pol, m.code, q.mz, q.rel, m.M});
        }
        g_cands2.push_back(mc);
    }
}

}  // namespace
}  // namespace ms

extern "C" {

MS_API const char* ms_adduct_name(int code) {
    return code >= 0 && code < MS_ADDUCT_COUNT ? ms::ADD[code].name : "";
}

MS_API const char* ms_additive_name(int code) {
    return code >= 0 && code < MS_ADDITIVE_COUNT ? ms::ADDITIVES[code].name : "";
}

MS_API int ms_neutral_masses2(const double* pmz, const double* pit, long np, const double* nmz, const double* nit,
                              long nn, double tol, double min_rel, int max_ions, int top, const ms_mass_ion** ions,
                              long* nions, const ms_mass_cand2** cands, long* ncands, double* base_mz, int* explained) {
    return ms::guarded([&] {
        using namespace ms;
        check_args(np, pmz, pit, nn, nmz, nit, tol, min_rel, max_ions, top,
                   ions && nions && cands && ncands && base_mz && explained);
        g_ions.clear();
        g_cands2.clear();
        *ions = g_ions.data();
        *nions = 0;
        *cands = g_cands2.data();
        *ncands = 0;
        neutral_masses(pmz, pit, np, nmz, nit, nn, tol, min_rel, max_ions, top, g_res);
        flatten(g_res);
        *ions = g_ions.data();
        *nions = (long)g_ions.size();
        *cands = g_cands2.data();
        *ncands = (long)g_cands2.size();
        for (int p = 0; p < 2; p++) {
            base_mz[p] = g_res.base_mz[p];
            explained[p] = g_res.explained[p];
        }
    });
}

MS_API int ms_neutral_masses(const double* pmz, const double* pit, long np, const double* nmz, const double* nit,
                             long nn, double tol, double min_rel, int max_ions, int top, const ms_mass_ion** ions,
                             long* nions, const ms_mass_cand** cands, long* ncands) {
    return ms::guarded([&] {
        using namespace ms;
        check_args(np, pmz, pit, nn, nmz, nit, tol, min_rel, max_ions, top, ions && nions && cands && ncands);
        g_ions.clear();
        g_cands.clear();
        *ions = g_ions.data();
        *nions = 0;
        *cands = g_cands.data();
        *ncands = 0;
        neutral_masses(pmz, pit, np, nmz, nit, nn, tol, min_rel, max_ions, top, g_res);
        flatten(g_res);
        for (const ms_mass_cand2& m : g_cands2) g_cands.push_back({m.mass, m.score, m.both, m.assumed, m.first, m.n});
        *ions = g_ions.data();
        *nions = (long)g_ions.size();
        *cands = g_cands.data();
        *ncands = (long)g_cands.size();
    });
}

}  // extern "C"

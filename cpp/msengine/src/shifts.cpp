// Mass shift finder on the masses of a deconvolution result (MS Analysis 3.35): which masses differ by a
// known modification or adduct, by k tags plus at most two other shifts, which differences from the
// reference match nothing (unknown), and the degree of conjugation. Port of the reference in
// ms_shifts.py (find_shifts_py): the same steps, the same arithmetic in the same order, so that both give
// the same numbers (tests/check_shifts.py compares them). Compiled without contraction of multiply and
// add (CMakeLists.txt), as the reference computes.
//
//  1. Masses sorted (stable: by mass, then index). Heights and areas: not finite or negative count as 0.
//     The tolerance of two masses a, b: tol_da + tol_ppm x 1e-6 x max(a, b).
//  2. Isotope peaks (MS_SHIFT_COLLAPSE only): in the order of mass, a mass j is the isotope peak of a
//     lighter mass i when |m_j - m_i - n ISO| <= tol (n = 1 or 2, ISO = 1.00335483507); of several, the
//     smallest n, then the smallest error, then the nearest i. It belongs to the species of i (the first
//     mass of its chain): the species' height is the largest of its chain, its area the sum (in the order
//     of mass). Pairs of kind 3 (from the first mass of the chain) record it; isotope peaks take no
//     further part.
//  3. Species below min_rel x the height of the tallest species take no further part (status 5).
//     Reference: the given one (its species, also below min_rel), else the tallest species (the first in
//     input order of equal heights).
//  4. Explanations: sign x (one shift), or sign x (k x tag + up to max_extra shifts, a shift may repeat,
//     two shifts that cancel (|s1 + s2| <= 0.001 Da, e.g. water and loss of water) are left out), k = 1
//     .. kmax, one tag kind at a time; with MS_SHIFT_ISOTOPES each also + iso x ISO, iso = -1 or +1.
//     Value: base = k x tag (0 without a tag), base += s1, base += s2 (s1 <= s2), value = sign x base +
//     iso x ISO. Terms: 1 for the tags (any k), 1 per shift, 1 per isotope offset.
//  5. A difference d (tolerance tol) is explained by the candidates with |d - value| <= tol. The best has
//     the lowest cost = terms + |d - value| / (0.5 tol) (each half tolerance of error counts as one more
//     term), then the smaller error, then sign +, tag, k, s1, s2, iso (ascending). The alternative is the
//     next one in this order whose value differs from the best's by more than 1e-6 Da (one of the same
//     value is the same explanation: water and the reversed loss of water, a tag + double oxidation and a
//     tag + 2 oxidations); n_alt counts those of another value.
//  6. Explained pairs: every two species a, b (a lighter or equal) with |m_b - m_a| > tol whose
//     difference is explained (it then is in both directions; each direction has its own best
//     explanation: a loss as written going down, an addition going up).
//  7. Compositions, from the reference outwards (its composition is empty): a species j next to a placed
//     species i (an explained pair) can have the composition of i plus the explanation from i to j:
//     signed counts of each tag and shift and the isotope offset; its value is the sum counts x masses,
//     tags first, then the shifts, each in index order, then iso x ISO; its error (m_j - m_ref) - value
//     must be within the tolerance of (ref, j), every tag count within -kmax .. kmax, the shift counts at
//     most max_extra + 1 together, the isotope offset -1 to 1; its cost is terms (1 per tag kind used, 1
//     per shift, 1 per isotope) + |error| / (0.5 tol). Of the compositions a species can have, the cheapest is kept
//     (then the smaller error, then the parent placed first); repeatedly the species with the cheapest
//     one (then the smaller error, then its sorted position) is placed, until none is left: the others
//     are unknown.
//  8. Pairs: the explained pairs and the unknown pairs (the reference and an unknown species whose
//     difference nothing explains), oriented from the parent: of two placed species the one placed first,
//     a placed species before an unknown one, of two unknown ones the one nearer in mass to the reference
//     (equal: the lighter). Order: by the sorted position of the lighter, then of the heavier; the
//     isotope pairs follow.
//  9. Links (drawn by default): a placed species links to its nearest consistent parent (placed before
//     it, with an explained pair, its composition plus the explanation is the species' composition;
//     nearest in mass, then sorted position); an unknown species to the parent of its explained pair with
//     the fewest terms (then nearest, smaller error, sorted position), else to the reference with its
//     unknown pair.
// 10. Degree of conjugation, per tag: the species (reference and placed ones, in input order) with k of
//     this tag and no other tag, 0 <= k <= kmax, summed by height and by area: % = 100 x sum_k / total,
//     tags per molecule = (sum of k x sum_k) / total.
// 11. Reference for tags (only with an automatic reference, ref < 0, and at least one tag): when the
//     tallest species is a conjugate (common at high labelling), the unmodified species lies below it and
//     its tag counts are negative. Among the reference and the placed species (status 0 to 2), those with
//     every tag count <= 0 and at least one < 0: if there are any, the one with the lowest total tag count,
//     then the fewest other terms (|shift counts| + |isotope offset|), then the tallest (species height),
//     then input order, becomes the reference and steps 3 to 10 run again with it (ref_used reports it).
#include "common.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

namespace ms {
namespace {

const double ISO = 1.00335483507;  // 13C - 12C
const double NaN = std::numeric_limits<double>::quiet_NaN();
const long MAX_N = 2000, MAX_SHIFTS = 64, MAX_TAGS = 3;
const int MAX_K = 20;

struct Cand {
    double v;
    int terms, sign, tag, k, s1, s2, iso;
};

// the order of explanations (step 5): cost = terms + error / (0.5 tol), error, sign +, tag, k, s1, s2, iso
bool key_less(const Cand& a, double ea, const Cand& b, double eb, double tol) {
    const double ca = a.terms + ea / (0.5 * tol), cb = b.terms + eb / (0.5 * tol);
    if (ca != cb) return ca < cb;
    if (ea != eb) return ea < eb;
    const int sa = a.sign > 0 ? 0 : 1, sb = b.sign > 0 ? 0 : 1;
    if (sa != sb) return sa < sb;
    if (a.tag != b.tag) return a.tag < b.tag;
    if (a.k != b.k) return a.k < b.k;
    if (a.s1 != b.s1) return a.s1 < b.s1;
    if (a.s2 != b.s2) return a.s2 < b.s2;
    return a.iso < b.iso;
}

struct Explained {
    bool ok = false;
    long best = -1, alt = -1;  // indexes into the candidates
    double err = NaN, alt_err = NaN;
    int n_alt = 0;
};

struct Pair {
    long from, to;
    double delta, value, error;
    int kind, sign, tag, k, s1, s2, iso, n_alt;
    int alt_tag, alt_k, alt_s1, alt_s2, alt_sign, alt_iso;
    double alt_error;
    int ref, link;
};

struct Comp {
    std::vector<int> tc, sc;
    int iso = 0;
    bool operator==(const Comp& o) const { return tc == o.tc && sc == o.sc && iso == o.iso; }
};

class Finder {
public:
    Finder(const double* mass, const double* height, const double* area, long n, const double* shifts, long ns,
           const double* tags, long nt, double tol_da, double tol_ppm, long ref, int kmax, int max_extra,
           int flags, double min_rel)
        : m_(mass, mass + n), n_(n), s_(shifts, shifts + ns), t_(tags, tags + nt), ns_(ns), nt_(nt),
          tol_da_(tol_da), ppm_(tol_ppm * 1e-6), ref_in_(ref), kmax_(kmax), max_extra_(max_extra),
          flags_(flags), min_rel_(min_rel) {
        h_.resize(n);
        a_.resize(n);
        for (long i = 0; i < n; i++) {
            const double hv = height[i];
            h_[i] = std::isfinite(hv) && hv > 0 ? hv : 0.0;
            const double av = area ? area[i] : 0.0;
            a_[i] = std::isfinite(av) && av > 0 ? av : 0.0;
        }
    }

    void run();

    // results
    long r = -1;
    std::vector<Pair> pairs;
    std::vector<ms_shift_species> species;
    std::vector<int> comp;      // n x (nt + ns + 1)
    std::vector<double> conj;   // nt x (2 kmax + 5)

private:
    std::vector<double> m_, h_, a_;
    long n_;
    std::vector<double> s_, t_;
    long ns_, nt_;
    double tol_da_, ppm_;
    long ref_in_;
    int kmax_, max_extra_, flags_;
    double min_rel_;
    std::vector<Cand> cands_;
    std::vector<double> values_;

    double tol_of(double x, double y) const { return tol_da_ + ppm_ * (x > y ? x : y); }
    void candidates();
    Explained explain(double d, double tol) const;
    double canon(const Comp& c) const;
    static int terms(const Comp& c);
};

void Finder::candidates() {
    cands_.clear();
    std::vector<int> isos = {0};
    if (flags_ & MS_SHIFT_ISOTOPES) isos = {0, -1, 1};
    for (int sign : {1, -1}) {
        for (int iso : isos) {
            for (long a = 0; a < ns_; a++) {
                const double base = s_[a];
                double v = sign > 0 ? base : -base;
                if (iso) v = v + iso * ISO;
                cands_.push_back({v, 1 + std::abs(iso), sign, -1, 0, (int)a, -1, iso});
            }
            for (long ti = 0; ti < nt_; ti++) {
                for (int k = 1; k <= kmax_; k++) {
                    std::vector<std::pair<int, int>> extras = {{-1, -1}};
                    if (max_extra_ >= 1)
                        for (long a = 0; a < ns_; a++) extras.push_back({(int)a, -1});
                    if (max_extra_ >= 2)  // (two shifts that cancel, e.g. water and loss of water: none)
                        for (long a = 0; a < ns_; a++)
                            for (long b = a; b < ns_; b++)
                                if (std::fabs(s_[a] + s_[b]) > 1e-3) extras.push_back({(int)a, (int)b});
                    for (const auto& e : extras) {
                        double base = k * t_[ti];
                        if (e.first >= 0) base = base + s_[e.first];
                        if (e.second >= 0) base = base + s_[e.second];
                        double v = sign > 0 ? base : -base;
                        if (iso) v = v + iso * ISO;
                        cands_.push_back({v, 1 + (e.first >= 0) + (e.second >= 0) + std::abs(iso), sign, (int)ti, k,
                                          e.first, e.second, iso});
                    }
                }
            }
        }
    }
    std::stable_sort(cands_.begin(), cands_.end(), [](const Cand& x, const Cand& y) { return x.v < y.v; });
    values_.resize(cands_.size());
    for (size_t i = 0; i < cands_.size(); i++) values_[i] = cands_[i].v;
}

Explained Finder::explain(double d, double tol) const {
    Explained out;
    const double margin = 1e-9 * (std::fabs(d) + 1.0);
    const long i0 = (long)(std::lower_bound(values_.begin(), values_.end(), d - tol - margin) - values_.begin());
    const long i1 = (long)(std::upper_bound(values_.begin(), values_.end(), d + tol + margin) - values_.begin());
    std::vector<long> found;
    std::vector<double> errs;
    for (long i = i0; i < i1; i++) {
        const double e = std::fabs(d - cands_[i].v);
        if (e <= tol) {
            found.push_back(i);
            errs.push_back(e);
        }
    }
    if (found.empty()) return out;
    long bi = -1, ai = -1;
    for (size_t q = 0; q < found.size(); q++)
        if (bi < 0 || key_less(cands_[found[q]], errs[q], cands_[found[bi]], errs[bi], tol)) bi = (long)q;
    const double vb = cands_[found[bi]].v;
    int n_alt = 0;
    for (size_t q = 0; q < found.size(); q++) {
        if ((long)q == bi || !(std::fabs(cands_[found[q]].v - vb) > 1e-6)) continue;  // same value: same explanation
        n_alt++;
        if (ai < 0 || key_less(cands_[found[q]], errs[q], cands_[found[ai]], errs[ai], tol)) ai = (long)q;
    }
    out.ok = true;
    out.best = found[bi];
    out.err = d - cands_[out.best].v;
    out.n_alt = n_alt;
    if (ai >= 0) {
        out.alt = found[ai];
        out.alt_err = d - cands_[out.alt].v;
    }
    return out;
}

double Finder::canon(const Comp& c) const {
    double v = 0.0;
    for (long i = 0; i < nt_; i++)
        if (c.tc[i]) v = v + c.tc[i] * t_[i];
    for (long i = 0; i < ns_; i++)
        if (c.sc[i]) v = v + c.sc[i] * s_[i];
    if (c.iso) v = v + c.iso * ISO;
    return v;
}

int Finder::terms(const Comp& c) {
    int t = 0;
    for (int x : c.tc) t += x != 0;
    for (int x : c.sc) t += std::abs(x);
    return t + std::abs(c.iso);
}

void Finder::run() {
    const long n = n_;
    // 1. order
    std::vector<long> order(n);
    for (long i = 0; i < n; i++) order[i] = i;
    std::stable_sort(order.begin(), order.end(), [&](long x, long y) { return m_[x] < m_[y]; });
    std::vector<long> pos(n);
    for (long p = 0; p < n; p++) pos[order[p]] = p;
    std::vector<double> sh(h_), sa(a_);
    std::vector<long> root(n);
    for (long i = 0; i < n; i++) root[i] = i;
    std::vector<Pair> iso_pairs;
    // 2. isotope peaks
    if (flags_ & MS_SHIFT_COLLAPSE) {
        for (long pj = 0; pj < n; pj++) {
            const long j = order[pj];
            const double mj = m_[j];
            bool have = false;
            int bn = 0;
            double be = 0;
            long bd = 0, bi = -1;
            for (long pi = pj - 1; pi >= 0; pi--) {
                const long i = order[pi];
                const double d = mj - m_[i];
                const double tl = tol_of(m_[i], mj);
                if (d > 2.0 * ISO + tl) break;
                for (int nn = 1; nn <= 2; nn++) {
                    const double e = std::fabs(d - nn * ISO);
                    if (e <= tl) {
                        const long dist = pj - pi;
                        if (!have || nn < bn || (nn == bn && (e < be || (e == be && dist < bd)))) {
                            have = true;
                            bn = nn;
                            be = e;
                            bd = dist;
                            bi = i;
                        }
                        break;
                    }
                }
            }
            if (have) {
                const long r0 = root[bi];
                root[j] = r0;
                if (h_[j] > sh[r0]) sh[r0] = h_[j];
                sa[r0] = sa[r0] + a_[j];
                const double d = mj - m_[r0];
                const int nn = (int)std::nearbyint(d / ISO);
                Pair p{};
                p.from = r0;
                p.to = j;
                p.delta = d;
                p.value = nn * ISO;
                p.error = d - nn * ISO;
                p.kind = 3;
                p.sign = 1;
                p.tag = -1;
                p.k = 0;
                p.s1 = p.s2 = -1;
                p.iso = nn;
                p.n_alt = 0;
                p.alt_tag = p.alt_s1 = p.alt_s2 = -1;
                p.alt_k = p.alt_sign = p.alt_iso = 0;
                p.alt_error = NaN;
                p.ref = p.link = 0;
                iso_pairs.push_back(p);
            }
        }
    }
    // 3. masses below min_rel x the tallest species are left out; the reference
    double top = 0.0;
    for (long p = 0; p < n; p++)
        if (root[order[p]] == order[p] && sh[order[p]] > top) top = sh[order[p]];
    std::vector<char> low(n, 0);
    for (long i = 0; i < n; i++)
        if (root[i] == i && sh[i] < min_rel_ * top) low[i] = 1;
    r = -1;
    if (ref_in_ >= 0) {
        r = root[ref_in_];
        low[r] = 0;
    } else {
        for (long i = 0; i < n; i++) {
            if (root[i] != i) continue;
            if (r < 0 || sh[i] > sh[r]) r = i;
        }
    }
    const long nc = nt_ + ns_ + 1;
    species.assign(n, ms_shift_species{});
    comp.assign(n * nc, 0);
    for (long i = 0; i < n; i++) {
        ms_shift_species& q = species[i];
        q.status = root[i] != i ? 4 : (low[i] ? 5 : 3);
        q.parent = root[i] != i ? root[i] : -1;
        q.link = -1;
        q.value = NaN;
        q.error = NaN;
        q.height = root[i] == i ? sh[i] : h_[i];
        q.area = root[i] == i ? sa[i] : a_[i];
    }
    const long nb = 2 * kmax_ + 5;
    conj.assign(nt_ * nb, NaN);
    pairs.clear();
    if (r < 0) {
        for (long ti = 0; ti < nt_; ti++) conj[ti * nb + nb - 1] = 0.0;
        pairs = iso_pairs;
        return;
    }
    const double mr = m_[r];
    // 4. explanations
    candidates();
    // 6. explained pairs (both directions)
    std::vector<long> sp;
    for (long p = 0; p < n; p++)
        if (root[order[p]] == order[p] && !low[order[p]]) sp.push_back(order[p]);
    struct PairOf {
        long a, b;  // lighter, heavier (sorted positions)
        int ok;
    };
    std::vector<PairOf> plist;
    std::vector<std::vector<long>> nbrs(n);
    std::vector<Explained> exs;
    std::vector<int> ex_of((size_t)n * (size_t)n, -1);  // (from, to) -> explanation of mass[to] - mass[from]
    const long nsp = (long)sp.size();
    for (long x = 0; x < nsp; x++) {
        const long ia = sp[x];
        const double ma = m_[ia];
        for (long y = x + 1; y < nsp; y++) {
            const long ib = sp[y];
            const double mb = m_[ib];
            const double tl = tol_of(ma, mb);
            if (std::fabs(mb - ma) <= tl) continue;
            const Explained e1 = explain(mb - ma, tl);
            Explained e2;
            if (e1.ok) e2 = explain(ma - mb, tl);
            if (!e1.ok || !e2.ok) {
                if (ia == r || ib == r) plist.push_back({ia, ib, 0});
                continue;
            }
            ex_of[(size_t)ia * n + ib] = (int)exs.size();
            exs.push_back(e1);
            ex_of[(size_t)ib * n + ia] = (int)exs.size();
            exs.push_back(e2);
            nbrs[ia].push_back(ib);
            nbrs[ib].push_back(ia);
            plist.push_back({ia, ib, 1});
        }
    }
    auto counts_of = [&](long from, long to, Comp& out) {
        const Cand& cb = cands_[exs[ex_of[(size_t)from * n + to]].best];
        out.tc.assign(nt_, 0);
        out.sc.assign(ns_, 0);
        if (cb.tag >= 0) out.tc[cb.tag] += cb.sign * cb.k;
        if (cb.s1 >= 0) out.sc[cb.s1] += cb.sign;
        if (cb.s2 >= 0) out.sc[cb.s2] += cb.sign;
        out.iso = cb.iso;
    };
    // 7. compositions, from the reference outwards
    std::vector<long> rank(n, -1);
    std::vector<int> status(n, -1);
    std::vector<Comp> cmp(n);
    for (long i = 0; i < n; i++) {
        cmp[i].tc.assign(nt_, 0);
        cmp[i].sc.assign(ns_, 0);
    }
    long nrank = 0;
    rank[r] = nrank++;
    status[r] = 0;
    species[r].status = 0;
    species[r].value = 0.0;
    species[r].error = 0.0;
    struct Best {
        bool have = false;
        double cost = 0, ae = 0;
        long orank = 0, orig = -1;
        Comp c;
        double v = 0, e = 0;
    };
    std::vector<Best> best(n);
    Comp c1, c;
    auto relax = [&](long i) {
        for (long j : nbrs[i]) {
            if (rank[j] >= 0) continue;
            counts_of(i, j, c1);
            c.tc.resize(nt_);
            c.sc.resize(ns_);
            for (long q = 0; q < nt_; q++) c.tc[q] = cmp[i].tc[q] + c1.tc[q];
            for (long q = 0; q < ns_; q++) c.sc[q] = cmp[i].sc[q] + c1.sc[q];
            c.iso = cmp[i].iso + c1.iso;
            bool simple = std::abs(c.iso) <= 1;
            int nsh = 0;
            for (long q = 0; q < nt_; q++)
                if (std::abs(c.tc[q]) > kmax_) simple = false;
            for (long q = 0; q < ns_; q++) nsh += std::abs(c.sc[q]);
            if (!simple || nsh > max_extra_ + 1) continue;  // more complex than a direct explanation can be
            const double mj = m_[j];
            const double v = canon(c);
            const double e = (mj - mr) - v;
            const double tl = tol_of(mr, mj);
            if (std::fabs(e) > tl) continue;
            const double cost = terms(c) + std::fabs(e) / (0.5 * tl);
            const double ae = std::fabs(e);
            Best& b = best[j];
            if (!b.have || cost < b.cost || (cost == b.cost && (ae < b.ae || (ae == b.ae && rank[i] < b.orank)))) {
                b.have = true;
                b.cost = cost;
                b.ae = ae;
                b.orank = rank[i];
                b.orig = i;
                b.c = c;
                b.v = v;
                b.e = e;
            }
        }
    };
    relax(r);
    for (;;) {
        long j = -1;
        for (long q = 0; q < n; q++) {
            if (!best[q].have || rank[q] >= 0) continue;
            if (j < 0) {
                j = q;
                continue;
            }
            const Best& a = best[q];
            const Best& b = best[j];
            if (a.cost < b.cost || (a.cost == b.cost && (a.ae < b.ae || (a.ae == b.ae && pos[q] < pos[j])))) j = q;
        }
        if (j < 0) break;
        Best& b = best[j];
        rank[j] = nrank++;
        cmp[j] = b.c;
        status[j] = b.orig == r ? 1 : 2;
        species[j].status = status[j];
        species[j].parent = b.orig;
        species[j].value = b.v;
        species[j].error = b.e;
        b.have = false;
        relax(j);
    }
    for (long i : sp)
        if (status[i] < 0) {
            status[i] = 3;
            species[i].status = 3;
            species[i].parent = -1;
        }
    // 8. pairs, oriented from the parent
    std::vector<int> pair_at((size_t)n * (size_t)n, -1);
    for (const PairOf& L : plist) {
        const long ia = L.a, ib = L.b;
        const long ra = rank[ia], rb = rank[ib];
        long p, ch;
        if (ra >= 0 && rb >= 0) {
            p = ra < rb ? ia : ib;
            ch = ra < rb ? ib : ia;
        } else if (ra >= 0) {
            p = ia;
            ch = ib;
        } else if (rb >= 0) {
            p = ib;
            ch = ia;
        } else {
            const double da = std::fabs(m_[ia] - mr), db = std::fabs(m_[ib] - mr);
            p = db < da ? ib : ia;
            ch = db < da ? ia : ib;
        }
        Pair pr{};
        pr.from = p;
        pr.to = ch;
        pr.delta = m_[ch] - m_[p];
        pr.alt_tag = pr.alt_s1 = pr.alt_s2 = -1;
        pr.alt_k = pr.alt_sign = pr.alt_iso = 0;
        pr.alt_error = NaN;
        pr.ref = p == r || ch == r;
        pr.link = 0;
        if (!L.ok) {
            if (status[ch] != 3) continue;  // the species is explained through another mass
            pr.value = NaN;
            pr.error = NaN;
            pr.kind = 2;
            pr.sign = 1;
            pr.tag = -1;
            pr.k = 0;
            pr.s1 = pr.s2 = -1;
            pr.iso = 0;
            pr.n_alt = 0;
            pairs.push_back(pr);
            continue;
        }
        const Explained& ex = exs[ex_of[(size_t)p * n + ch]];
        const Cand& cb = cands_[ex.best];
        pr.value = cb.v;
        pr.error = ex.err;
        pr.kind = cb.tag >= 0 ? 1 : 0;
        pr.sign = cb.sign;
        pr.tag = cb.tag;
        pr.k = cb.k;
        pr.s1 = cb.s1;
        pr.s2 = cb.s2;
        pr.iso = cb.iso;
        pr.n_alt = ex.n_alt;
        if (ex.alt >= 0) {
            const Cand& ca = cands_[ex.alt];
            pr.alt_tag = ca.tag;
            pr.alt_k = ca.k;
            pr.alt_s1 = ca.s1;
            pr.alt_s2 = ca.s2;
            pr.alt_sign = ca.sign;
            pr.alt_iso = ca.iso;
            pr.alt_error = ex.alt_err;
        }
        pair_at[(size_t)p * n + ch] = (int)pairs.size();
        pairs.push_back(pr);
    }
    // 9. links
    for (long j : sp) {
        if (j == r) continue;
        const double mj = m_[j];
        bool have = false;
        long bk = -1;
        if (status[j] == 1 || status[j] == 2) {
            double bdist = 0;
            long bpos = 0;
            for (long i : nbrs[j]) {
                if (rank[i] < 0 || rank[i] >= rank[j]) continue;
                counts_of(i, j, c1);
                c.tc.resize(nt_);
                c.sc.resize(ns_);
                for (long q = 0; q < nt_; q++) c.tc[q] = cmp[i].tc[q] + c1.tc[q];
                for (long q = 0; q < ns_; q++) c.sc[q] = cmp[i].sc[q] + c1.sc[q];
                c.iso = cmp[i].iso + c1.iso;
                if (!(c == cmp[j])) continue;
                const double dist = std::fabs(mj - m_[i]);
                if (!have || dist < bdist || (dist == bdist && pos[i] < bpos)) {
                    have = true;
                    bdist = dist;
                    bpos = pos[i];
                    bk = pair_at[(size_t)i * n + j];
                }
            }
        } else {
            int bterms = 0;
            double bdist = 0, berr = 0;
            long bpos = 0;
            for (long i : nbrs[j]) {
                const int kp = pair_at[(size_t)i * n + j];
                if (kp < 0) continue;
                counts_of(i, j, c1);
                const int tt = terms(c1);
                const double dist = std::fabs(mj - m_[i]);
                const double ee = std::fabs(pairs[kp].error);
                if (!have || tt < bterms ||
                    (tt == bterms && (dist < bdist || (dist == bdist && (ee < berr || (ee == berr && pos[i] < bpos)))))) {
                    have = true;
                    bterms = tt;
                    bdist = dist;
                    berr = ee;
                    bpos = pos[i];
                    bk = kp;
                }
            }
            if (!have) {
                for (size_t q = 0; q < pairs.size(); q++) {
                    if (pairs[q].kind == 2 && pairs[q].to == j) {
                        have = true;
                        bk = (long)q;
                        break;
                    }
                }
            }
        }
        if (have && bk >= 0) {
            pairs[bk].link = 1;
            species[j].link = bk;
        }
    }
    for (long i = 0; i < n; i++) {
        if (status[i] == 1 || status[i] == 2) {
            for (long q = 0; q < nt_; q++) comp[i * nc + q] = cmp[i].tc[q];
            for (long q = 0; q < ns_; q++) comp[i * nc + nt_ + q] = cmp[i].sc[q];
            comp[i * nc + nc - 1] = cmp[i].iso;
        }
    }
    pairs.insert(pairs.end(), iso_pairs.begin(), iso_pairs.end());
    // 10. degree of conjugation
    for (long ti = 0; ti < nt_; ti++) {
        std::vector<double> H(kmax_ + 1, 0.0), A(kmax_ + 1, 0.0);
        long cnt = 0;
        for (long j = 0; j < n; j++) {
            if (status[j] < 0 || status[j] > 2) continue;
            const int kk = cmp[j].tc[ti];
            bool other = false;
            for (long x = 0; x < nt_; x++)
                if (x != ti && cmp[j].tc[x]) other = true;
            if (kk < 0 || kk > kmax_ || other) continue;
            H[kk] = H[kk] + sh[j];
            A[kk] = A[kk] + sa[j];
            cnt++;
        }
        double* out = &conj[ti * nb];
        for (int part = 0; part < 2; part++) {
            const std::vector<double>& S = part == 0 ? H : A;
            double tot = 0.0, wk = 0.0;
            for (int k = 0; k <= kmax_; k++) {
                tot = tot + S[k];
                wk = wk + k * S[k];
            }
            for (int k = 0; k <= kmax_; k++) out[part * (kmax_ + 1) + k] = tot > 0 ? 100.0 * S[k] / tot : NaN;
            out[2 * (kmax_ + 1) + part] = tot > 0 ? wk / tot : NaN;
        }
        out[nb - 1] = (double)cnt;
    }
}

void check_args(const double* mass, const double* height, long n, const double* shifts, long ns,
                const double* tags, long nt, double tol_da, double tol_ppm, long ref, int kmax, int max_extra,
                int flags, double min_rel, bool outs) {
    if (!outs) throw std::invalid_argument("ms_mass_shifts: an output pointer is NULL");
    if (n < 0 || n > MAX_N) throw std::invalid_argument("ms_mass_shifts: at most 2000 masses");
    if (n > 0 && (!mass || !height)) throw std::invalid_argument("ms_mass_shifts: mass or height is NULL");
    for (long i = 0; i < n; i++)
        if (!std::isfinite(mass[i]) || !(mass[i] > 0))
            throw std::invalid_argument("ms_mass_shifts: every mass must be a positive number");
    if (ns < 0 || ns > MAX_SHIFTS) throw std::invalid_argument("ms_mass_shifts: at most 64 shifts");
    if (nt < 0 || nt > MAX_TAGS) throw std::invalid_argument("ms_mass_shifts: at most 3 tags");
    if (ns > 0 && !shifts) throw std::invalid_argument("ms_mass_shifts: shifts is NULL");
    if (nt > 0 && !tags) throw std::invalid_argument("ms_mass_shifts: tags is NULL");
    for (long i = 0; i < ns; i++)
        if (!std::isfinite(shifts[i]) || shifts[i] == 0 || std::fabs(shifts[i]) > 1e5)
            throw std::invalid_argument("ms_mass_shifts: a shift must be a number other than 0, at most 100000 Da");
    for (long i = 0; i < nt; i++)
        if (!std::isfinite(tags[i]) || !(tags[i] > 0) || tags[i] > 1e6)
            throw std::invalid_argument("ms_mass_shifts: a tag mass must be above 0, at most 1000000 Da");
    if (!std::isfinite(tol_da) || !std::isfinite(tol_ppm) || tol_da < 0 || tol_ppm < 0 || tol_da > 1000 ||
        tol_ppm > 1e5 || !(tol_da > 0 || tol_ppm > 0))
        throw std::invalid_argument("ms_mass_shifts: the tolerance must be above 0 (at most 1000 Da and 100000 ppm)");
    if (ref < -1 || ref >= n) throw std::invalid_argument("ms_mass_shifts: ref must be -1 or the index of a mass");
    if (kmax < 1 || kmax > MAX_K) throw std::invalid_argument("ms_mass_shifts: kmax must be 1 to 20");
    if (max_extra < 0 || max_extra > 2) throw std::invalid_argument("ms_mass_shifts: max_extra must be 0 to 2");
    if (flags < 0 || flags > 3) throw std::invalid_argument("ms_mass_shifts: unknown flags");
    if (!(min_rel >= 0.0 && min_rel <= 1.0)) throw std::invalid_argument("ms_mass_shifts: min_rel must be 0 to 1");
}

thread_local std::vector<ms_shift_pair> g_pairs;
thread_local std::vector<ms_shift_species> g_species;
thread_local std::vector<int> g_comp;
thread_local std::vector<double> g_conj;

}  // namespace
}  // namespace ms

extern "C" {

MS_API int ms_mass_shifts(const double* mass, const double* height, const double* area, long n,
                          const double* shifts, long ns, const double* tags, long nt, double tol_da, double tol_ppm,
                          long ref, int kmax, int max_extra, int flags, double min_rel, const ms_shift_pair** pairs,
                          long* npairs, const ms_shift_species** species, const int** comp, const double** conj,
                          long* ref_used) {
    return ms::guarded([&] {
        using namespace ms;
        check_args(mass, height, n, shifts, ns, tags, nt, tol_da, tol_ppm, ref, kmax, max_extra, flags, min_rel,
                   pairs && npairs && species && comp && conj && ref_used);
        *npairs = 0;
        *ref_used = -1;
        g_pairs.clear();
        g_species.clear();
        g_comp.clear();
        g_conj.clear();
        *pairs = g_pairs.data();
        *species = g_species.data();
        *comp = g_comp.data();
        *conj = g_conj.data();
        Finder f0(mass, height, area, n, shifts, ns, tags, nt, tol_da, tol_ppm, ref, kmax, max_extra, flags, min_rel);
        f0.run();
        // 11. reference for tags: the lightest of a tag series below the tallest species
        std::unique_ptr<Finder> f1;
        if (ref < 0 && nt > 0) {
            const long nc = nt + ns + 1;
            long best = -1, bt = 0, bo = 0;
            double bh = 0;
            for (long j = 0; j < n; j++) {
                const int st = f0.species[j].status;
                if (st < 0 || st > 2) continue;
                const int* cj = &f0.comp[j * nc];
                long tot = 0, other = 0;
                bool neg = false, pos = false;
                for (long q = 0; q < nt; q++) {
                    tot += cj[q];
                    if (cj[q] > 0) pos = true;
                    if (cj[q] < 0) neg = true;
                }
                if (pos || !neg) continue;
                for (long q = 0; q < ns; q++) other += std::abs(cj[nt + q]);
                other += std::abs(cj[nc - 1]);
                const double hh = f0.species[j].height;
                if (best < 0 || tot < bt || (tot == bt && (other < bo || (other == bo && hh > bh)))) {
                    best = j;
                    bt = tot;
                    bo = other;
                    bh = hh;
                }
            }
            if (best >= 0) {
                f1.reset(new Finder(mass, height, area, n, shifts, ns, tags, nt, tol_da, tol_ppm, best, kmax, max_extra,
                                    flags, min_rel));
                f1->run();
            }
        }
        const Finder& f = f1 ? *f1 : f0;
        for (const Pair& p : f.pairs) {
            ms_shift_pair q;
            q.from = p.from;
            q.to = p.to;
            q.delta = p.delta;
            q.value = p.value;
            q.error = p.error;
            q.kind = p.kind;
            q.sign = p.sign;
            q.tag = p.tag;
            q.k = p.k;
            q.s1 = p.s1;
            q.s2 = p.s2;
            q.iso = p.iso;
            q.n_alt = p.n_alt;
            q.alt_tag = p.alt_tag;
            q.alt_k = p.alt_k;
            q.alt_s1 = p.alt_s1;
            q.alt_s2 = p.alt_s2;
            q.alt_sign = p.alt_sign;
            q.alt_iso = p.alt_iso;
            q.alt_error = p.alt_error;
            q.ref = p.ref;
            q.link = p.link;
            g_pairs.push_back(q);
        }
        g_species = f.species;
        g_comp = f.comp;
        g_conj = f.conj;
        *pairs = g_pairs.data();
        *npairs = (long)g_pairs.size();
        *species = g_species.data();
        *comp = g_comp.data();
        *conj = g_conj.data();
        *ref_used = f.r;
    });
}

}  // extern "C"

// Maximum entropy charge state deconvolution (written for MS Analysis from the
// published principle: Skilling and Bryan, Mon. Not. R. astr. Soc. 1984;
// Ferrige et al., Rapid Commun. Mass Spectrom. 1992; the Python reference of
// the 2.7 algorithm is ms_deconv.maxent_deconvolute, kept as the fallback).
//
// Design (2.8)
//   grid      the m/z axis of the fit is a uniform grid in ln(m/z) (TOF: the
//             peak width grows with m/z) or in m/z (resolution < 0), six
//             points per peak width, restricted to the windows the masses
//             reach at the chosen charges; the data are averaged into the
//             grid cells (area kept). The peak width of the linear grid is
//             the median over the tallest well separated peaks.
//   model     every mass gives a peak at each charge, weighted by a charge
//             envelope c(z | M): one sparse matrix (grid x masses, two entries
//             per mass and charge, CSR, 10 bytes per entry) and its transpose,
//             followed by a Gaussian blur along the grid (scipy
//             gaussian_filter1d, constant mode, truncate 4). The envelope is
//             defined at nodes across the mass range and interpolated linearly
//             in ln M between them (one node = one envelope shared by all
//             masses, the 2.7 model); each active mass stores its lower node
//             and the interpolation weight (6 bytes). The nodes start log
//             spaced (one per factor 1.5 in mass, at most 12); after the first
//             pass they are put at the ends of the mass range and at the
//             species found (local maxima of the smoothed mass intensity above
//             1 % of the largest, 2 % apart), so that every species carries
//             its own envelope.
//   noise     sigma_i^2 = s0^2 + ybar_i / g in the normalised units of the
//             gridded data: a noise floor and a counting (Poisson like) term,
//             both measured on the raw spectrum (variance of the second
//             differences against the mean in peak free windows, regressed as
//             a + b mean) and converted to the grid (the variance of a cell
//             average of a linear interpolant is 2/3 of a raw point's).
//             ybar is the larger of the three point running mean of the data,
//             the data seen through the peak shape, and the fit of the
//             previous round (Pearson's form), so that the flanks the model
//             puts under a peak are not weighted with the floor alone. After
//             the first pass the two terms are re-estimated from the
//             residuals in bins of intensity (the measured values are lower
//             bounds), which adds what the model cannot fit at the weak
//             peaks. Points the minimum intensity removed count as data
//             below the threshold, not as zero.
//   envelope  refitted after each optimisation by weighted non negative least
//             squares on the normal equations of the node columns (accelerated
//             projected gradient; Lawson Hanson for a single charge) with a
//             curvature penalty across the charges (15 or more), a difference
//             penalty between neighbouring nodes (a node without signal
//             follows its neighbours) and the minimum intensity ratio of
//             Ferrige et al.: a charge state's neighbours carry at least a
//             third of its intensity (weaker for z 1 and 2), so that no node
//             explains a stretch of the spectrum with one charge alone.
//   fit       u = ln f is optimised with an L-BFGS (two loop recursion, m = 10,
//             More-Thuente line search, projection onto the bounds, the
//             entropy metric 1 / f as the initial inverse Hessian) on
//             chi2 / 2 - alpha x entropy(f | flat default).
//   alpha     the classic criterion: alpha is searched (slope limited steps,
//             then a secant on the bracket) until chi squared per signal
//             point lies within 10 % of the target (1.0); when the best fit
//             stays above the target (the peak shape is not the model's) the
//             search ends at that best fit, when it falls below the target
//             the best fit is kept as well (chi squared = N would shrink the
//             strong peaks by a quarter on data the model explains to within
//             the noise). Bounded by the rounds parameter, at least
//             min_rounds; the progress callback counts the rounds.
//   output    the fitted peak of every charge state moved to the mass axis and
//             added (data units), then peak picking.
// OpenMP parallelism is over rows of the grid / masses; every grid point is
// summed by one thread in a fixed order, so results do not depend on the
// thread count.
#include "e2_internal.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <cstdint>
#include <memory>
#ifdef _OPENMP
#include <omp.h>
#endif

using ms::e2::Peak;

namespace {

// ============================================================ small helpers
inline int omp_threads() { return ms::threads(); }

template <class F> void parallel_range(size_t n, F&& f) {
    // f(i0, i1) over disjoint contiguous blocks
    int nt = std::max(1, omp_threads());
    if (n < 4096 || nt == 1) { f((size_t)0, n); return; }
#ifdef _OPENMP
#pragma omp parallel num_threads(nt)
    {
        int t = omp_get_thread_num(), T = omp_get_num_threads();
        size_t a = n * (size_t)t / T, b = n * (size_t)(t + 1) / T;
        if (b > a) f(a, b);
    }
#else
    f((size_t)0, n);
#endif
}

double dot(const double* a, const double* b, size_t n) {
    double s = 0;
    for (size_t i = 0; i < n; i++) s += a[i] * b[i];
    return s;
}

// np.interp(xq, xp, fp, left=0, right=0)
double interp0(double xq, const double* xp, const double* fp, size_t n) {
    if (n == 0 || xq < xp[0] || xq > xp[n - 1]) return 0.0;
    size_t j = std::upper_bound(xp, xp + n, xq) - xp;
    if (j == 0) return fp[0];
    if (j >= n) return fp[n - 1];
    j--;
    double dx = xp[j + 1] - xp[j];
    if (dx == 0) return fp[j];
    return (fp[j + 1] - fp[j]) / dx * (xq - xp[j]) + fp[j];
}

// np.interp with end values (clamped)
double interp_clamp(double xq, const double* xp, const double* fp, size_t n) {
    if (xq <= xp[0]) return fp[0];
    if (xq >= xp[n - 1]) return fp[n - 1];
    size_t j = std::upper_bound(xp, xp + n, xq) - xp - 1;
    double dx = xp[j + 1] - xp[j];
    if (dx == 0) return fp[j];
    return (fp[j + 1] - fp[j]) / dx * (xq - xp[j]) + fp[j];
}

// np.arange(start, stop, step) values
std::vector<double> arange(double start, double stop, double step) {
    double len = std::ceil((stop - start) / step);
    size_t n = len > 0 ? (size_t)len : 0;
    std::vector<double> v(n);
    for (size_t i = 0; i < n; i++) v[i] = start + (double)i * step;
    return v;
}

// binary dilation with the 3 point structure, `iters` iterations
void dilate(std::vector<char>& m, long iters) {
    size_t n = m.size();
    if (n == 0 || iters <= 0) return;
    std::vector<char> out(n, 0);
    long last = -1;   // last true index seen scanning left to right
    for (size_t i = 0; i < n; i++) {
        if (m[i]) last = (long)i;
        if (last >= 0 && (long)i - last <= iters) out[i] = 1;
    }
    long next = -1;
    for (long i = (long)n - 1; i >= 0; i--) {
        if (m[i]) next = i;
        if (next >= 0 && next - i <= iters) out[i] = 1;
    }
    m.swap(out);
}

// ============================================================ Gaussian blur
// scipy.ndimage.gaussian_filter1d(v, sigma, mode="constant", truncate=4.0)
struct Blur {
    int radius = 0;
    std::vector<double> k;   // k[t], t = 0..radius (symmetric)
    void init(double sigma, double truncate = 4.0) {
        radius = (int)(truncate * sigma + 0.5);
        std::vector<double> phi(2 * radius + 1);
        double s2 = sigma * sigma;
        for (int i = -radius; i <= radius; i++) phi[i + radius] = std::exp(-0.5 / s2 * (double)(i * i));
        double sum = ms::e2::np_sum(phi.data(), phi.size());
        k.assign(radius + 1, 0.0);
        for (int t = 0; t <= radius; t++) k[t] = phi[radius + t] / sum;
    }
    // out[i - i0] = blurred v at i for i in [i0, i1); v has n points, zero outside
    void apply_part(const double* v, size_t n, size_t i0, size_t i1, double* out) const {
        int r = radius;
        for (size_t i = i0; i < i1; i++) {
            double s = v[i] * k[0];
            if (i >= (size_t)r && i + r < n) {
                for (int t = r; t >= 1; t--) s += (v[i - t] + v[i + t]) * k[t];
            } else {
                for (int t = r; t >= 1; t--) {
                    double l = (long)i - t >= 0 ? v[i - t] : 0.0;
                    double h = i + t < n ? v[i + t] : 0.0;
                    s += (l + h) * k[t];
                }
            }
            out[i - i0] = s;
        }
    }
    void apply(const double* v, double* out, size_t n) const {
        parallel_range(n, [&](size_t a, size_t b) { apply_part(v, n, a, b, out + a); });
    }
};

// ============================================================ charge model
// A (grid x active masses): two entries per mass and charge (linear
// interpolation between the two grid points around (m + z a)/z), stored as
// CSR by row and as its transpose by column. Per entry: column (int32),
// charge index (uint16), weight (float32), 10 bytes; the envelope weight
// c(z | M) is applied on the fly from the node table, so a new envelope
// costs nothing. Per active mass: the lower envelope node (uint16) and the
// weight of the upper one (float32).
struct Model {
    size_t n = 0;        // grid points
    int K = 0;           // charges
    size_t na = 0;       // active masses
    std::vector<int> act;        // active mass -> mass index
    double sig = 1;
    Blur blur;
    // A (grid x active): rows in order, within a row by charge then mass
    std::vector<size_t> a_ptr; std::vector<int> a_col; std::vector<uint16_t> a_k; std::vector<float> a_w;
    // transpose (active x grid)
    std::vector<size_t> t_ptr; std::vector<int> t_col; std::vector<uint16_t> t_k; std::vector<float> t_w;
    // envelope nodes: Q nodes log spaced over the mass range, C is (Q + 1) x K
    // (row Q repeats row Q - 1 so that the interpolation never reads past the end)
    int Q = 1;
    double lnm0 = 0, lnm1 = 0;        // ln of the first and last mass
    std::vector<double> node_mass;    // Q
    std::vector<double> coord_cp;     // the grid coordinate (copy), for the column supports
    bool uselog = false;
    double adduct = 1.0;
    std::vector<int> charges_cp;
    std::vector<uint16_t> nq;         // active mass -> lower node
    std::vector<float> nt;            // weight of the upper node (0 when Q == 1)
    std::vector<double> C;
    std::vector<double> tmp;          // scratch of length n

    // bytes the model takes for nnz entries
    static double bytes_for(size_t nnz, size_t n, size_t na) {
        return 2.0 * nnz * (sizeof(int) + sizeof(uint16_t) + sizeof(float)) + 8.0 * (n + na + 2) + 8.0 * n + 6.0 * na;
    }

    inline double env(size_t j, int k) const {
        const double* c0 = C.data() + (size_t)nq[j] * K;
        return c0[k] + nt[j] * (c0[K + k] - c0[k]);
    }
    // node position of a mass: lower node and weight of the upper one (linear in ln M between the nodes)
    void node_of(double mass, int& q, double& t) const {
        if (Q <= 1) { q = 0; t = 0; return; }
        if (mass <= node_mass.front()) { q = 0; t = 0; return; }
        if (mass >= node_mass.back()) { q = Q - 2; t = 1; return; }
        q = (int)(std::upper_bound(node_mass.begin(), node_mass.end(), mass) - node_mass.begin()) - 1;
        q = std::min(std::max(q, 0), Q - 2);
        double l0 = std::log(node_mass[q]), l1 = std::log(node_mass[q + 1]);
        t = l1 > l0 ? (std::log(mass) - l0) / (l1 - l0) : 0.0;
        t = std::min(1.0, std::max(0.0, t));
    }
    // new nodes (sorted masses, the first and last at the ends of the mass range): every
    // active mass is reassigned and the envelope table is read off the old nodes
    void set_nodes(const std::vector<double>& nodes, const std::vector<double>& masses) {
        std::vector<double> cnew((size_t)nodes.size() * K);
        for (size_t q = 0; q < nodes.size(); q++)
            for (int k = 0; k < K; k++) cnew[q * K + k] = env_at(nodes[q], k);
        Q = (int)nodes.size();
        node_mass = nodes;
        for (size_t j = 0; j < na; j++) {
            int q; double t;
            node_of(masses[act[j]], q, t);
            nq[j] = (uint16_t)q; nt[j] = (float)t;
        }
        set_envelope(cnew);
    }
    // envelope at any mass (for the output)
    double env_at(double mass, int k) const {
        int q; double t;
        node_of(mass, q, t);
        const double* c0 = C.data() + (size_t)q * K;
        return c0[k] + t * (c0[K + k] - c0[k]);
    }

    void build(const std::vector<double>& coord, double d, const std::vector<double>& masses,
               const std::vector<int>& charges, double a, double fwhm, bool log, double budget, int nodes) {
        n = coord.size();
        K = (int)charges.size();
        if (K > 65535) throw std::invalid_argument("too many charge states");
        sig = std::max(0.3, fwhm / d / 2.3548);
        blur.init(sig);
        coord_cp = coord; uselog = log; adduct = a; charges_cp = charges;
        const size_t nm = masses.size();
        std::vector<char> active(nm, 0);
        std::vector<size_t> row_cnt(n + 1, 0), col_cnt(nm, 0);
        // pass 1: count entries per row and per mass (positions increase with the mass)
        auto walk = [&](int k, auto&& fn) {
            const double z = charges[k];
            long i0 = -1;
            for (size_t j = 0; j < nm; j++) {
                double pos = (masses[j] + z * a) / z;
                if (log) pos = std::log(std::max(pos, 1e-9));
                while (i0 + 1 < (long)n && coord[i0 + 1] <= pos) i0++;
                if (i0 < 0 || i0 >= (long)n - 1) continue;
                if (!(coord[i0 + 1] - coord[i0] < 1.5 * d)) continue;   // across a gap
                const double w1 = std::min(1.0, std::max(0.0, (pos - coord[i0]) / d));
                fn(j, (size_t)i0, w1);
            }
        };
        size_t nnz = 0;
        for (int k = 0; k < K; k++)
            walk(k, [&](size_t j, size_t i0, double) { active[j] = 1; row_cnt[i0]++; row_cnt[i0 + 1]++; col_cnt[j] += 2; nnz += 2; });
        std::vector<int> colmap(nm, -1);
        act.clear();
        for (size_t j = 0; j < nm; j++) if (active[j]) { colmap[j] = (int)act.size(); act.push_back((int)j); }
        na = act.size();
        {
            const double need = bytes_for(nnz, n, na) + 230.0 * na + 16.0 * K * n;
            if (need > budget) {
                char b[600];
                snprintf(b, sizeof b, "This deconvolution needs about %s of memory (%zu masses x %d charges on %zu grid points); this computer can give it %s. Use a larger mass step or narrower mass and charge ranges (isotope peaks need 0.01 to 0.02 Da, envelopes 0.5 to 1 Da).",
                         ms::fmt_bytes(need).c_str(), na, K, n, ms::fmt_bytes(budget).c_str());
                throw std::invalid_argument(b);
            }
        }
        a_ptr.assign(n + 1, 0);
        for (size_t i = 0; i < n; i++) a_ptr[i + 1] = a_ptr[i] + row_cnt[i];
        const size_t nac = std::max<size_t>(na, 1);
        t_ptr.assign(nac + 1, 0);
        for (size_t j = 0, q = 0; j < nm; j++) if (active[j]) { t_ptr[q + 1] = t_ptr[q] + col_cnt[j]; q++; }
        a_col.resize(nnz); a_k.resize(nnz); a_w.resize(nnz);
        t_col.resize(nnz); t_k.resize(nnz); t_w.resize(nnz);
        std::vector<size_t> a_fill(a_ptr.begin(), a_ptr.end() - 1), t_fill(t_ptr.begin(), t_ptr.end() - 1);
        // pass 2: fill (rows: by charge then mass; columns: by charge, then grid point)
        for (int k = 0; k < K; k++)
            walk(k, [&](size_t j, size_t i0, double w1) {
                const int cj = colmap[j];
                size_t p = a_fill[i0]++;
                a_col[p] = cj; a_k[p] = (uint16_t)k; a_w[p] = (float)(1.0 - w1);
                p = a_fill[i0 + 1]++;
                a_col[p] = cj; a_k[p] = (uint16_t)k; a_w[p] = (float)w1;
                p = t_fill[cj]++;
                t_col[p] = (int)i0; t_k[p] = (uint16_t)k; t_w[p] = (float)(1.0 - w1);
                p = t_fill[cj]++;
                t_col[p] = (int)i0 + 1; t_k[p] = (uint16_t)k; t_w[p] = (float)w1;
            });
        tmp.assign(n, 0.0);
        // envelope nodes, log spaced to start with
        Q = std::max(1, nodes);
        lnm0 = std::log(std::max(masses.front(), 1e-9)); lnm1 = std::log(std::max(masses.back(), 1e-9));
        if (!(lnm1 > lnm0)) Q = 1;
        node_mass.resize(Q);
        for (int q = 0; q < Q; q++) node_mass[q] = Q > 1 ? std::exp(lnm0 + (lnm1 - lnm0) * q / (Q - 1)) : 0.5 * (masses.front() + masses.back());
        nq.assign(nac, 0); nt.assign(nac, 0.0f);
        for (size_t j = 0; j < na; j++) {
            int q; double t;
            node_of(masses[act[j]], q, t);
            nq[j] = (uint16_t)q; nt[j] = (float)t;
        }
        std::vector<double> c0((size_t)Q * K, 1.0 / std::max(K, 1));
        set_envelope(c0);
    }

    // C (Q x K) -> the padded table
    void set_envelope(const std::vector<double>& cqk) {
        C.assign((size_t)(Q + 1) * K, 0.0);
        std::copy(cqk.begin(), cqk.begin() + (size_t)Q * K, C.begin());
        std::copy(C.begin() + (size_t)(Q - 1) * K, C.begin() + (size_t)Q * K, C.begin() + (size_t)Q * K);
    }
    std::vector<double> envelope() const { return std::vector<double>(C.begin(), C.begin() + (size_t)Q * K); }

    // m/z spectrum (grid) of the masses f with the envelope
    void forward(const double* f, double* out) {
        parallel_range(n, [&](size_t a, size_t b) {
            for (size_t i = a; i < b; i++) {
                double s = 0;
                for (size_t p = a_ptr[i]; p < a_ptr[i + 1]; p++) s += env(a_col[p], a_k[p]) * a_w[p] * f[a_col[p]];
                tmp[i] = s;
            }
        });
        blur.apply(tmp.data(), out, n);
    }

    // gradient with respect to the masses of a (weighted) residual r on the grid
    void adjoint_f(const double* r, double* out) {
        blur.apply(r, tmp.data(), n);
        const double* g = tmp.data();
        parallel_range(na, [&](size_t a, size_t b) {
            for (size_t j = a; j < b; j++) {
                double s = 0;
                for (size_t p = t_ptr[j]; p < t_ptr[j + 1]; p++) s += env(j, t_k[p]) * t_w[p] * g[t_col[p]];
                out[j] = s;
            }
        });
    }

    // unweighted m/z spectrum of each charge state: K x n (row major), for the output
    void columns(const double* f, std::vector<double>& out) {
        out.assign((size_t)K * n, 0.0);
        std::vector<double> raw((size_t)K * n, 0.0);
        parallel_range(n, [&](size_t a, size_t b) {
            for (size_t i = a; i < b; i++)
                for (size_t p = a_ptr[i]; p < a_ptr[i + 1]; p++)
                    raw[(size_t)a_k[p] * n + i] += a_w[p] * f[a_col[p]];
        });
        for (int k = 0; k < K; k++) blur.apply(raw.data() + (size_t)k * n, out.data() + (size_t)k * n, n);
    }

    // Normal equations of the envelope fit: the column of node q and charge k
    // is col_qk = blur(sum over the entries of charge k of the node weight x
    // A x f); G[(q,k),(q',k')] = sum_i w_i col_qk(i) col_q'k'(i), b = sum_i
    // w_i col_qk(i) y_i. Computed in blocks of the grid (with the halo the
    // blur reaches), so the Q x K columns are never held whole; every sum
    // runs in a fixed order (results independent of the thread count).
    // grid index range the column of node q and charge k can be non zero in: the
    // masses between the neighbouring nodes at that charge, widened by the blur
    std::pair<long, long> support(int q, int k) const {
        const double m0 = node_mass[q > 0 ? q - 1 : 0], m1 = node_mass[q + 1 < Q ? q + 1 : Q - 1];
        const double z = charges_cp[k];
        double c0 = (m0 + z * adduct) / z, c1 = (m1 + z * adduct) / z;
        if (uselog) { c0 = std::log(std::max(c0, 1e-9)); c1 = std::log(std::max(c1, 1e-9)); }
        long i0 = std::lower_bound(coord_cp.begin(), coord_cp.end(), c0) - coord_cp.begin();
        long i1 = std::upper_bound(coord_cp.begin(), coord_cp.end(), c1) - coord_cp.begin();
        i0 = std::max(0L, i0 - blur.radius - 2); i1 = std::min((long)n, i1 + blur.radius + 2);
        return {i0, i1};
    }

    void node_normal(const double* f, const double* y, const double* w, std::vector<double>& G, std::vector<double>& b) const {
        const int P = Q * K;
        G.assign((size_t)P * P, 0.0); b.assign(P, 0.0);
        const int r = blur.radius;
        const size_t B = std::max<size_t>(1024, std::min<size_t>(32768, (size_t)(1.5e6 / std::max(P, 1))));   // block buffers of about 24 MB
        std::vector<std::pair<long, long>> sup(P);
        for (int q = 0; q < Q; q++) for (int k = 0; k < K; k++) sup[(size_t)q * K + k] = Q > 1 ? support(q, k) : std::pair<long, long>{0L, (long)n};
        // only pairs of columns whose supports overlap
        std::vector<std::pair<int, int>> pairs;
        for (int p = 0; p < P; p++) for (int q = p; q < P; q++)
            if (std::min(sup[p].second, sup[q].second) > std::max(sup[p].first, sup[q].first)) pairs.push_back({p, q});
        std::vector<double> gv(pairs.size(), 0.0), raw, col;
        for (size_t a = 0; a < n; a += B) {
            const size_t b1 = std::min(n, a + B);
            const size_t ea = a >= (size_t)r ? a - r : 0, eb = std::min(n, b1 + r);
            const size_t L = eb - ea, M = b1 - a;
            raw.assign((size_t)P * L, 0.0);
            parallel_range(L, [&](size_t i0, size_t i1) {
                for (size_t ii = i0; ii < i1; ii++) {
                    const size_t i = ea + ii;
                    for (size_t p = a_ptr[i]; p < a_ptr[i + 1]; p++) {
                        const int j = a_col[p], k = a_k[p];
                        const double v = a_w[p] * f[j], t = nt[j];
                        const size_t q = nq[j];
                        raw[(q * K + k) * L + ii] += v * (1.0 - t);
                        if (t > 0) raw[((q + 1) * K + k) * L + ii] += v * t;
                    }
                }
            });
            col.assign((size_t)P * M, 0.0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(omp_threads())
#endif
            for (int p = 0; p < P; p++)
                if (sup[p].second > (long)a && sup[p].first < (long)b1)
                    blur.apply_part(raw.data() + (size_t)p * L, L, a - ea, b1 - ea, col.data() + (size_t)p * M);
            const double* ya = y + a;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(omp_threads())
#endif
            for (long t = 0; t < (long)pairs.size(); t++) {
                const int p = pairs[t].first, q = pairs[t].second;
                const long lo = std::max({(long)a, sup[p].first, sup[q].first}), hi = std::min({(long)b1, sup[p].second, sup[q].second});
                if (hi <= lo) continue;
                const double* cp = col.data() + (size_t)p * M - a;
                const double* cq = col.data() + (size_t)q * M - a;
                double s = 0;
                for (long i = lo; i < hi; i++) s += w[i] * cp[i] * cq[i];
                gv[t] += s;
            }
            for (int p = 0; p < P; p++) {
                const long lo = std::max((long)a, sup[p].first), hi = std::min((long)b1, sup[p].second);
                if (hi <= lo) continue;
                const double* cp = col.data() + (size_t)p * M - a;
                double s = 0;
                for (long i = lo; i < hi; i++) s += w[i] * cp[i] * ya[i - a];
                b[p] += s;
            }
        }
        for (size_t t = 0; t < pairs.size(); t++) {
            G[(size_t)pairs[t].first * P + pairs[t].second] = gv[t];
            G[(size_t)pairs[t].second * P + pairs[t].first] = gv[t];
        }
    }
};

// ============================================================ NNLS (Lawson Hanson)
// argmin |A x - b|, x >= 0; A is m x n row major
std::vector<double> nnls(const std::vector<double>& A, const std::vector<double>& b, int m, int n) {
    std::vector<double> x(n, 0.0), w(n), s(n);
    std::vector<char> P(n, 0);
    int maxiter = 3 * n;
    auto residual_grad = [&]() {   // w = A^T (b - A x)
        std::vector<double> r(m);
        for (int i = 0; i < m; i++) {
            double t = b[i];
            for (int j = 0; j < n; j++) t -= A[i * n + j] * x[j];
            r[i] = t;
        }
        for (int j = 0; j < n; j++) {
            double t = 0;
            for (int i = 0; i < m; i++) t += A[i * n + j] * r[i];
            w[j] = t;
        }
    };
    auto solve_P = [&]() {         // least squares on the passive columns
        std::vector<int> idx;
        for (int j = 0; j < n; j++) if (P[j]) idx.push_back(j);
        int np = (int)idx.size();
        std::fill(s.begin(), s.end(), 0.0);
        if (np == 0) return;
        std::vector<double> Ap((size_t)m * np);
        for (int i = 0; i < m; i++) for (int q = 0; q < np; q++) Ap[(size_t)i * np + q] = A[i * n + idx[q]];
        std::vector<double> sp;
        ms::e2::lstsq_small(Ap, b, m, np, sp);
        for (int q = 0; q < np; q++) s[idx[q]] = sp[q];
    };
    int iter = 0;
    double amax = 0;
    for (double v : A) amax = std::max(amax, std::fabs(v));
    double tol = 10 * std::numeric_limits<double>::epsilon() * amax * std::max(m, n);
    while (true) {
        residual_grad();
        int jbest = -1; double wbest = 0;
        for (int j = 0; j < n; j++) if (!P[j] && w[j] > wbest) { wbest = w[j]; jbest = j; }
        if (jbest < 0 || wbest <= tol) break;
        P[jbest] = 1;
        while (true) {
            if (++iter > maxiter) return x;
            solve_P();
            bool feasible = true;
            for (int j = 0; j < n; j++) if (P[j] && s[j] <= 0) { feasible = false; break; }
            if (feasible) { x = s; break; }
            double alpha = std::numeric_limits<double>::infinity();
            int jalpha = -1;
            for (int j = 0; j < n; j++) if (P[j] && s[j] <= 0) {
                double den = x[j] - s[j];
                double a = den != 0 ? x[j] / den : 0.0;
                if (a < alpha) { alpha = a; jalpha = j; }
            }
            for (int j = 0; j < n; j++) x[j] += alpha * (s[j] - x[j]);
            if (jalpha >= 0) { x[jalpha] = 0.0; P[jalpha] = 0; }
            for (int j = 0; j < n; j++) if (P[j] && x[j] <= 0) { P[j] = 0; x[j] = 0.0; }
            bool any = false; for (int j = 0; j < n; j++) any |= P[j];
            if (!any) break;
        }
    }
    return x;
}

// argmin 1/2 x^T G x - b^T x + R(x), x >= 0 (G symmetric positive semidefinite)
// by an accelerated projected gradient with restarts, warm started from x0.
// R is the one sided ratio penalty of the charge envelope (Ferrige et al.:
// a charge state's neighbours must carry at least `ratio` of its intensity;
// a species never shows one charge state alone): for every node and
// neighbouring charges z, z + 1, lam_z (max(0, r x_z - x_z+1)^2 +
// max(0, r x_z+1 - x_z)^2) with lam_z = pen[z] of that node. Also used for
// the plain problem (pen empty) when the Lawson Hanson solver would be slow.
std::vector<double> nnls_qp(const std::vector<double>& G, const std::vector<double>& b, int n, std::vector<double> x0,
                            int Q = 1, int K = 1, double ratio = 0.0, const std::vector<double>& pen = std::vector<double>()) {
    std::vector<double> x(n, 0.0), y(n), g(n), xp(n);
    const bool use_pen = ratio > 0 && !pen.empty() && Q * K == n;
    if ((int)x0.size() == n) {
        // scale the warm start to the least squares optimum along itself
        double num = 0, den = 0;
        for (int i = 0; i < n; i++) { double gi = 0; for (int j = 0; j < n; j++) gi += G[(size_t)i * n + j] * x0[j]; num += b[i] * x0[i]; den += gi * x0[i]; }
        double sc = den > 0 && num > 0 ? num / den : 0.0;
        for (int i = 0; i < n; i++) x[i] = std::max(0.0, x0[i] * sc);
    }
    // largest eigenvalue (power iteration) for the step
    std::vector<double> v(n, 1.0), Gv(n);
    double lam = 0;
    for (int it = 0; it < 40; it++) {
        double nv = 0; for (double t : v) nv += t * t; nv = std::sqrt(nv);
        if (nv <= 0) break;
        for (auto& t : v) t /= nv;
        for (int i = 0; i < n; i++) { double sm = 0; for (int j = 0; j < n; j++) sm += G[(size_t)i * n + j] * v[j]; Gv[i] = sm; }
        lam = 0; for (int i = 0; i < n; i++) lam += v[i] * Gv[i];
        v = Gv;
    }
    if (!(lam > 0)) return x;
    double lpen = 0;
    if (use_pen) for (double t : pen) lpen = std::max(lpen, t);
    const double step = 1.0 / (1.05 * lam + 4.0 * (1.0 + ratio * ratio) * lpen);
    auto add_pen_grad = [&](const std::vector<double>& z, std::vector<double>& grad) {
        if (!use_pen) return;
        for (int q = 0; q < Q; q++)
            for (int k = 0; k + 1 < K; k++) {
                const size_t i0 = (size_t)q * K + k, i1 = i0 + 1;
                const double l = pen[i0];
                if (!(l > 0)) continue;
                double d1 = std::max(0.0, ratio * z[i0] - z[i1]);   // the upper neighbour too small
                double d2 = std::max(0.0, ratio * z[i1] - z[i0]);   // the lower neighbour too small
                grad[i0] += l * (2.0 * ratio * d1 - 2.0 * d2);
                grad[i1] += l * (-2.0 * d1 + 2.0 * ratio * d2);
            }
    };
    y = x; xp = x;
    double tk = 1.0;
    for (int it = 0; it < 3000; it++) {
        for (int i = 0; i < n; i++) { double sm = -b[i]; for (int j = 0; j < n; j++) sm += G[(size_t)i * n + j] * y[j]; g[i] = sm; }
        add_pen_grad(y, g);
        double change = 0, norm = 0, gy = 0;
        for (int i = 0; i < n; i++) {
            double xn = std::max(0.0, y[i] - step * g[i]);
            gy += (xn - xp[i]) * (y[i] - xn);
            change += (xn - xp[i]) * (xn - xp[i]); norm += xn * xn;
            x[i] = xn;
        }
        if (gy > 0) {   // momentum points uphill: restart
            tk = 1.0; y = x; xp = x; continue;
        }
        double tn = 0.5 * (1.0 + std::sqrt(1.0 + 4.0 * tk * tk));
        double beta = (tk - 1.0) / tn;
        for (int i = 0; i < n; i++) y[i] = x[i] + beta * (x[i] - xp[i]);
        xp = x; tk = tn;
        if (change <= 1e-12 * std::max(norm, 1e-300)) break;
    }
    return x;
}

// Charge envelope at the nodes (Q x K, row major) from the normal equations:
// weighted NNLS with
//  - a curvature penalty across the charges (15 or more: ESI envelopes of
//    proteins are smooth, and an envelope spiking at a few charges lets
//    masses at M x z'/z (ghosts) explain the data),
//  - a difference penalty between neighbouring nodes (a node whose masses
//    carry no signal follows its neighbours instead of going wild),
//  - the minimum intensity ratio between neighbouring charge states (Ferrige
//    et al. 1992): a species never shows a single charge state, so a node
//    cannot explain a stretch of the spectrum with one charge alone. Weighted
//    by the data curvature of the two charges, which makes the spike dearer
//    than any fit it could buy.
std::vector<double> envelope_fit(std::vector<double> G, const std::vector<double>& b, int Q, int K, int zfirst,
                                 const std::vector<double>& warm, int smooth = 15, double lam_node = 0.1,
                                 double ratio = 0.33, double lam_ratio = 3.0) {
    const int P = Q * K;
    double tr = 0;
    for (int p = 0; p < P; p++) tr += G[(size_t)p * P + p];
    const double scale = tr / std::max(P, 1);
    std::vector<double> pen;
    if (ratio > 0 && K >= 2) {
        // weaker between the lowest charges: a small molecule is often seen singly charged alone
        pen.assign(P, 0.0);
        for (int q = 0; q < Q; q++)
            for (int k = 0; k + 1 < K; k++) {
                const size_t i0 = (size_t)q * K + k, i1 = i0 + 1;
                const int z = zfirst + k;
                const double strength = z >= 3 ? 1.0 : z == 2 ? 0.5 : 0.25;
                pen[i0] = strength * lam_ratio * std::max({G[i0 * P + i0], G[i1 * P + i1], 1e-3 * scale});
            }
    }
    if (K >= smooth && K > 2) {
        // d2^T d2 with rows (1, -2, 1) within every node
        for (int q = 0; q < Q; q++)
            for (int i = 0; i < K - 2; i++) {
                double row[3] = {1.0, -2.0, 1.0};
                for (int p = 0; p < 3; p++) for (int t = 0; t < 3; t++)
                    G[(size_t)(q * K + i + p) * P + (q * K + i + t)] += 0.1 * scale * row[p] * row[t];
            }
    }
    if (Q > 1 && lam_node > 0) {
        for (int q = 0; q + 1 < Q; q++)
            for (int k = 0; k < K; k++) {
                const size_t i0 = (size_t)q * K + k, i1 = i0 + K;
                G[i0 * P + i0] += lam_node * scale; G[i1 * P + i1] += lam_node * scale;
                G[i0 * P + i1] -= lam_node * scale; G[i1 * P + i0] -= lam_node * scale;
            }
    }
    double dmax = 0;
    for (int p = 0; p < P; p++) dmax = std::max(dmax, G[(size_t)p * P + p]);
    for (int p = 0; p < P; p++) G[(size_t)p * P + p] += 1e-12 * std::max(dmax, 1e-300);
    if (P > 60 || !pen.empty()) return nnls_qp(G, b, P, warm, Q, K, ratio, pen);
    // Cholesky G = L L^T (lower); a larger ridge when it fails
    std::vector<double> L((size_t)P * P, 0.0);
    for (int attempt = 0; attempt < 6; attempt++) {
        bool ok = true;
        std::fill(L.begin(), L.end(), 0.0);
        for (int i = 0; i < P && ok; i++) {
            for (int j = 0; j <= i; j++) {
                double sm = G[(size_t)i * P + j];
                for (int p = 0; p < j; p++) sm -= L[(size_t)i * P + p] * L[(size_t)j * P + p];
                if (i == j) {
                    if (!(sm > 0)) { ok = false; break; }
                    L[(size_t)i * P + i] = std::sqrt(sm);
                } else {
                    L[(size_t)i * P + j] = sm / L[(size_t)j * P + j];
                }
            }
        }
        if (ok) break;
        for (int p = 0; p < P; p++) G[(size_t)p * P + p] += std::pow(10.0, attempt - 10) * std::max(dmax, 1e-300);
    }
    // z = L^-1 b
    std::vector<double> z(P);
    for (int i = 0; i < P; i++) {
        double sm = b[i];
        for (int p = 0; p < i; p++) sm -= L[(size_t)i * P + p] * z[p];
        z[i] = L[(size_t)i * P + i] != 0 ? sm / L[(size_t)i * P + i] : 0.0;
    }
    // nnls(L^T, z)
    std::vector<double> LT((size_t)P * P);
    for (int i = 0; i < P; i++) for (int j = 0; j < P; j++) LT[(size_t)i * P + j] = L[(size_t)j * P + i];
    return nnls(LT, z, P, P);
}

// ============================================================ objective
// chi2 / 2 - alpha x entropy, in u = ln f; w are the weights 1 / sigma_i^2
struct Objective {
    Model* model; const std::vector<double>* y; const std::vector<double>* w; double alpha, dflt, log_dflt;
    std::vector<double> f, r, adj, fit;
    size_t evals = 0;
    double operator()(const std::vector<double>& u, std::vector<double>& grad) {
        size_t na = u.size(), n = y->size();
        f.resize(na); r.resize(n); adj.resize(na); fit.resize(n);
        for (size_t j = 0; j < na; j++) f[j] = std::exp(u[j]);
        model->forward(f.data(), fit.data());
        const double* yy = y->data(); const double* ww = w->data();
        double chi2 = 0;
        for (size_t i = 0; i < n; i++) { double d = fit[i] - yy[i]; chi2 += ww[i] * d * d; r[i] = ww[i] * d; }
        chi2 *= 0.5;
        double ent = 0;
        for (size_t j = 0; j < na; j++) ent += f[j] - dflt - f[j] * (u[j] - log_dflt);
        model->adjoint_f(r.data(), adj.data());
        grad.resize(na);
        for (size_t j = 0; j < na; j++) grad[j] = f[j] * (adj[j] + alpha * (u[j] - log_dflt));
        evals++;
        return chi2 - alpha * ent;
    }
};

// ---- More-Thuente line search (MINPACK-2 dcsrch / dcstep, as L-BFGS-B uses it)
struct Dcsrch {
    double ftol = 1e-3, gtol = 0.9, xtol = 0.1, stpmin = 0.0, stpmax = 1e10;
    bool brackt = false; int stage = 1;
    double finit = 0, ginit = 0, gtest = 0, width = 0, width1 = 0;
    double stx = 0, fx = 0, gx = 0, sty = 0, fy = 0, gy = 0, stmin = 0, stmax = 0;
    enum Task { FG, CONV, WARN, ERR };

    static double sgn(double v) { return v > 0 ? 1.0 : (v < 0 ? -1.0 : 0.0); }

    static void dcstep(double& stx, double& fx, double& dx, double& sty, double& fy, double& dy, double& stp,
                       double fp, double dp, bool& brackt, double stpmin, double stpmax) {
        double sgnd = sgn(dp) * sgn(dx);
        double stpf;
        if (fp > fx) {
            double theta = 3.0 * (fx - fp) / (stp - stx) + dx + dp;
            double s = std::max({std::fabs(theta), std::fabs(dx), std::fabs(dp)});
            double gamma = s * std::sqrt((theta / s) * (theta / s) - (dx / s) * (dp / s));
            if (stp < stx) gamma = -gamma;
            double p = (gamma - dx) + theta;
            double q = ((gamma - dx) + gamma) + dp;
            double r = p / q;
            double stpc = stx + r * (stp - stx);
            double stpq = stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / 2.0) * (stp - stx);
            if (std::fabs(stpc - stx) <= std::fabs(stpq - stx)) stpf = stpc;
            else stpf = stpc + (stpq - stpc) / 2.0;
            brackt = true;
        } else if (sgnd < 0.0) {
            double theta = 3.0 * (fx - fp) / (stp - stx) + dx + dp;
            double s = std::max({std::fabs(theta), std::fabs(dx), std::fabs(dp)});
            double gamma = s * std::sqrt((theta / s) * (theta / s) - (dx / s) * (dp / s));
            if (stp > stx) gamma = -gamma;
            double p = (gamma - dp) + theta;
            double q = ((gamma - dp) + gamma) + dx;
            double r = p / q;
            double stpc = stp + r * (stx - stp);
            double stpq = stp + (dp / (dp - dx)) * (stx - stp);
            if (std::fabs(stpc - stp) > std::fabs(stpq - stp)) stpf = stpc;
            else stpf = stpq;
            brackt = true;
        } else if (std::fabs(dp) < std::fabs(dx)) {
            double theta = 3.0 * (fx - fp) / (stp - stx) + dx + dp;
            double s = std::max({std::fabs(theta), std::fabs(dx), std::fabs(dp)});
            double gamma = s * std::sqrt(std::max(0.0, (theta / s) * (theta / s) - (dx / s) * (dp / s)));
            if (stp > stx) gamma = -gamma;
            double p = (gamma - dp) + theta;
            double q = (gamma + (dx - dp)) + gamma;
            double r = p / q;
            double stpc;
            if (r < 0 && gamma != 0) stpc = stp + r * (stx - stp);
            else if (stp > stx) stpc = stpmax;
            else stpc = stpmin;
            double stpq = stp + (dp / (dp - dx)) * (stx - stp);
            if (brackt) {
                if (std::fabs(stpc - stp) < std::fabs(stpq - stp)) stpf = stpc;
                else stpf = stpq;
                if (stp > stx) stpf = std::min(stp + 0.66 * (sty - stp), stpf);
                else stpf = std::max(stp + 0.66 * (sty - stp), stpf);
            } else {
                if (std::fabs(stpc - stp) > std::fabs(stpq - stp)) stpf = stpc;
                else stpf = stpq;
                stpf = std::min(stpmax, std::max(stpmin, stpf));
            }
        } else {
            if (brackt) {
                double theta = 3.0 * (fp - fy) / (sty - stp) + dy + dp;
                double s = std::max({std::fabs(theta), std::fabs(dy), std::fabs(dp)});
                double gamma = s * std::sqrt((theta / s) * (theta / s) - (dy / s) * (dp / s));
                if (stp > sty) gamma = -gamma;
                double p = (gamma - dp) + theta;
                double q = ((gamma - dp) + gamma) + dy;
                double r = p / q;
                stpf = stp + r * (sty - stp);
            } else if (stp > stx) stpf = stpmax;
            else stpf = stpmin;
        }
        if (fp > fx) { sty = stp; fy = fp; dy = dp; }
        else {
            if (sgnd < 0) { sty = stx; fy = fx; dy = dx; }
            stx = stp; fx = fp; dx = dp;
        }
        stp = stpf;
    }

    Task start(double& stp, double f, double g) {
        if (stp < stpmin || stp > stpmax || g >= 0) return ERR;
        brackt = false; stage = 1; finit = f; ginit = g; gtest = ftol * ginit;
        width = stpmax - stpmin; width1 = width / 0.5;
        stx = 0; fx = finit; gx = ginit; sty = 0; fy = finit; gy = ginit;
        stmin = 0; stmax = stp + 4.0 * stp;
        return FG;
    }

    // f, g: function and directional derivative at stp; on FG evaluate at the new stp
    Task iterate(double& stp, double f, double g) {
        const double xtrapl = 1.1, xtrapu = 4.0;
        double ftest = finit + stp * gtest;
        if (stage == 1 && f <= ftest && g >= 0) stage = 2;
        Task task = FG;
        if (brackt && (stp <= stmin || stp >= stmax)) task = WARN;
        if (brackt && stmax - stmin <= xtol * stmax) task = WARN;
        if (stp == stpmax && f <= ftest && g <= gtest) task = WARN;
        if (stp == stpmin && (f > ftest || g >= gtest)) task = WARN;
        if (f <= ftest && std::fabs(g) <= gtol * (-ginit)) task = CONV;
        if (task != FG) return task;
        if (stage == 1 && f <= fx && f > ftest) {
            double fm = f - stp * gtest, fxm = fx - stx * gtest, fym = fy - sty * gtest;
            double gm = g - gtest, gxm = gx - gtest, gym = gy - gtest;
            dcstep(stx, fxm, gxm, sty, fym, gym, stp, fm, gm, brackt, stmin, stmax);
            fx = fxm + stx * gtest; fy = fym + sty * gtest; gx = gxm + gtest; gy = gym + gtest;
        } else {
            dcstep(stx, fx, gx, sty, fy, gy, stp, f, g, brackt, stmin, stmax);
        }
        if (brackt) {
            if (std::fabs(sty - stx) >= 0.66 * width1) stp = stx + 0.5 * (sty - stx);
            width1 = width; width = std::fabs(sty - stx);
        }
        if (brackt) { stmin = std::min(stx, sty); stmax = std::max(stx, sty); }
        else { stmin = stp + xtrapl * (stp - stx); stmax = stp + xtrapu * (stp - stx); }
        if (!(stp >= stpmin)) stp = stpmin;     // also catches NaN
        if (stp > stpmax) stp = stpmax;
        if ((brackt && (stp <= stmin || stp >= stmax)) || (brackt && stmax - stmin <= xtol * stmax)) stp = stx;
        return FG;
    }
};

// L-BFGS with bounds: two loop recursion (m pairs), projected search direction,
// More-Thuente line search (ftol 1e-3, gtol 0.9, xtol 0.1, 20 evaluations at
// most, as L-BFGS-B). The initial inverse Hessian is the diagonal
// gamma / f, with f taken at the start of the call (the entropy metric of
// Skilling and Bryan: in u = ln f a step proportional to dQ/df, so that a
// mass still at the default model grows as fast as a strong one; a plain
// L-BFGS in u moves it by f x dQ/df and needs a hundred iterations to raise
// it from the default). The first direction is that scaled gradient with a
// unit step for the strongest pull. Stops on maxiter, the projected gradient
// (pgtol) or the relative decrease (ftol).
void lbfgs_b(Objective& obj, std::vector<double>& x, double lo, double hi, int maxiter, int m = 10,
             double ftol = 2.220446049250313e-09, double pgtol = 1e-5, int maxls = 20) {
    size_t n = x.size();
    auto clip = [&](double v) { return std::min(hi, std::max(lo, v)); };
    for (auto& v : x) v = clip(v);
    std::vector<double> g, gn, xn(n), d(n), hd(n);
    double fx = obj(x, g);
    std::vector<std::vector<double>> S, Y;
    std::vector<double> rho;
    int restarts = 0;
    auto set_diag = [&]() { for (size_t j = 0; j < n; j++) hd[j] = 1.0 / std::max(obj.f[j], 1e-300); };
    auto first_direction = [&]() {
        double gmax = 0;
        for (size_t j = 0; j < n; j++) gmax = std::max(gmax, std::fabs(hd[j] * g[j]));
        double sc = gmax > 0 ? 1.0 / gmax : 1.0;
        for (size_t j = 0; j < n; j++) d[j] = clip(x[j] - sc * hd[j] * g[j]) - x[j];
    };
    for (int it = 0; it < maxiter; it++) {
        double pgmax = 0;
        for (size_t i = 0; i < n; i++) {
            double gi = g[i];
            if ((x[i] <= lo && gi > 0) || (x[i] >= hi && gi < 0)) gi = 0;
            pgmax = std::max(pgmax, std::fabs(gi));
        }
        if (pgmax <= pgtol) { if (getenv("MS_DEBUG")) fprintf(stderr, "  lbfgs pgtol stop it %d f %.6g\n", it, fx); break; }
        if (it == 0) set_diag();   // one metric per call (a fixed variable scaling keeps the pairs consistent)
        int h = (int)S.size();
        if (h == 0) {
            first_direction();
        } else {
            std::vector<double> q = g;
            std::vector<double> a(h);
            for (int i = h - 1; i >= 0; i--) {
                a[i] = rho[i] * dot(S[i].data(), q.data(), n);
                for (size_t j = 0; j < n; j++) q[j] -= a[i] * Y[i][j];
            }
            // gamma for the diagonal initial matrix: s'y / y'Hy
            double yhy = 0;
            for (size_t j = 0; j < n; j++) yhy += Y[h - 1][j] * hd[j] * Y[h - 1][j];
            double gamma = dot(S[h - 1].data(), Y[h - 1].data(), n) / std::max(yhy, 1e-300);
            for (size_t j = 0; j < n; j++) q[j] *= gamma * hd[j];
            for (int i = 0; i < h; i++) {
                double b = rho[i] * dot(Y[i].data(), q.data(), n);
                for (size_t j = 0; j < n; j++) q[j] += S[i][j] * (a[i] - b);
            }
            for (size_t j = 0; j < n; j++) d[j] = clip(x[j] - q[j]) - x[j];
        }
        double dg = dot(d.data(), g.data(), n);
        if (!(dg < 0)) {
            S.clear(); Y.clear(); rho.clear();
            first_direction();
            dg = dot(d.data(), g.data(), n);
            if (!(dg < 0)) break;
        }
        // largest step inside the box
        double stpmx = 1e10;
        for (size_t j = 0; j < n; j++) {
            if (d[j] > 0) stpmx = std::min(stpmx, (hi - x[j]) / d[j]);
            else if (d[j] < 0) stpmx = std::min(stpmx, (lo - x[j]) / d[j]);
        }
        // L-BFGS-B (lnsrlb): with constraints the first iteration's step is at most 1
        if (it == 0) stpmx = 1.0;
        Dcsrch ls;
        ls.stpmax = std::max(stpmx, 0.0);
        double stp = std::min(1.0, ls.stpmax);
        Dcsrch::Task task = ls.start(stp, fx, dg);
        double fn = fx;
        bool ok = false;
        int nev = 0;
        if (task == Dcsrch::FG) {
            while (nev < maxls) {
                for (size_t j = 0; j < n; j++) xn[j] = clip(x[j] + stp * d[j]);
                fn = obj(xn, gn);
                nev++;
                double stp_eval = stp;
                double gd = dot(gn.data(), d.data(), n);
                if (!std::isfinite(fn)) { fn = 1e300; gd = 1e300; }
                task = ls.iterate(stp, fn, gd);
                if (task == Dcsrch::CONV || task == Dcsrch::WARN) {
                    if (stp != stp_eval) {   // the best point of the search: evaluate it
                        for (size_t j = 0; j < n; j++) xn[j] = clip(x[j] + stp * d[j]);
                        fn = obj(xn, gn);
                        nev++;
                    }
                    ok = std::isfinite(fn) && fn < fx;
                    break;
                }
                if (task == Dcsrch::ERR) break;
            }
        }
        if (!ok) {
            // no acceptable step: restart from the gradient once, else stop
            if (getenv("MS_DEBUG")) fprintf(stderr, "  lbfgs line search failed it %d f %.6g (restarts %d)\n", it, fx, restarts);
            if (S.empty() || restarts >= 2) break;
            S.clear(); Y.clear(); rho.clear();
            restarts++;
            continue;
        }
        std::vector<double> s(n), yk(n);
        for (size_t j = 0; j < n; j++) { s[j] = xn[j] - x[j]; yk[j] = gn[j] - g[j]; }
        double sy = dot(s.data(), yk.data(), n);
        if (sy > 2.220446049250313e-16 * dot(yk.data(), yk.data(), n)) {
            if ((int)S.size() >= m) { S.erase(S.begin()); Y.erase(Y.begin()); rho.erase(rho.begin()); }
            S.push_back(std::move(s)); Y.push_back(std::move(yk)); rho.push_back(1.0 / sy);
        }
        double fold = fx;
        x.swap(xn); g.swap(gn); fx = fn;
        if ((fold - fx) / std::max({std::fabs(fold), std::fabs(fx), 1.0}) <= ftol) { if (getenv("MS_DEBUG")) fprintf(stderr, "  lbfgs ftol stop it %d f %.6g\n", it, fx); break; }
        if (getenv("MS_DEBUG") && it == maxiter - 1) fprintf(stderr, "  lbfgs maxiter f %.6g\n", fx);
    }
}

}  // namespace
// ============================================================ spectrum preparation
// (shared with unidec.cpp / isodec.cpp through e2_internal.h)
namespace ms {
namespace e2 {

Spec restrict_spec(const double* mz, const double* it, long n, double lo, double hi) {
    std::vector<size_t> idx;
    idx.reserve(n);
    for (long i = 0; i < n; i++) {
        if (!std::isfinite(mz[i]) || !std::isfinite(it[i])) continue;
        if (!(mz[i] > 0)) continue;
        if (it[i] < 0) throw std::invalid_argument("Spectrum intensities must not be negative");
        idx.push_back((size_t)i);
    }
    std::stable_sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return mz[a] < mz[b]; });
    Spec s;
    bool range = hi > lo;
    for (size_t i : idx) {
        if (range && (mz[i] < lo || mz[i] > hi)) continue;
        s.x.push_back(mz[i]); s.y.push_back(it[i]);
    }
    return s;
}

// scipy 1-d filters with mode "nearest": window [i - k/2, i - k/2 + k - 1]
template <class Op> void filter_nearest(const std::vector<double>& in, std::vector<double>& out, long k, Op op, double init) {
    long n = (long)in.size();
    out.assign(n, 0.0);
    long s1 = k / 2;
    for (long i = 0; i < n; i++) {
        double v = init;
        for (long j = i - s1; j < i - s1 + k; j++) {
            long jj = std::min(n - 1, std::max(0L, j));
            v = op(v, in[jj]);
        }
        out[i] = v;
    }
}

void subtract_baseline(Spec& s, double factor) {
    size_t n = s.x.size();
    if (n < 20) return;
    double width = factor * ms::e2::estimate_peak_width(s.x.data(), s.y.data(), n);
    double dx = ms::e2::median_diff(s.x.data(), n);
    long k = (long)std::max(5.0, std::min((double)(n / 4), ms::e2::py_round(width / std::max(dx, 1e-12))));
    std::vector<double> b1, b2, b3;
    filter_nearest(s.y, b1, k, [](double a, double b) { return std::min(a, b); }, std::numeric_limits<double>::infinity());
    filter_nearest(b1, b2, k, [](double a, double b) { return std::max(a, b); }, -std::numeric_limits<double>::infinity());
    long ku = std::max(3L, k / 2);
    // uniform filter: running sum
    b3.assign(n, 0.0);
    {
        long s1 = ku / 2;
        for (long i = 0; i < (long)n; i++) {
            double v = 0;
            for (long j = i - s1; j < i - s1 + ku; j++) v += b2[std::min((long)n - 1, std::max(0L, j))];
            b3[i] = v / ku;
        }
    }
    for (size_t i = 0; i < n; i++) s.y[i] = std::max(0.0, s.y[i] - std::min(b3[i], s.y[i]));
}

// scipy.signal.find_peaks(y, height=thr, plateau_size=(1, None)) + peak_prominences bases
// -> mask of points of the peaks reaching thr, number of peaks (ms_deconv.peak_mask)
std::pair<std::vector<char>, long> peak_mask(const Spec& s, double thr) {
    const std::vector<double>& x = s.x; const std::vector<double>& y = s.y;
    long n = (long)y.size();
    std::vector<char> keep(n, 0);
    if (n < 3) {
        long c = 0;
        for (long i = 0; i < n; i++) if (y[i] >= thr) { keep[i] = 1; c++; }
        return {keep, c};
    }
    std::vector<long> idx;
    long i = 1, imax = n - 1;
    while (i < imax) {
        if (y[i - 1] < y[i]) {
            long ahead = i + 1;
            while (ahead < imax && y[ahead] == y[i]) ahead++;
            if (y[ahead] < y[i]) {
                long mid = (i + ahead - 1) / 2;
                if (y[mid] >= thr) idx.push_back(mid);
                i = ahead;
            }
        }
        i++;
    }
    if (idx.empty()) return {keep, 0};
    double fw = ms::e2::estimate_peak_width(x.data(), y.data(), n);
    long i0 = std::max_element(y.begin(), y.end()) - y.begin();
    long lo = std::max(0L, i0 - 50), hi = std::min(n - 1, i0 + 50);
    double dx = hi > lo ? ms::e2::median_diff(x.data() + lo, hi - lo + 1) : ms::e2::median_diff(x.data(), n);
    double n_fw = dx > 0 ? fw / dx : 5.0;
    long wlen = (long)std::max(7.0, std::min(10 * n_fw, (double)n)) | 1;
    for (long pk : idx) {
        long imin = 0, imx = n - 1;
        if (wlen >= 2) { imin = std::max(pk - wlen / 2, imin); imx = std::min(pk + wlen / 2, imx); }
        long lb = pk, rb = pk;
        double lmin = y[pk];
        for (long j = pk; j >= imin && y[j] <= y[pk]; j--) if (y[j] < lmin) { lmin = y[j]; lb = j; }
        double rmin = y[pk];
        for (long j = pk; j <= imx && y[j] <= y[pk]; j++) if (y[j] < rmin) { rmin = y[j]; rb = j; }
        for (long j = lb; j <= rb; j++) keep[j] = 1;
    }
    return {keep, (long)idx.size()};
}

std::string fmt_int(double v) {
    char buf[64];
    double a = std::fabs(v);
    if (a >= 1e5) {
        snprintf(buf, sizeof buf, "%.2e", v);
        std::string s = buf;
        size_t p = s.find("e+0");
        if (p != std::string::npos) s.replace(p, 3, "e");
        p = s.find("e+");
        if (p != std::string::npos) s.replace(p, 2, "e");
        return s;
    }
    if (a >= 100) snprintf(buf, sizeof buf, "%.0f", v);
    else snprintf(buf, sizeof buf, "%.3g", v);
    return buf;
}

}  // namespace e2
}  // namespace ms

namespace {

double averagine_sigma(double mass) {
    double n = mass / 111.1254;
    struct El { double cnt; std::vector<std::pair<double, double>> iso; };
    std::vector<El> els = {
        {4.9384 * n, {{0, 0.9893}, {1.00336, 0.0107}}},
        {7.7583 * n, {{0, 0.999885}, {1.00628, 0.000115}}},
        {1.3577 * n, {{0, 0.99636}, {0.99704, 0.00364}}},
        {1.4773 * n, {{0, 0.99757}, {1.00422, 0.00038}, {2.00425, 0.00205}}},
        {0.0417 * n, {{0, 0.9499}, {0.99939, 0.0075}, {1.99580, 0.0425}}}};
    double var = 0;
    for (auto& e : els) {
        double m1 = 0, m2 = 0;
        for (auto& p : e.iso) { m1 += p.first * p.second; m2 += p.first * p.first * p.second; }
        var += e.cnt * (m2 - m1 * m1);
    }
    return std::sqrt(var);
}

// FWHM (m/z) for the linear grid: the median over the tallest well separated
// peaks (up to five, at least two widths apart), as the resolving power of the
// log grid is measured; the tallest peak alone is a poor sample on binned
// quadrupole data, whose apparent widths vary with the bin the centroid falls in.
double peak_width_median(const double* x, const double* y, size_t n) {
    double w0 = ms::e2::estimate_peak_width(x, y, n);
    if (n < 5 || !(w0 > 0)) return w0;
    std::vector<size_t> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::stable_sort(order.begin(), order.end(), [y](size_t a, size_t b) { return y[a] > y[b]; });
    double top = y[order[0]];
    std::vector<double> used, vals;
    for (size_t t = 0; t < std::min<size_t>(5000, n); t++) {
        size_t i = order[t];
        if (y[i] < 0.05 * top || vals.size() >= 5) break;
        bool near = false;
        for (double u : used) if (std::fabs(x[i] - u) < 2.0 * w0) { near = true; break; }
        if (near) continue;
        used.push_back(x[i]);
        double w = ms::e2::fwhm_at(x, y, n, i);
        if (w > 0 && std::isfinite(w)) vals.push_back(w);
    }
    return vals.empty() ? w0 : ms::e2::median(vals);
}

}  // namespace
// ============================================================ result object
struct ms_maxent_result {
    std::vector<double> mass_x, mass_y, fit_x, fit_y, z, zdist;
    std::vector<double> pk_mass, pk_height, pk_area, pk_apex;
    std::vector<int> pk_niso;
    double r2 = 0, chi2 = 0;
    int rounds = 0;
    std::string notes;
};

namespace {

// ============================================================ noise model
// sigma_i^2 = s0^2 + ybar_i / g + (eps ybar_i)^2 in the normalised units of
// the gridded data (base peak = 1); ybar is the three point running mean of
// the data (an estimate of the expected intensity that does not favour the
// points where the noise happened to dip).
struct NoiseModel {
    double s0 = 1e-3;     // floor
    double ginv = 0;      // 1 / g, variance per unit of intensity (counting)
    double eps = 0;       // relative mismatch of the peak shape
    double F = 1.0;       // variance of a grid cell relative to a raw point (for the notes)
    bool from_data = false;
    double var(double ybar) const { return s0 * s0 + ginv * ybar + eps * eps * ybar * ybar; }
};

// Weighted least squares of v ~ sum_p x_p basis_p over a few points with the
// non negative parameters chosen by trying every active set (3 parameters at
// most). basis is nb x np row major; returns x (np).
std::vector<double> nnls_tiny(const std::vector<double>& basis, const std::vector<double>& v, const std::vector<double>& wt, int nb, int np) {
    std::vector<double> best(np, 0.0);
    double best_sse = std::numeric_limits<double>::infinity();
    for (int mask = 0; mask < (1 << np); mask++) {
        std::vector<int> idx;
        for (int p = 0; p < np; p++) if (mask & (1 << p)) idx.push_back(p);
        int m = (int)idx.size();
        std::vector<double> x(np, 0.0);
        if (m > 0) {
            std::vector<double> A((size_t)m * m, 0.0), rhs(m, 0.0);
            for (int i = 0; i < nb; i++)
                for (int a = 0; a < m; a++) {
                    rhs[a] += wt[i] * basis[(size_t)i * np + idx[a]] * v[i];
                    for (int c = 0; c < m; c++) A[(size_t)a * m + c] += wt[i] * basis[(size_t)i * np + idx[a]] * basis[(size_t)i * np + idx[c]];
                }
            // Gaussian elimination with partial pivoting
            bool ok = true;
            for (int c = 0; c < m && ok; c++) {
                int piv = c;
                for (int r = c + 1; r < m; r++) if (std::fabs(A[(size_t)r * m + c]) > std::fabs(A[(size_t)piv * m + c])) piv = r;
                if (!(std::fabs(A[(size_t)piv * m + c]) > 1e-300)) { ok = false; break; }
                if (piv != c) { for (int t = 0; t < m; t++) std::swap(A[(size_t)c * m + t], A[(size_t)piv * m + t]); std::swap(rhs[c], rhs[piv]); }
                for (int r = c + 1; r < m; r++) {
                    double fct = A[(size_t)r * m + c] / A[(size_t)c * m + c];
                    for (int t = c; t < m; t++) A[(size_t)r * m + t] -= fct * A[(size_t)c * m + t];
                    rhs[r] -= fct * rhs[c];
                }
            }
            if (!ok) continue;
            std::vector<double> sol(m);
            for (int c = m - 1; c >= 0; c--) {
                double s = rhs[c];
                for (int t = c + 1; t < m; t++) s -= A[(size_t)c * m + t] * sol[t];
                sol[c] = s / A[(size_t)c * m + c];
            }
            bool feasible = true;
            for (int a = 0; a < m; a++) { if (!(sol[a] >= 0) || !std::isfinite(sol[a])) feasible = false; x[idx[a]] = sol[a]; }
            if (!feasible) continue;
        }
        double sse = 0;
        for (int i = 0; i < nb; i++) {
            double fit = 0;
            for (int p = 0; p < np; p++) fit += basis[(size_t)i * np + p] * x[p];
            sse += wt[i] * (v[i] - fit) * (v[i] - fit);
        }
        if (sse < best_sse) { best_sse = sse; best = x; }
    }
    return best;
}

// Floor and counting term from the raw spectrum: in windows of W points
// without peaks (judged by the second differences: their rms against their
// robust scale), the variance of the second differences / 6 against the mean
// intensity is regressed as a + b x mean (a, b >= 0), after taking medians
// in eight bins of the mean. The result is converted to the gridded data
// (the grid cells average F raw spacings; the variance of an average of
// a linear interpolant over a cell is about 2/3 of the raw one).
bool noise_from_data(const ms::e2::Spec& s, const std::vector<char>& keep, double scale, double dx_raw, double d_mz, NoiseModel& nm) {
    const size_t n = s.y.size();
    const size_t W = n >= 4000 ? 16 : 8;
    if (n < 20 * W || !(scale > 0)) return false;
    std::vector<double> mw, vw, d2;
    for (size_t a = 0; a + W <= n; a += W) {
        bool ok = true; double mean = 0;
        for (size_t i = a; i < a + W; i++) {
            if ((!keep.empty() && !keep[i]) || !(s.y[i] > 0)) { ok = false; break; }
            mean += s.y[i];
        }
        if (!ok) continue;
        mean /= W;
        d2.clear();
        for (size_t i = a + 1; i + 1 < a + W; i++) d2.push_back(s.y[i - 1] - 2.0 * s.y[i] + s.y[i + 1]);
        double med = ms::e2::median(d2);
        std::vector<double> ad(d2.size());
        for (size_t i = 0; i < d2.size(); i++) ad[i] = std::fabs(d2[i] - med);
        double sd = 1.4826 * ms::e2::median(ad) / std::sqrt(6.0);
        double rms = 0;
        for (double v : d2) rms += v * v;
        rms = std::sqrt(rms / d2.size() / 6.0);
        if (!(sd > 0) || rms > 2.5 * sd) continue;   // a peak inside the window
        mw.push_back(mean); vw.push_back(sd * sd);
    }
    if (mw.size() < 40) return false;
    std::vector<size_t> order(mw.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return mw[a] < mw[b]; });
    const int nb = 8;
    std::vector<double> bm, bv, bw;
    for (int b = 0; b < nb; b++) {
        size_t i0 = order.size() * b / nb, i1 = order.size() * (b + 1) / nb;
        if (i1 <= i0) continue;
        std::vector<double> m1, v1;
        for (size_t t = i0; t < i1; t++) { m1.push_back(mw[order[t]]); v1.push_back(vw[order[t]]); }
        bm.push_back(ms::e2::median(m1)); bv.push_back(ms::e2::median(v1)); bw.push_back((double)(i1 - i0));
    }
    if (bm.size() < 3) return false;
    std::vector<double> basis(bm.size() * 2);
    for (size_t i = 0; i < bm.size(); i++) { basis[i * 2] = 1.0; basis[i * 2 + 1] = bm[i]; }
    // relative least squares: every bin counts by its number of windows
    std::vector<double> wt(bm.size());
    for (size_t i = 0; i < bm.size(); i++) wt[i] = bw[i] / std::max(bv[i] * bv[i], 1e-300);
    std::vector<double> ab = nnls_tiny(basis, bv, wt, (int)bm.size(), 2);
    if (!(ab[0] >= 0) || !(ab[1] >= 0)) return false;
    const double F = (d_mz > 0 && dx_raw > 0) ? std::min(2.0 / 3.0, dx_raw / d_mz) : 2.0 / 3.0;
    nm.s0 = std::sqrt(F * ab[0]) / scale;
    nm.ginv = F * ab[1] / scale;
    nm.eps = 0;
    nm.F = F;
    nm.from_data = true;
    return true;
}

// Re-estimation from the residuals of a fit: in bins of ybar (log spaced
// above 3 s0), the robust variance of the residuals (1.4826 x median |r|)^2
// against ybar; the excess over `base` (the terms measured on the data are
// lower bounds) is fitted as da + db ybar + c ybar^2, all non negative,
// in relative least squares weighted by the points per bin. c (the shape
// term) is only fitted for kind 3: with it the few points of the strongest
// peaks count for nothing in chi squared and the entropy shrinks those
// peaks by a third (measured on prot_spec), so it is not the default.
NoiseModel noise_from_residuals(const std::vector<double>& ybar, const std::vector<double>& r, const std::vector<char>& signal,
                                const NoiseModel& base, int kind, std::string* dbg) {
    const size_t n = ybar.size();
    const double lo = std::max(3.0 * base.s0, 1e-6);
    const int nb = 10;
    std::vector<double> edges;
    edges.push_back(0.0);
    for (int b = 0; b <= nb; b++) edges.push_back(lo * std::pow(1.0 / lo, (double)b / nb));
    std::vector<double> ym, vm, cnt;
    std::vector<std::vector<double>> absr(edges.size() - 1);
    std::vector<double> ysum(edges.size() - 1, 0.0);
    for (size_t i = 0; i < n; i++) {
        if (!signal[i]) continue;
        double yb = ybar[i];
        size_t k = std::upper_bound(edges.begin(), edges.end(), yb) - edges.begin();
        if (k == 0) continue;
        k--;
        if (k >= absr.size()) k = absr.size() - 1;
        absr[k].push_back(std::fabs(r[i]));
        ysum[k] += yb;
    }
    for (size_t k = 0; k < absr.size(); k++) {
        if (absr[k].size() < 30) continue;
        double med = ms::e2::median(absr[k]);
        double sd = 1.4826 * med;
        ym.push_back(ysum[k] / absr[k].size()); vm.push_back(sd * sd); cnt.push_back((double)absr[k].size());
    }
    NoiseModel nm = base;
    if (ym.size() < 2) {
        // too few populated bins (a few kept peaks after the minimum intensity): one robust scale of all signal residuals
        std::vector<double> all;
        for (size_t i = 0; i < n; i++) if (signal[i]) all.push_back(std::fabs(r[i]));
        if (all.size() >= 10) {
            double sd = 1.4826 * ms::e2::median(all);
            double mean_var = 0; size_t cnt = 0;
            for (size_t i = 0; i < n; i++) if (signal[i]) { mean_var += base.var(ybar[i]); cnt++; }
            mean_var /= std::max<size_t>(cnt, 1);
            if (sd * sd > mean_var) nm.s0 = std::sqrt(base.s0 * base.s0 + sd * sd - mean_var);
        }
        return nm;
    }
    const int np = kind == 3 ? 3 : 2;
    std::vector<double> basis(ym.size() * np), ex(ym.size()), wt(ym.size());
    for (size_t i = 0; i < ym.size(); i++) {
        basis[i * np] = 1.0; basis[i * np + 1] = ym[i];
        if (np == 3) basis[i * np + 2] = ym[i] * ym[i];
        ex[i] = vm[i] - base.var(ym[i]);
        wt[i] = cnt[i] / std::max(vm[i] * vm[i], 1e-300);
    }
    std::vector<double> x = nnls_tiny(basis, ex, wt, (int)ym.size(), np);
    nm.s0 = std::sqrt(base.s0 * base.s0 + std::max(0.0, x[0]));
    nm.ginv = base.ginv + std::max(0.0, x[1]);
    nm.eps = np == 3 ? std::sqrt(std::max(0.0, x[2])) : 0.0;
    if (dbg) {
        char b[256];
        for (size_t i = 0; i < ym.size(); i++) {
            snprintf(b, sizeof b, "    bin ybar %.3g n %.0f sd %.3g (model %.3g -> %.3g)\n", ym[i], cnt[i], std::sqrt(vm[i]), std::sqrt(base.var(ym[i])), std::sqrt(nm.var(ym[i])));
            *dbg += b;
        }
    }
    return nm;
}

// ============================================================ alpha search
// Points (ln alpha, ln(chi2n / target)) of one pass. The next ln alpha: a
// regula falsi step (Illinois rule: the retained end's value is halved when
// the same end is kept twice) once chi squared has been seen on both sides
// of the target, else a slope limited step towards it.
struct AlphaSearch {
    std::vector<std::pair<double, double>> pts;
    bool have_hi = false, have_lo = false;      // hi: chi2 above the target, lo: below
    double la_hi = 0, lc_hi = 0, la_lo = 0, lc_lo = 0;
    int side = 0;                               // the end kept at the last update (+1 hi, -1 lo)
    void reset() { pts.clear(); have_hi = have_lo = false; side = 0; }
    void add(double la, double lc) {
        pts.push_back({la, lc});
        if (lc > 0) {
            have_hi = true; la_hi = la; lc_hi = lc;
            if (have_lo) { if (side == -1) lc_lo *= 0.5; side = -1; }   // the lo end is kept
        } else {
            have_lo = true; la_lo = la; lc_lo = lc;
            if (have_hi) { if (side == 1) lc_hi *= 0.5; side = 1; }     // the hi end is kept
        }
    }
    double bracket_width() const {
        if (!have_hi || !have_lo) return std::numeric_limits<double>::infinity();
        return std::fabs(la_hi - la_lo);
    }
    double next() const {
        const auto& last = pts.back();
        double la = last.first, lc = last.second;
        if (have_hi && have_lo) {
            double t = lc_hi / (lc_hi - lc_lo);           // regula falsi in ln alpha
            t = std::min(0.9, std::max(0.1, t));          // stay inside the bracket
            return la_hi + t * (la_lo - la_hi);
        }
        double slope = 0.5;                                // d ln chi2 / d ln alpha, typical
        if (pts.size() >= 2) {
            const auto& prev = pts[pts.size() - 2];
            double dla = la - prev.first, dlc = lc - prev.second;
            if (std::fabs(dla) > 1e-6) {
                double sl = dlc / dla;
                // on a plateau (chi squared does not react to alpha) take the largest step
                slope = sl < 0.05 ? 0.0 : std::min(2.0, std::max(0.1, sl));
            }
        }
        if (!(slope > 0)) return la + (lc > 0 ? -4.0 : 4.0);   // plateau: a factor 55 per round
        double step = -lc / slope;
        step = std::min(3.0, std::max(-3.0, step));
        return la + step;
    }
};

// ============================================================ driver
void run_maxent(const double* mz, const double* it, long n0, const ms_maxent_params& p, ms_maxent_result& out,
                ms_progress_fn cb, void* user) {
    using namespace ms::e2;
    const bool debug = getenv("MS_DEBUG") != nullptr;
    // ---- check_grid
    double lo = p.mass_lo, hi = p.mass_hi, step = p.mass_step;
    double z0 = p.z_lo, z1 = p.z_hi;
    if (!std::isfinite(lo) || !std::isfinite(hi) || !std::isfinite(step))
        throw std::invalid_argument("Mass, charge and step values must be finite numbers");
    if (lo <= 0 || z0 < 1 || z1 < z0) throw std::invalid_argument("Use positive masses and whole-number charges in increasing order");
    if (step <= 0) throw std::invalid_argument("The mass step must be larger than 0");
    if (step < 0.001) {
        char b[256];
        snprintf(b, sizeof b, "A mass step of %g Da is finer than any mass spectrometer resolves (it would only make the calculation extremely slow); use 0.005 Da or more", step);
        throw std::invalid_argument(b);
    }
    if (hi <= lo) throw std::invalid_argument("The upper mass must be larger than the lower mass");
    double nm_est = (hi - lo) / step;
    if (nm_est > 2e9) throw std::invalid_argument("Mass range / mass step gives more than two billion mass points");
    const double budget = ms::memory_budget();
    {
        // the model alone (two entries of 10 bytes per mass and charge, twice): a first check before any allocation
        const double need = 40.0 * nm_est * (z1 - z0 + 1);
        if (need > budget) {
            char b[600];
            snprintf(b, sizeof b, "Mass range %g to %g Da every %g Da with %d charges gives %.0f x %d grid points and would need about %s of memory; this computer can give it %s. Use a larger mass step (isotope peaks need 0.01 to 0.02 Da, envelopes 0.5 to 1 Da) or narrower mass and charge ranges.",
                     lo, hi, step, (int)(z1 - z0 + 1), nm_est, (int)(z1 - z0 + 1), ms::fmt_bytes(need).c_str(), ms::fmt_bytes(budget).c_str());
            throw std::invalid_argument(b);
        }
    }
    int rounds = p.rounds > 0 ? p.rounds : 10;
    int iterations = p.iterations > 0 ? p.iterations : 50;
    int min_rounds = p.min_rounds > 0 ? p.min_rounds : 4;
    if (p.rounds < 0 || p.iterations < 0 || p.min_rounds < 0) throw std::invalid_argument("rounds must be a positive whole number");
    if (!std::isfinite(p.peak_width) || p.peak_width < 0) throw std::invalid_argument("peak_width must be finite and non-negative");
    const int noise_kind = p.noise_model >= 0 && p.noise_model <= 3 ? p.noise_model : 0;
    const double chi2_target = p.chi2_target > 0 && std::isfinite(p.chi2_target) ? p.chi2_target : 1.0;
    if (p.envelope_nodes < 0 || p.envelope_nodes > 64) throw std::invalid_argument("envelope_nodes must be between 0 (automatic) and 64");

    // ---- run_method: restriction, minimum intensity, baseline
    Spec spec = restrict_spec(mz, it, n0, p.mz_lo, p.mz_hi);
    size_t n = spec.x.size();
    double base = n ? *std::max_element(spec.y.begin(), spec.y.end()) : 0.0;
    double thr = 0;
    if (p.min_intensity > 0) thr = p.min_intensity_pct ? p.min_intensity / 100.0 * base : p.min_intensity;
    std::vector<char> keep;
    long npk = 0;
    if (thr > 0 && n) {
        auto pm = peak_mask(spec, thr);
        keep = pm.first; npk = pm.second;
    } else {
        keep.assign(n, 0);
        for (size_t i = 0; i < n; i++) keep[i] = spec.y[i] > 0;
        if (n > 2) for (size_t i = 1; i + 1 < n; i++) if (spec.y[i] > spec.y[i - 1] && spec.y[i] >= spec.y[i + 1] && spec.y[i] > 0) npk++;
    }
    long points_used = 0;
    for (char k : keep) points_used += k;
    std::string summary;
    if (n) {
        char b[512];
        snprintf(b, sizeof b, "m/z %.2f to %.2f", spec.x.front(), spec.x.back());
        summary = b;
        if (thr > 0) {
            double pct = base > 0 ? 100.0 * thr / base : 0.0;
            snprintf(b, sizeof b, ", peaks from %s up (%.3g %% of the base peak %s): %ld peaks kept whole, %ld of %ld points",
                     fmt_int(thr).c_str(), pct, fmt_int(base).c_str(), npk, points_used, (long)n);
            summary += b;
        } else {
            snprintf(b, sizeof b, ", no minimum intensity (base peak %s)", fmt_int(base).c_str());
            summary += b;
        }
    }
    if (p.baseline) subtract_baseline(spec, p.baseline_width > 0 ? p.baseline_width : 15.0);
    if (thr > 0) {
        bool anypos = false;
        for (size_t i = 0; i < n; i++) { if (!keep[i]) spec.y[i] = 0.0; else if (spec.y[i] > 0) anypos = true; }
        if (npk < 1 || points_used < 5 || !anypos) {
            throw std::invalid_argument("No data above the minimum intensity (" + fmt_int(thr) + "; the base peak in the m/z range is " +
                                        fmt_int(base) + "). Lower the minimum intensity.");
        }
    }
    // ---- _maxent
    if (n < 10) throw std::invalid_argument("Too few data points in the m/z range");
    {
        bool anypos = false;
        for (double v : spec.y) if (v > 0) { anypos = true; break; }
        if (!anypos) throw std::invalid_argument("No positive signal remains for maximum entropy deconvolution");
    }
    std::vector<double> x(n);
    for (size_t i = 0; i < n; i++) x[i] = std::max(spec.x[i], 1e-6);
    int zlo = p.z_lo, zhi = p.z_hi;
    std::vector<int> charges;
    for (int z = std::max(1, zlo); z <= std::max(zlo, zhi); z++) charges.push_back(z);
    int K = (int)charges.size();
    std::vector<double> masses = arange(lo, hi + step / 2, step);
    if (masses.size() < 3) throw std::invalid_argument("The mass range is shorter than three mass steps");
    bool resolved = p.resolved_isotopes != 0;
    double res_power = p.resolution > 0 ? p.resolution : 0.0;
    bool use_log = p.resolution >= 0;
    double a = p.adduct_mass != 0 && std::isfinite(p.adduct_mass) ? p.adduct_mass : PROTON;
    double limit = 1.5e6;
    std::string note_shape;
    std::vector<double> coord, grid, y, ymask;
    Model model;
    double d, fwhm, rel = 0;
    std::string iso_mode = resolved ? "resolved" : "envelope";
    // raw mask of the points the minimum intensity removed (1 = removed), averaged into the grid cells like the data
    std::vector<double> removed;
    if (thr > 0) { removed.assign(n, 0.0); for (size_t i = 0; i < n; i++) removed[i] = keep[i] ? 0.0 : 1.0; }
    auto cell_average = [&](const std::vector<double>& xs, const std::vector<double>& vals, const std::vector<double>& centres, double width, std::vector<double>& outv) {
        std::vector<double> cum(n, 0.0);
        for (size_t i = 1; i < n; i++) cum[i] = cum[i - 1] + 0.5 * (vals[i] + vals[i - 1]) * (xs[i] - xs[i - 1]);
        outv.resize(centres.size());
        for (size_t i = 0; i < centres.size(); i++)
            outv[i] = (interp_clamp(centres[i] + width / 2, xs.data(), cum.data(), n) - interp_clamp(centres[i] - width / 2, xs.data(), cum.data(), n)) / width;
    };
    // envelope nodes: automatic from the mass range (one per factor 1.5 in mass, at most 12), limited by the cost of the envelope fit
    auto node_count = [&](size_t ng) {
        int Qn = p.envelope_nodes;
        if (Qn <= 0) {
            double ratio = masses.back() / std::max(masses.front(), 1e-9);
            Qn = 1 + (int)std::floor(std::log(std::max(ratio, 1.0)) / std::log(1.5) + 1e-9);
            Qn = std::max(1, std::min(12, Qn));
        }
        // the normal equations cost (Q K)^2 x grid points per envelope fit
        int qcost = (int)std::floor(std::sqrt(3e9 / std::max<double>(ng, 1.0)) / std::max(K, 1));
        if (p.envelope_nodes <= 0) Qn = std::max(1, std::min(Qn, qcost));
        return Qn;
    };
    if (use_log) {
        if (res_power <= 0) res_power = estimate_resolution(spec.x.data(), spec.y.data(), n);
        rel = 1.0 / res_power;
        std::string switched;
        if (iso_mode == "resolved") {
            double fwhm_lo = masses[0] / res_power;
            if (step > 0.8 * fwhm_lo) {
                iso_mode = "envelope";
                char b[512];
                snprintf(b, sizeof b, " (the mass step of %g Da is too coarse for resolved isotope peaks, %.3f Da wide at %g Da: isotope envelopes were fitted instead; use a step of %s Da and a narrower mass range for resolved isotopes)",
                         step, fwhm_lo, masses[0], fwhm_lo > 0.1 ? "0.01 to 0.05" : "0.01");
                switched = b;
            }
        }
        char b[512];
        if (iso_mode == "envelope") {
            double mc = 0.5 * (masses.front() + masses.back());
            rel = std::sqrt(rel * rel + std::pow(2.3548 * averagine_sigma(mc) / mc, 2));
            snprintf(b, sizeof b, "isotope envelopes (average masses), resolving power %.0f", res_power);
            note_shape = std::string(b) + switched;
        } else {
            snprintf(b, sizeof b, "resolved isotopes, resolving power %.0f", res_power);
            note_shape = b;
        }
        std::vector<double> lx(n);
        for (size_t i = 0; i < n; i++) lx[i] = std::log(x[i]);
        double du = rel / 6.0;
        // windows
        std::vector<std::pair<double, double>> iv;
        for (int z : charges) {
            double p0 = (masses.front() + z * a) / z, p1 = (masses.back() + z * a) / z;
            p0 = std::log(std::max(p0, 1e-9)); p1 = std::log(std::max(p1, 1e-9));
            p0 = std::max(p0 - 8 * rel, lx.front()); p1 = std::min(p1 + 8 * rel, lx.back());
            if (p1 > p0) iv.push_back({p0, p1});
        }
        std::sort(iv.begin(), iv.end());
        std::vector<std::pair<double, double>> win;
        for (auto& w : iv) {
            if (!win.empty() && w.first <= win.back().second) win.back().second = std::max(win.back().second, w.second);
            else win.push_back(w);
        }
        double est = 0;
        for (auto& w : win) est += std::max(0.0, std::ceil((w.second - w.first) / du + 0.5));
        if (est > limit) throw std::invalid_argument("The spectrum needs too many grid points at this resolution: narrow the ranges or use a lower resolving power");
        for (auto& w : win) { auto part = arange(w.first, w.second + du / 2, du); coord.insert(coord.end(), part.begin(), part.end()); }
        if (coord.size() > limit) throw std::invalid_argument("The spectrum needs too many grid points at this resolution");
        if (coord.size() < 10) throw std::invalid_argument("None of these charges puts the masses inside the m/z range of the spectrum");
        cell_average(lx, spec.y, coord, du, y);
        if (!removed.empty()) cell_average(lx, removed, coord, du, ymask);
        grid.resize(coord.size());
        for (size_t i = 0; i < coord.size(); i++) grid[i] = std::exp(coord[i]);
        d = du;
        model.build(coord, du, masses, charges, a, rel, true, budget, node_count(coord.size()));
        fwhm = rel * median(grid);
    } else {
        fwhm = p.peak_width > 0 ? p.peak_width : 0.0;
        if (fwhm <= 0) fwhm = peak_width_median(spec.x.data(), spec.y.data(), n);
        d = fwhm / 6.0;
        std::vector<std::pair<double, double>> iv;
        for (int z : charges) {
            double p0 = (masses.front() + z * a) / z, p1 = (masses.back() + z * a) / z;
            p0 = std::max(p0 - 8 * fwhm, x.front()); p1 = std::min(p1 + 8 * fwhm, x.back());
            if (p1 > p0) iv.push_back({p0, p1});
        }
        std::sort(iv.begin(), iv.end());
        std::vector<std::pair<double, double>> win;
        for (auto& w : iv) {
            if (!win.empty() && w.first <= win.back().second) win.back().second = std::max(win.back().second, w.second);
            else win.push_back(w);
        }
        double est = 0;
        for (auto& w : win) est += std::max(0.0, std::ceil((w.second - w.first) / d + 0.5));
        if (est > limit) throw std::invalid_argument("The spectrum needs too many grid points: narrow the ranges or use a wider peak width");
        for (auto& w : win) { auto part = arange(w.first, w.second + d / 2, d); coord.insert(coord.end(), part.begin(), part.end()); }
        if (coord.size() > limit) throw std::invalid_argument("The spectrum needs too many grid points");
        if (coord.size() < 10) throw std::invalid_argument("None of these charges puts the masses inside the m/z range of the spectrum");
        grid = coord;
        cell_average(x, spec.y, coord, d, y);
        if (!removed.empty()) cell_average(x, removed, coord, d, ymask);
        model.build(coord, d, masses, charges, a, fwhm, false, budget, node_count(coord.size()));
        char b[128];
        snprintf(b, sizeof b, "peak width %.3g m/z", fwhm);
        note_shape = b;
    }
    if (model.na < 1) throw std::invalid_argument("None of these charges puts the masses inside the m/z range of the spectrum");
    const size_t ng = y.size();
    double ymax = *std::max_element(y.begin(), y.end());
    double scale = ymax > 0 ? ymax : 1.0;
    if (!(ymax > 0)) throw std::invalid_argument("No positive signal falls inside the selected mass and charge ranges");
    for (auto& v : y) v /= scale;
    // three point running mean (expected intensity for the noise model)
    std::vector<double> ybar(ng);
    for (size_t i = 0; i < ng; i++) {
        double s = y[i]; int c = 1;
        if (i > 0) { s += y[i - 1]; c++; }
        if (i + 1 < ng) { s += y[i + 1]; c++; }
        ybar[i] = std::max(0.0, s / c);
    }
    // points the minimum intensity emptied: the data there are below the threshold, not zero
    std::vector<double> yeff = ybar;
    if (!ymask.empty()) {
        const double thr_n = thr / scale;
        for (size_t i = 0; i < ng; i++) if (ymask[i] > 0.5) yeff[i] = std::max(yeff[i], thr_n);
    }
    // ---- noise model: floor and counting term from the data
    double noise_old = 0.01;   // the 2.7 estimate: successive differences of the gridded data
    if (ng > 5) {
        std::vector<double> dy(ng - 1);
        for (size_t i = 0; i + 1 < ng; i++) dy[i] = y[i + 1] - y[i];
        double md = median(dy);
        std::vector<double> ad(ng - 1);
        for (size_t i = 0; i + 1 < ng; i++) ad[i] = std::fabs(dy[i] - md);
        noise_old = 1.4826 * median(ad) / std::sqrt(2.0);
    }
    NoiseModel noise;
    {
        double dx_raw = median_diff(spec.x.data(), n);
        double d_mz = use_log ? d * median(grid) : d;
        bool ok = noise_kind != 2 && noise_from_data(spec, thr > 0 ? keep : std::vector<char>(), scale, dx_raw, d_mz, noise);
        if (!ok) { noise.s0 = noise_old; noise.ginv = 0; noise.eps = 0; noise.from_data = false; }
        if (noise_kind == 2) noise.s0 = std::max(noise_old, 2e-3);
        else noise.s0 = std::max(noise.s0, 1e-4);
    }
    const NoiseModel noise_data = noise;   // lower bounds for the re-estimation
    // weights 1 / sigma_i^2 with the variance of the larger of the data and the
    // fit of the previous round (Pearson's form): where the model has intensity
    // but the data none, e.g. on the flanks of a peak whose shape is not the
    // model's, the counting term follows the model instead of leaving the floor
    // alone, which would otherwise make the flanks dearer than the whole peak
    std::vector<double> w(ng), yvar = yeff, fitv(ng, 0.0), yblur(ng);
    model.blur.apply(yeff.data(), yblur.data(), ng);   // the data seen through the peak shape: the flanks a peak reaches
    for (size_t i = 0; i < ng; i++) yeff[i] = std::max(yeff[i], yblur[i]);
    bool have_fit = false;
    auto set_weights = [&]() {
        for (size_t i = 0; i < ng; i++) {
            yvar[i] = have_fit ? std::max(yeff[i], fitv[i]) : yeff[i];
            w[i] = 1.0 / noise.var(yvar[i]);
        }
    };
    set_weights();
    // the stopping test uses the points with signal only: on a TOF spectrum
    // most of the grid is empty baseline, which would make the average chi
    // squared look fine while weak peaks are not fitted yet
    std::vector<char> signal(ng, 0);
    for (size_t i = 0; i < ng; i++) signal[i] = ybar[i] > 3.0 * noise.s0;
    dilate(signal, std::max(1L, (long)(3 * model.sig)));
    long n_sig = 0;
    for (char s : signal) n_sig += s;
    n_sig = std::max(n_sig, 1L);
    double total = std::max(np_sum(y.data(), ng), 1e-9);
    double dflt = total / masses.size() * 0.01;
    size_t na = model.na;
    double log_dflt = std::log(dflt);
    std::vector<double> u(na, log_dflt), f(na), G, bvec;
    // starting envelope: equal weights for the charges whose window holds
    // signal, per node (a least squares fit at this point would be sparse and
    // the masses would settle on a single charge each: ghosts at M x z'/z)
    {
        for (size_t j = 0; j < na; j++) f[j] = std::exp(u[j]);
        model.node_normal(f.data(), y.data(), w.data(), G, bvec);
        std::vector<double> c((size_t)model.Q * K, 1.0 / K);
        for (int q = 0; q < model.Q; q++) {
            double cmax = 0;
            for (int k = 0; k < K; k++) cmax = std::max(cmax, bvec[(size_t)q * K + k]);
            if (cmax > 0) {
                double s = 0;
                for (int k = 0; k < K; k++) { double v = bvec[(size_t)q * K + k] > 0.01 * cmax ? 1.0 : 1e-6; c[(size_t)q * K + k] = v; s += v; }
                for (int k = 0; k < K; k++) c[(size_t)q * K + k] /= s;
            }
        }
        model.set_envelope(c);
    }
    // starting alpha: a fraction of the strongest pull of the data on a mass
    // (the gradient of chi squared / 2 at the flat start)
    double alpha;
    {
        std::vector<double> wy(ng), pull(na);
        for (size_t i = 0; i < ng; i++) wy[i] = w[i] * y[i];
        model.adjoint_f(wy.data(), pull.data());
        double pmax = 0;
        for (size_t j = 0; j < na; j++) pmax = std::max(pmax, std::fabs(pull[j]));
        alpha = pmax > 0 ? pmax / 20.0 : 1.0;
    }
    double lo_u = log_dflt - 40.0, hi_u = std::log(total) + 2.0;
    std::vector<std::pair<double, double>> history;
    std::vector<double> f_prev;
    Objective obj;
    obj.model = &model; obj.y = &y; obj.w = &w; obj.dflt = dflt; obj.log_dflt = log_dflt;
    auto chi2_signal = [&]() {
        double rs = 0;
        for (size_t i = 0; i < ng; i++) if (signal[i]) { double r = fitv[i] - y[i]; rs += w[i] * r * r; }
        return rs / n_sig;
    };
    // one envelope refit: weighted NNLS at the nodes, then one common scale moved into f
    auto refit_envelope = [&]() {
        model.node_normal(f.data(), y.data(), w.data(), G, bvec);
        std::vector<double> cn = envelope_fit(G, bvec, model.Q, K, charges.front(), model.envelope());
        double s = 0;
        for (double v : cn) s += v;
        s /= model.Q;
        if (s > 0) {
            for (auto& v : cn) v /= s;
            double ls = std::log(s);
            for (auto& v : u) v += ls;
            model.set_envelope(cn);
        }
    };
    // envelope nodes at the species found so far (automatic mode): the ends of the
    // mass range plus the local maxima of the smoothed mass intensity above 1 % of
    // the largest, at least 2 % apart in mass, at most 10 (and within the cost cap)
    auto species_nodes = [&]() {
        const size_t nm = masses.size();
        std::vector<double> I(nm, 0.0);
        for (size_t j = 0; j < na; j++) {
            double cs = 0;
            for (int k = 0; k < K; k++) cs += model.env(j, k);
            I[model.act[j]] = f[j] * cs;
        }
        const long h = std::max(2L, (long)(0.0005 * 0.5 * (masses.front() + masses.back()) / step));
        std::vector<double> cum(nm + 1, 0.0);
        for (size_t i = 0; i < nm; i++) cum[i + 1] = cum[i] + I[i];
        std::vector<double> S(nm);
        double smax = 0;
        for (long i = 0; i < (long)nm; i++) {
            long a0 = std::max(0L, i - h), b0 = std::min((long)nm, i + h + 1);
            S[i] = (cum[b0] - cum[a0]) / (b0 - a0);
            smax = std::max(smax, S[i]);
        }
        std::vector<std::pair<double, double>> cand;   // (height, mass)
        for (long i = 1; i + 1 < (long)nm; i++)
            if (S[i] > S[i - 1] && S[i] >= S[i + 1] && S[i] >= 0.01 * smax) cand.push_back({S[i], masses[i]});
        std::sort(cand.begin(), cand.end(), [](const std::pair<double, double>& a, const std::pair<double, double>& b) { return a.first > b.first; });
        std::vector<double> nodes = {masses.front(), masses.back()};
        const int qcost = (int)std::floor(std::sqrt(3e9 / std::max<double>(ng, 1.0)) / std::max(K, 1));
        const int qmax = std::max(2, std::min(12, qcost));
        for (auto& c : cand) {
            if ((int)nodes.size() >= qmax) break;
            bool near = false;
            for (double m : nodes) if (std::fabs(std::log(c.second / m)) < 0.02) { near = true; break; }
            if (!near) nodes.push_back(c.second);
        }
        std::sort(nodes.begin(), nodes.end());
        return nodes;
    };
    // The search for alpha. First pass: alpha falls (slope limited steps in ln
    // alpha towards the target) until chi squared per signal point improves by
    // less than 10 % per step or lands within the target band; then the noise
    // model is re-estimated from the residuals and the envelope nodes are put at
    // the species found. Final pass: alpha falls until chi squared reaches the
    // target band (the classic chi squared = N), or stops improving (3 %: the
    // plateau P of the best fit, above the target when the peak shape is not
    // the model's), or undershoots the target (P below it: data the model
    // explains to within the noise). With a cap c > 0 the target is then raised
    // to min(target, (1 + c) P) and a secant search on the bracket climbs back
    // to it: the most entropic spectrum that fits within c of the best. The
    // default keeps the best fit (c = 0): on a synthetic two protein spectrum
    // with exact Poisson noise, chi squared = N shrank the main peak by 27 %
    // and the weaker species by 38 %, a 10 % cap by 11 % and 17 %, without
    // fewer spurious features than the best fit.
    const double ltol = std::log(1.10);       // within 10 % of the target (descent)
    const double ltol_fine = std::log(1.05);  // within 5 % of the capped target (ascent)
    double cap = 0.0;                         // the fit may be this much worse than the best in chi squared
    if (getenv("MS_MAXENT_CAP")) cap = atof(getenv("MS_MAXENT_CAP"));   // development
    int max_passes = noise_kind == 2 ? 1 : 2;
    if (getenv("MS_MAXENT_PASSES")) max_passes = std::max(1, atoi(getenv("MS_MAXENT_PASSES")));   // development
    int pass = 0;
    bool reached = false, capped = false;
    double target_cur = chi2_target, plateau = std::numeric_limits<double>::quiet_NaN();
    std::string stop_reason;
    AlphaSearch search;
    std::vector<std::pair<double, double>> pass_pts;   // (ln alpha, chi2n) of the current pass
    double la = std::log(alpha);
    ms::tick(cb, user, 0, rounds);
    for (int k = 0; k < rounds; k++) {
        alpha = std::exp(la);
        obj.alpha = alpha;
        if (have_fit) set_weights();
        for (int sweep = 0; sweep < 2; sweep++) {
            lbfgs_b(obj, u, lo_u, hi_u, iterations);
            for (size_t j = 0; j < na; j++) { u[j] = std::min(hi_u, std::max(lo_u, u[j])); f[j] = std::exp(u[j]); }
            refit_envelope();
        }
        for (size_t j = 0; j < na; j++) f[j] = std::exp(u[j]);
        model.forward(f.data(), fitv.data());
        have_fit = true;
        double chi2n = chi2_signal();
        double lc = std::log(std::max(chi2n, 1e-300) / target_cur);
        history.push_back({alpha, chi2n});
        search.add(la, lc);
        pass_pts.push_back({la, chi2n});
        if (debug) {
            fprintf(stderr, "round %d pass %d alpha %.6g chi2n %.4f (target %.3f) evals %zu bracket %.3g\n", k + 1, pass + 1, alpha, chi2n, target_cur, obj.evals, search.bracket_width());
            for (int q = 0; q < model.Q; q++) {
                fprintf(stderr, "  node %d (%.0f Da) c:", q, model.node_mass[q]);
                for (int kk = 0; kk < K; kk++) fprintf(stderr, " %.3f", model.C[(size_t)q * K + kk]);
                fprintf(stderr, "\n");
            }
        }
        ms::tick(cb, user, k + 1, rounds);
        double change = 1.0;
        if (!f_prev.empty()) {
            double num = 0, den = 0;
            for (size_t j = 0; j < na; j++) { num += std::fabs(f[j] - f_prev[j]); den += std::fabs(f[j]); }
            change = num / std::max(den, 1e-30);
        }
        f_prev = f;
        const bool enough = k + 1 >= min_rounds;
        const bool final_pass = pass + 1 >= max_passes;
        const double tol_now = capped ? ltol_fine : ltol;
        const bool on_target = std::fabs(lc) <= tol_now;
        // chi squared stopped improving although alpha fell (the plateau of the best fit)
        bool flat = false, stalled = false;
        if (search.pts.size() >= 2) {
            const auto& prev = search.pts[search.pts.size() - 2];
            const double need = final_pass ? 0.03 : 0.10;
            if (prev.first - la >= 0.5 && prev.second - lc < need && prev.second - lc > -need) flat = true;
            if (std::fabs(prev.first - la) >= 0.3 && change < 0.01) stalled = true;
        }
        if (la < std::log(history.front().first) - 25.0 || la > std::log(history.front().first) + 25.0) flat = true;
        bool narrow = search.bracket_width() < 0.05;
        if (narrow && !on_target && k + 1 < rounds) {
            // the bracket was built from fits in different states of convergence: keep the last point only
            search.reset(); search.add(la, lc);
            narrow = false;
        }
        const bool last_round = k + 1 == rounds;
        if (!final_pass) {
            // first pass: a reasonable fit is enough; then re-estimate and re-node
            const bool pass_end = (on_target && enough) || (flat && lc > -ltol) || (stalled && enough) || narrow || last_round;
            if (!pass_end) { la = on_target ? la : search.next(); continue; }
            if (debug) fprintf(stderr, "  pass %d ends: target %d flat %d stalled %d narrow %d\n", pass + 1, on_target, flat, stalled, narrow);
            if (last_round) { reached = on_target; stop_reason = on_target ? "target" : "rounds"; break; }
            std::vector<double> r(ng);
            for (size_t i = 0; i < ng; i++) r[i] = fitv[i] - y[i];
            std::string dbg;
            NoiseModel nm = noise_from_residuals(yvar, r, signal, noise_data, noise_kind, debug ? &dbg : nullptr);
            double dmax = 0;
            for (double yb : {0.0, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0}) dmax = std::max(dmax, std::fabs(0.5 * std::log(nm.var(yb) / noise.var(yb))));
            if (debug) fprintf(stderr, "  noise re-estimation: s0 %.3g -> %.3g, 1/g %.3g -> %.3g, eps %.3g -> %.3g (max change %.2f)\n%s",
                               noise.s0, nm.s0, noise.ginv, nm.ginv, noise.eps, nm.eps, dmax, dbg.c_str());
            pass++;
            bool renoded = false;
            if (p.envelope_nodes <= 0) {
                std::vector<double> nodes = species_nodes();
                if (nodes != model.node_mass && nodes.size() >= 2) {
                    model.set_nodes(nodes, masses);
                    renoded = true;
                    if (debug) { fprintf(stderr, "  envelope nodes now at:"); for (double m : nodes) fprintf(stderr, " %.1f", m); fprintf(stderr, "\n"); }
                }
            }
            noise = nm;
            set_weights();
            search.reset(); pass_pts.clear();
            if (renoded) continue;   // with new nodes the fit changes a lot at the same alpha: measure it first
            // the same fit under the new weights is the first point of the final pass
            chi2n = chi2_signal();
            lc = std::log(std::max(chi2n, 1e-300) / target_cur);
            history.back().second = chi2n;
            search.add(la, lc); pass_pts.push_back({la, chi2n});
            // within the band already: probe a lower alpha, the best fit may lie well below the target
            la = std::fabs(lc) <= ltol ? la - 2.0 : search.next();
            continue;
        }
        // final pass
        if (std::isnan(plateau)) {
            // descent: until the fit no longer improves or undershoots the target
            bool found = false;
            if (flat && lc > -ltol) {
                plateau = chi2n;
                for (auto& q : pass_pts) plateau = std::min(plateau, q.second);
                found = true;
            } else if (lc < -ltol) {
                plateau = chi2n;
                found = true;
            }
            if (found) {
                const double t2 = std::min(chi2_target, (1.0 + cap) * plateau);
                if (t2 < target_cur) {
                    target_cur = t2; capped = true;
                    search.reset();
                    for (auto& q : pass_pts) search.add(q.first, std::log(std::max(q.second, 1e-300) / target_cur));
                    lc = std::log(std::max(chi2n, 1e-300) / target_cur);
                    if (debug) fprintf(stderr, "  best fit chi2n %.4f: target now %.4f\n", plateau, target_cur);
                }
                const bool ok_now = std::fabs(lc) <= (capped ? ltol_fine : ltol) || (cap <= 0 && lc <= ltol);
                if (plateau > chi2_target / (1.0 + cap) && lc > ltol && (enough || last_round)) { reached = false; stop_reason = "unattainable"; break; }   // the best fit stays above the target
                if ((ok_now && enough) || last_round || (stalled && enough)) { reached = ok_now; stop_reason = ok_now ? "target" : stalled ? "stalled" : "rounds"; break; }
                if (ok_now || plateau > chi2_target / (1.0 + cap)) continue;   // min_rounds not reached yet: once more at this alpha
                la = search.next();
                continue;
            }
            if ((on_target && enough) || last_round || (stalled && enough)) { reached = on_target; stop_reason = on_target ? "target" : stalled ? "stalled" : "rounds"; break; }
            la = on_target ? la : search.next();
            continue;
        }
        // ascent to the capped target
        if (on_target || last_round || (narrow && enough) || (stalled && enough)) {
            reached = on_target;
            stop_reason = on_target ? "target" : narrow ? "narrow" : stalled ? "stalled" : "rounds";
            break;
        }
        la = search.next();
    }
    for (size_t j = 0; j < na; j++) f[j] = std::exp(u[j]);
    model.forward(f.data(), fitv.data());
    const double chi2_final = chi2_signal();
    // ---- _to_mass_axis
    std::vector<double> cols;
    model.columns(f.data(), cols);
    double mstep = masses.size() > 1 ? masses[1] - masses[0] : 1.0;
    double mmean = 0;
    for (double m : masses) mmean += m;
    mmean /= masses.size();
    double fwhm_mass = use_log ? model.sig * 2.3548 * d * mmean : model.sig * 2.3548 * d * (double)charges.front();
    double out_step = std::min(mstep, std::max(fwhm_mass / 6.0, 1e-4));
    // the zero charge axis at the mass step (or a sixth of the peak width): a
    // resolved isotope peak needs several points to be drawn and picked; up to
    // 2.7 the axis was capped at 400000 points, so 1000 to 50000 Da every 0.01 Da
    // got a point every 0.12 Da (isotope peaks of two points, apexes and labels
    // off the drawn maxima). The cap is now 6 million points (about 150 MB of
    // working arrays); the output is thinned away from the peaks afterwards.
    const size_t max_axis = 6000000;
    std::vector<double> mgrid = arange(masses.front(), masses.back() + out_step / 2, out_step);
    if (mgrid.size() > max_axis) {
        mgrid.resize(max_axis);
        double a0 = masses.front(), b0 = masses.back();
        for (size_t i = 0; i < max_axis; i++) mgrid[i] = a0 + (b0 - a0) * (double)i / (double)(max_axis - 1);
    }
    size_t nmg = mgrid.size();
    std::vector<double> mtotal(nmg, 0.0), h(nmg), zh(K, 0.0);
    for (int k = 0; k < K; k++) {
        double z = charges[k];
        const double* col = cols.data() + (size_t)k * ng;
        parallel_range(nmg, [&](size_t a0, size_t b0) {
            for (size_t i = a0; i < b0; i++) {
                double pos = (mgrid[i] + z * a) / z;
                double xq = use_log ? std::log(std::max(pos, 1e-9)) : pos;
                h[i] = interp0(xq, coord.data(), col, ng) * model.env_at(mgrid[i], k) * scale;
            }
        });
        double hm = 0;
        for (size_t i = 0; i < nmg; i++) { mtotal[i] += h[i]; hm = std::max(hm, h[i]); }
        zh[k] = hm;
    }
    // ---- peaks
    std::vector<Peak> peaks;
    double thresh = p.peak_threshold;
    if (use_log && iso_mode == "resolved") {
        peaks = group_isotopes(mgrid.data(), mtotal.data(), nmg, thresh);
    } else {
        double window = std::max(p.peak_window, 2.0 * step);
        if (use_log) window = std::max(window, 0.5 * model.sig * 2.3548 * d * mmean);
        peaks = pick_mass_peaks(mgrid.data(), mtotal.data(), nmg, window, thresh);
    }
    // ---- r2
    {
        double mean = 0;
        for (double v : y) mean += v;
        mean /= ng;
        double ss = 0, sr = 0;
        for (size_t i = 0; i < ng; i++) { ss += (y[i] - mean) * (y[i] - mean); sr += (y[i] - fitv[i]) * (y[i] - fitv[i]); }
        out.r2 = ss > 0 ? 1.0 - sr / ss : 0.0;
    }
    // ---- notes
    {
        char b[1400];
        std::string nz;
        if (noise_kind == 2) {
            snprintf(b, sizeof b, "noise: one sigma of %.2g %% of the base peak", 100.0 * noise.s0);
            nz = b;
        } else {
            snprintf(b, sizeof b, "noise: floor %.2g %% of the base peak", 100.0 * noise.s0);
            nz = b;
            if (noise.ginv > 0) { snprintf(b, sizeof b, ", counting noise as for %s counts at the base peak", fmt_int(noise.F / noise.ginv).c_str()); nz += b; }
            if (noise.eps > 0) { snprintf(b, sizeof b, ", peak shape mismatch %.1f %%", 100.0 * noise.eps); nz += b; }
            if (!noise_data.from_data) nz += " (floor measured on the gridded data)";
        }
        std::string env;
        if (model.Q <= 1) env = "charge envelope shared by all masses";
        else {
            env = "charge envelope at " + std::to_string(model.Q) + " nodes (";
            for (int q = 0; q < model.Q; q++) { snprintf(b, sizeof b, "%s%.0f", q ? ", " : "", model.node_mass[q]); env += b; }
            env += " Da)";
        }
        std::string tgt;
        if (capped && reached && cap > 0) snprintf(b, sizeof b, "best fit %.2f, kept within %.0f %% of it", plateau, 100.0 * cap);
        else if (!std::isnan(plateau) && plateau < chi2_target / 1.1 && reached) snprintf(b, sizeof b, "best fit, below the target %g", chi2_target);
        else if (reached) snprintf(b, sizeof b, "target %g reached", chi2_target);
        else if (stop_reason == "unattainable") snprintf(b, sizeof b, "target %g not reachable with this peak shape", chi2_target);
        else snprintf(b, sizeof b, "target %g not reached", chi2_target);
        tgt = b;
        snprintf(b, sizeof b, "Maximum entropy: charge %d to %d, mass %g to %g Da every %g Da, %s; chi squared per point %.2f after %d round%s (%s); %s; %s",
                 charges.front(), charges.back(), masses.front(), masses.back(), step, note_shape.c_str(),
                 chi2_final, (int)history.size(), history.size() > 1 ? "s" : "", tgt.c_str(), nz.c_str(), env.c_str());
        out.notes = b;
        if (!summary.empty()) out.notes += "; " + summary;
    }
    out.chi2 = chi2_final;
    out.rounds = (int)history.size();
    // ---- thinned mass axis
    {
        size_t limit_pts = 120000;
        double mx = *std::max_element(mtotal.begin(), mtotal.end());
        if (nmg <= limit_pts || mx <= 0) {
            out.mass_x = mgrid; out.mass_y = mtotal;
        } else {
            std::vector<char> keepm(nmg, 0);
            // full sampling where the spectrum is visible (0.1 % of the tallest peak is below half a
            // pixel of a linear plot); every kk-th point elsewhere (at most 120000 of them)
            for (size_t i = 0; i < nmg; i++) keepm[i] = mtotal[i] >= 1e-3 * mx;
            dilate(keepm, std::max(2L, (long)(nmg / 20000)));
            size_t kk = (size_t)std::ceil((double)nmg / (double)limit_pts);
            for (size_t i = 0; i < nmg; i += kk) keepm[i] = 1;
            keepm[0] = keepm[nmg - 1] = 1;
            for (size_t i = 0; i < nmg; i++) if (keepm[i]) { out.mass_x.push_back(mgrid[i]); out.mass_y.push_back(mtotal[i]); }
        }
    }
    // ---- fit with NaN gaps
    {
        double nan = std::numeric_limits<double>::quiet_NaN();
        for (size_t i = 0; i < ng; i++) {
            if (i > 0 && coord[i] - coord[i - 1] > 1.5 * d) { out.fit_x.push_back(nan); out.fit_y.push_back(nan); }
            out.fit_x.push_back(grid[i]); out.fit_y.push_back(fitv[i] * scale);
        }
    }
    for (int k = 0; k < K; k++) { out.z.push_back(charges[k]); out.zdist.push_back(zh[k]); }
    for (auto& q : peaks) {
        out.pk_mass.push_back(q.mass); out.pk_height.push_back(q.height); out.pk_area.push_back(q.area);
        out.pk_apex.push_back(q.apex); out.pk_niso.push_back(q.n_iso);
    }
}

}  // namespace

extern "C" {

MS_API int ms_maxent(const double* mz, const double* it, long n, const ms_maxent_params* p,
                     ms_maxent_result** out, ms_progress_fn cb, void* user) {
    if (out) *out = nullptr;
    return ms::guarded([&] {
        if (!mz || !it || n < 0 || !p || !out) throw std::invalid_argument("ms_maxent: bad arguments");
        std::unique_ptr<ms_maxent_result> r(new ms_maxent_result);
        run_maxent(mz, it, n, *p, *r, cb, user);
        *out = r.release();
    });
}

MS_API void ms_maxent_free(ms_maxent_result* r) { delete r; }

MS_API long ms_maxent_mass(const ms_maxent_result* r, const double** mass, const double** it) {
    if (!r) return 0;
    if (mass) *mass = r->mass_x.data();
    if (it) *it = r->mass_y.data();
    return (long)r->mass_x.size();
}
MS_API long ms_maxent_fit(const ms_maxent_result* r, const double** mz, const double** it) {
    if (!r) return 0;
    if (mz) *mz = r->fit_x.data();
    if (it) *it = r->fit_y.data();
    return (long)r->fit_x.size();
}
MS_API long ms_maxent_zdist(const ms_maxent_result* r, const double** z, const double** it) {
    if (!r) return 0;
    if (z) *z = r->z.data();
    if (it) *it = r->zdist.data();
    return (long)r->z.size();
}
MS_API long ms_maxent_peaks(const ms_maxent_result* r, const double** mass, const double** height,
                            const double** area, const double** apex, const int** n_iso) {
    if (!r) return 0;
    if (mass) *mass = r->pk_mass.data();
    if (height) *height = r->pk_height.data();
    if (area) *area = r->pk_area.data();
    if (apex) *apex = r->pk_apex.data();
    if (n_iso) *n_iso = r->pk_niso.data();
    return (long)r->pk_mass.size();
}
MS_API double ms_maxent_r2(const ms_maxent_result* r) { return r ? r->r2 : 0.0; }
MS_API double ms_maxent_chi2(const ms_maxent_result* r) { return r ? r->chi2 : 0.0; }
MS_API int    ms_maxent_rounds(const ms_maxent_result* r) { return r ? r->rounds : 0; }
MS_API const char* ms_maxent_notes(const ms_maxent_result* r) { return r ? r->notes.c_str() : ""; }

}

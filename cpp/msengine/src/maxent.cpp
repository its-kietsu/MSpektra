// Maximum entropy charge state deconvolution (written for MSpektra from the
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
//   fit       u = ln f is optimised with an L-BFGS (compact form, m = 10,
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
#include <chrono>
#include <atomic>
#if defined(__x86_64__) || defined(__i386__)
#include <xmmintrin.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

using ms::e2::Peak;

namespace {

// ============================================================ small helpers
inline int omp_threads() { return ms::threads(); }

// time spent in the parts of a run (MS_PROFILE set: printed to stderr at its end)
enum { T_BUILD, T_FORWARD, T_ADJOINT, T_OBJECTIVE, T_LBFGS, T_NORMAL, T_ENVFIT, T_EXP, T_DIR, T_PAIR, T_OUTPUT, T_TOTAL, T_N };
double prof_t[T_N];
long prof_n[T_N];  // calls
const char* prof_name[T_N] = {"model build", "forward", "adjoint", "objective (all)", "lbfgs (all)",
                              "node normal", "envelope fit", "exp (objective)", "lbfgs direction",
                              "lbfgs pair", "output", "total"};
struct Timed {
    int k;
    std::chrono::steady_clock::time_point s;
    explicit Timed(int k) : k(k), s(std::chrono::steady_clock::now()) {}
    ~Timed() { prof_t[k] += std::chrono::duration<double>(std::chrono::steady_clock::now() - s).count(); prof_n[k]++; }
};

// The floating point settings of the calling thread, given to the worker threads: Python's thread runs the x87 unit
// in double precision (control word 027f), new threads in extended precision (037f), and msvcrt's exp uses it, so
// an exp on a worker thread gave other last bits and the results changed with the thread count.
struct FpState {
#if defined(__x86_64__) || defined(__i386__)
    unsigned short cw; unsigned csr;
    FpState() { __asm__ volatile("fnstcw %0" : "=m"(cw)); csr = _mm_getcsr() & ~0x3fu; }
    void apply() const { __asm__ volatile("fldcw %0" : : "m"(cw)); _mm_setcsr(csr); }
#else
    void apply() const {}
#endif
};

// exp of fdlibm (Sun Microsystems, e_exp.c; error below 1 ulp) in plain double arithmetic, for the masses
// f = exp(u): msvcrt's exp runs on the x87 unit, 29 ns a value against 8 ns, and was off the correctly
// rounded value more often (17 % against 10 % of values, never by more than 1 ulp)
inline double fexp(double x) {
    const double ln2hi = 6.93147180369123816490e-01, ln2lo = 1.90821492927058770002e-10,
                 invln2 = 1.44269504088896338700e+00, P1 = 1.66666666666666019037e-01,
                 P2 = -2.77777777770155933842e-03, P3 = 6.61375632143793436117e-05,
                 P4 = -1.65339022054652515390e-06, P5 = 4.13813679705723846039e-08;
    if (!(x == x)) return x;
    if (x > 7.09782712893383973096e+02) return HUGE_VAL;
    if (x < -7.45133219101941108420e+02) return 0.0;
    double hi = 0, lo = 0;
    int k = 0;
    const double ax = std::fabs(x);
    if (ax > 0.34657359027997264) {   // x = k ln 2 + r, |r| <= ln 2 / 2
        if (ax < 1.0397207708399179) { k = x > 0 ? 1 : -1; hi = x - k * ln2hi; lo = k * ln2lo; }
        else { k = (int)(invln2 * x + (x > 0 ? 0.5 : -0.5)); const double t = k; hi = x - t * ln2hi; lo = t * ln2lo; }
        x = hi - lo;
    } else if (ax < 3.725290298461914e-09) {
        return 1.0 + x;
    }
    const double t = x * x;
    const double c = x - t * (P1 + t * (P2 + t * (P3 + t * (P4 + t * P5))));
    if (k == 0) return 1.0 - ((x * c) / (c - 2.0) - x);
    double y = 1.0 - ((lo - (x * c) / (2.0 - c)) - hi);
    uint64_t b;
    std::memcpy(&b, &y, 8);
    if (k >= -1021) {   // times 2^k
        b += (uint64_t)(int64_t)k << 52;
        std::memcpy(&y, &b, 8);
        return y;
    }
    b += (uint64_t)(int64_t)(k + 1000) << 52;
    std::memcpy(&y, &b, 8);
    return y * 9.33263618503218878990e-302;   // 2^-1000
}

template <class F> void parallel_range(size_t n, F&& f) {
    // f(i0, i1) over disjoint contiguous blocks
    int nt = std::max(1, omp_threads());
    if (n < 4096 || nt == 1) { f((size_t)0, n); return; }
#ifdef _OPENMP
    const FpState fp;
#pragma omp parallel num_threads(nt)
    {
        int t = omp_get_thread_num(), T = omp_get_num_threads();
        fp.apply();
        size_t a = n * (size_t)t / T, b = n * (size_t)(t + 1) / T;
        if (b > a) f(a, b);
    }
#else
    f((size_t)0, n);
#endif
}

// Sum of f(i) over [0, n) on every core: fixed blocks of PBLOCK points are summed in parallel and the block sums
// added in block order, so the result does not depend on the thread count (one block: the plain serial sum).
// Below PMIN points the blocks are summed on one core (the same sum; starting threads would cost more).
constexpr size_t PBLOCK = 8192;
constexpr size_t PMIN = 65536;
template <class F> double psum(size_t n, F&& f) {
    const size_t nb = (n + PBLOCK - 1) / PBLOCK;
    if (nb <= 1) { double s = 0; for (size_t i = 0; i < n; i++) s += f(i); return s; }
    std::vector<double> part(nb, 0.0);
#ifdef _OPENMP
    const FpState fp;
#pragma omp parallel for schedule(static) num_threads(std::max(1, omp_threads())) if (n >= PMIN)
#endif
    for (long b = 0; b < (long)nb; b++) {
#ifdef _OPENMP
        fp.apply();
#endif
        const size_t i0 = (size_t)b * PBLOCK, i1 = std::min(n, i0 + PBLOCK);
        double s = 0;
        for (size_t i = i0; i < i1; i++) s += f(i);
        part[b] = s;
    }
    double s = 0;
    for (double v : part) s += v;
    return s;
}

double pdot(const double* a, const double* b, size_t n) { return psum(n, [&](size_t i) { return a[i] * b[i]; }); }

// m sums over [0, n) in one pass, in the fixed blocks of psum: f(i0, i1, acc) puts the m sums over the block
// [i0, i1) into acc (each summed in index order); the block sums are added in block order
template <class F> void psums(size_t n, int m, double* out, F&& f) {
    const size_t nb = std::max<size_t>(1, (n + PBLOCK - 1) / PBLOCK);
    std::vector<double> part(nb * (size_t)m, 0.0);
#ifdef _OPENMP
    const FpState fp;
#pragma omp parallel for schedule(static) num_threads(std::max(1, omp_threads())) if (n >= PMIN)
#endif
    for (long b = 0; b < (long)nb; b++) {
#ifdef _OPENMP
        fp.apply();
#endif
        const size_t i0 = (size_t)b * PBLOCK, i1 = std::min(n, i0 + PBLOCK);
        f(i0, i1, part.data() + (size_t)b * m);
    }
    for (int c = 0; c < m; c++) {
        double s = 0;
        for (size_t b = 0; b < nb; b++) s += part[b * (size_t)m + c];
        out[c] = s;
    }
}

// largest f(i) over [0, n) (start: the value for an empty range); any order gives the same maximum
template <class F> double pmax(size_t n, double start, F&& f) {
    const size_t nb = (n + PBLOCK - 1) / PBLOCK;
    std::vector<double> part(std::max<size_t>(nb, 1), start);
#ifdef _OPENMP
    const FpState fp;
#pragma omp parallel for schedule(static) num_threads(std::max(1, omp_threads())) if (n >= PMIN)
#endif
    for (long b = 0; b < (long)nb; b++) {
#ifdef _OPENMP
        fp.apply();
#endif
        const size_t i0 = (size_t)b * PBLOCK, i1 = std::min(n, i0 + PBLOCK);
        double m = start;
        for (size_t i = i0; i < i1; i++) m = std::max(m, f(i));
        part[b] = m;
    }
    double m = start;
    for (double v : part) m = std::max(m, v);
    return m;
}

// f(i0, i1) over [0, n) on every core (element by element work: the result is the same in any order); one core
// below PMIN points
template <class F> void pfor(size_t n, F&& f) {
    if (n < PMIN) { f((size_t)0, n); return; }
    parallel_range(n, f);
}

// f(b) for b in [0, nb) on every core, a block at a time to whichever core is free
template <class F> void pblocks(size_t nb, F&& f) {
#ifdef _OPENMP
    const FpState fp;
#pragma omp parallel for schedule(dynamic, 1) num_threads(std::max(1, omp_threads())) if (nb > 1)
#endif
    for (long b = 0; b < (long)nb; b++) {
#ifdef _OPENMP
        fp.apply();
#endif
        f((size_t)b);
    }
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
// A (grid x active masses): per mass and charge, the two grid points around
// (m + z a)/z with the weights of a linear interpolation. Stored by mass:
// each active mass has one entry per charge from its lowest to its highest
// charge inside the grid, the lower grid point (int32) and the weight of the
// upper one (float32; negative across a gap in the grid), 8 bytes. Both
// products read a 5 byte copy of the entries once, in order (bound by the memory
// speed). At a given charge the lower grid point never falls as the mass
// rises (for the products of one charge: rows_of). The
// envelope weight c(z | M) is applied on the fly from the node table, so a
// new envelope costs nothing. Per active mass: the lower envelope node
// (uint16) and the weight of the upper one (float32).
struct Model {
    size_t n = 0;        // grid points
    int K = 0;           // charges
    size_t na = 0;       // active masses
    std::vector<int> act;        // active mass -> mass index
    double sig = 1;
    Blur blur;
    struct Entry { int32_t i; float wh; };   // lower grid point, weight of the upper one (< 0: across a gap)
    std::vector<Entry> ent;
    std::vector<size_t> m_off;               // per active mass: entries ent[m_off[j] .. m_off[j + 1]),
    std::vector<uint16_t> m_k0;              // for the charges m_k0[j] on
    std::vector<size_t> e_lo, e_hi;          // per charge: the active masses lo .. hi - 1 it puts inside the grid
    // A fixed number of mass blocks with about equal entries (a fixed count, so that no sum depends on the thread
    // count): masses b_j[b] .. b_j[b + 1] - 1 reach grid rows b_r0[b] .. b_r1[b] - 1; their part of the forward
    // product is kept at part[b_p[b]]
    static constexpr size_t NBLK = 60;
    std::vector<size_t> b_j, b_r0, b_r1, b_p;
    std::vector<double> part;
    // EXPERIMENT (with MS_MEMSYS): the parts of the forward product by block and charge, over the rows that charge
    // reaches from the block (a block of masses lands in one narrow window per charge): window (b, k) is rows
    // w_r0 .. w_r1 - 1 at part_w[w_off]
    bool windowed = false;
    std::vector<size_t> w_r0, w_r1, w_off;
    std::vector<double> part_w;
    void make_windows() {
        const size_t nb = b_r0.size();
        w_r0.assign(nb * K, 0); w_r1.assign(nb * K, 0); w_off.assign(nb * K + 1, 0);
        pblocks(nb, [&](size_t b) {
            for (size_t j = b_j[b]; j < b_j[b + 1]; j++)
                for (size_t p = m_off[j]; p < m_off[j + 1]; p++) {
                    if (ent[p].wh < 0) continue;
                    const size_t q = b * K + m_k0[j] + (p - m_off[j]), i = (size_t)ent[p].i;
                    if (w_r1[q] == 0) { w_r0[q] = i; w_r1[q] = i + 2; }
                    else { w_r0[q] = std::min(w_r0[q], i); w_r1[q] = std::max(w_r1[q], i + 2); }
                }
        });
        for (size_t q = 0; q < nb * K; q++) w_off[q + 1] = w_off[q] + (w_r1[q] - w_r0[q]);
        part_w.assign(w_off[nb * K], 0.0);
        windowed = true;
    }
    // The entries again for the two products, 5 bytes each (they are bound by the memory speed): the weight and
    // the step in grid rows from the entry of the same charge at the mass before (255: the row is the next one
    // of x_row, for the first entry of a charge in a block and for larger steps). b_row[b K + k]: the row of
    // charge k before block b; b_x[b]: block b's first row in x_row
    std::vector<float> c_w;
    std::vector<uint8_t> c_d;
    std::vector<int32_t> x_row, b_row;
    std::vector<size_t> b_x;
    size_t nnz = 0;                          // non zero weights (for the profile)
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

    // bytes the model takes for ne entries
    static double bytes_for(size_t ne, size_t n, size_t na) {
        return (double)ne * (sizeof(Entry) + 5) + 8.0 * NBLK * n + 8.0 * n + 16.0 * na;
    }
    // the entry of mass j at charge k (k one of its charges)
    const Entry& entry(size_t j, int k) const { return ent[m_off[j] + (size_t)(k - m_k0[j])]; }
    // the masses j0 .. j1 - 1 whose lower grid point at charge k is in [lo, hi) (it never falls as the mass rises)
    std::pair<size_t, size_t> rows_of(int k, long lo, long hi) const {
        auto first = [&](long v) {
            size_t a = e_lo[k], b = e_hi[k];
            while (a < b) { const size_t m = a + (b - a) / 2; if (entry(m, k).i < v) a = m + 1; else b = m; }
            return a;
        };
        return {first(lo), first(hi)};
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
        Timed tm(T_BUILD);
        n = coord.size();
        K = (int)charges.size();
        if (K > 65535) throw std::invalid_argument("too many charge states");
        sig = std::max(0.3, fwhm / d / 2.3548);
        blur.init(sig);
        coord_cp = coord; uselog = log; adduct = a; charges_cp = charges;
        const size_t nm = masses.size();
        std::vector<char> active(nm, 0);
        // the masses ja .. jb - 1 inside the grid at charge k (positions increase with the mass): fn(j, i0, w1, ok),
        // ok false across a gap. The lower grid point of the first mass by a search, then by a cursor (the same
        // point a cursor from the first mass finds: a part of the masses can be walked on each core)
        auto walk = [&](int k, size_t ja, size_t jb, auto&& fn) {
            const double z = charges[k];
            long i0 = -1;
            if (ja < jb) {
                double pos = (masses[ja] + z * a) / z;
                if (log) pos = std::log(std::max(pos, 1e-9));
                i0 = (long)(std::upper_bound(coord.begin(), coord.end(), pos) - coord.begin()) - 1;
            }
            for (size_t j = ja; j < jb; j++) {
                double pos = (masses[j] + z * a) / z;
                if (log) pos = std::log(std::max(pos, 1e-9));
                while (i0 + 1 < (long)n && coord[i0 + 1] <= pos) i0++;
                if (i0 < 0 || i0 >= (long)n - 1) continue;
                if (!(coord[i0 + 1] - coord[i0] < 1.5 * d)) { fn(j, (size_t)i0, 0.0, false); continue; }
                const double w1 = std::min(1.0, std::max(0.0, (pos - coord[i0]) / d));
                fn(j, (size_t)i0, w1, true);
            }
        };
        // pass 1: the active masses, per charge the range of masses inside the grid and per mass its charges
        // (charges rise with k, so those of a mass inside the grid are one run of k)
        std::vector<size_t> m_lo(K, 0), m_hi(K, 0);
        std::vector<uint16_t> kmin(nm, 0xffff), kmax(nm, 0);
        nnz = 0;
        size_t ne_max = 0;
        for (int k = 0; k < K; k++) {   // each core a part of the masses (each mass is written by one)
            std::atomic<size_t> lo_k{SIZE_MAX}, hi_k{0}, nz{0};
            parallel_range(nm, [&](size_t ja, size_t jb) {
                size_t first = SIZE_MAX, last = 0, cnt = 0;
                walk(k, ja, jb, [&](size_t j, size_t, double, bool ok) {
                    if (first == SIZE_MAX) first = j;
                    last = j + 1;
                    if (kmin[j] == 0xffff) kmin[j] = (uint16_t)k;
                    kmax[j] = (uint16_t)k;
                    if (ok) { active[j] = 1; cnt += 2; }
                });
                if (first != SIZE_MAX) {
                    for (size_t v = lo_k.load(); first < v && !lo_k.compare_exchange_weak(v, first);) {}
                    for (size_t v = hi_k.load(); last > v && !hi_k.compare_exchange_weak(v, last);) {}
                }
                nz += cnt;
            });
            if (lo_k.load() != SIZE_MAX) { m_lo[k] = lo_k.load(); m_hi[k] = hi_k.load(); }
            nnz += nz.load();
            ne_max += m_hi[k] - m_lo[k];
        }
        std::vector<int> colmap(nm, -1);
        std::vector<size_t> before(nm + 1, 0);   // active masses below mass j
        act.clear();
        for (size_t j = 0; j < nm; j++) {
            before[j] = act.size();
            if (active[j]) { colmap[j] = (int)act.size(); act.push_back((int)j); }
        }
        before[nm] = act.size();
        na = act.size();
        {
            const double need = bytes_for(ne_max, n, na) + 230.0 * na + 16.0 * K * n;
            if (need > budget) {
                char b[600];
                snprintf(b, sizeof b, "This deconvolution needs about %s of memory (%zu masses x %d charges on %zu grid points); this computer can give it %s. Use a larger mass step or narrower mass and charge ranges (isotope peaks need 0.01 to 0.02 Da, envelopes 0.5 to 1 Da).",
                         ms::fmt_bytes(need).c_str(), na, K, n, ms::fmt_bytes(budget).c_str());
                throw std::invalid_argument(b);
            }
        }
        const size_t nac = std::max<size_t>(na, 1);
        e_lo.assign(K, 0); e_hi.assign(K, 0);
        for (int k = 0; k < K; k++) { e_lo[k] = before[m_lo[k]]; e_hi[k] = std::max(e_lo[k], before[m_hi[k]]); }
        m_off.assign(nac + 1, 0); m_k0.assign(nac, 0);
        for (size_t cj = 0; cj < na; cj++) {
            const size_t j = (size_t)act[cj];
            m_k0[cj] = kmin[j];
            m_off[cj + 1] = m_off[cj] + (size_t)(kmax[j] - kmin[j] + 1);
        }
        ent.assign(m_off[na], Entry{0, -1.0f});
        // pass 2: fill (every charge of an active mass is met)
        for (int k = 0; k < K; k++)
            parallel_range(nm, [&](size_t ja, size_t jb) {
                walk(k, ja, jb, [&](size_t j, size_t i0, double w1, bool ok) {
                    if (!active[j]) return;
                    Entry& e = ent[m_off[colmap[j]] + (size_t)(k - kmin[j])];
                    e.i = (int32_t)i0;
                    e.wh = ok ? (float)w1 : -1.0f;
                });
            });
        // the mass blocks: about equal entries each, and the grid rows each reaches
        const size_t ne = ent.size();
        const size_t nb = std::max<size_t>(1, std::min<size_t>(NBLK, ne / 262144));
        b_j.assign(nb + 1, na);
        b_j[0] = 0;
        for (size_t b = 1; b < nb; b++)
            b_j[b] = (size_t)(std::lower_bound(m_off.begin(), m_off.begin() + na, ne * b / nb) - m_off.begin());
        b_r0.assign(nb, 0); b_r1.assign(nb, 0); b_p.assign(nb + 1, 0);
        for (size_t b = 0; b < nb; b++) {
            long r0 = (long)n, r1 = 0;
            for (size_t p = m_off[b_j[b]]; p < m_off[b_j[b + 1]]; p++)
                if (ent[p].wh >= 0) { r0 = std::min(r0, (long)ent[p].i); r1 = std::max(r1, (long)ent[p].i + 2); }
            if (r1 <= r0) r0 = r1 = 0;
            b_r0[b] = (size_t)r0; b_r1[b] = (size_t)r1;
            b_p[b + 1] = b_p[b] + (b_r1[b] - b_r0[b]);
        }
        part.assign(b_p[nb], 0.0);
        // the 5 byte entries
        c_w.resize(ne); c_d.resize(ne); x_row.clear(); b_row.assign(nb * (size_t)K, 0); b_x.assign(nb + 1, 0);
        {
            std::vector<int64_t> prev(K, -1);   // the last row of each charge (-1: none yet)
            for (size_t b = 0; b < nb; b++) {
                for (int k = 0; k < K; k++) b_row[b * K + k] = (int32_t)std::max<int64_t>(prev[k], 0);
                b_x[b] = x_row.size();
                for (size_t j = b_j[b]; j < b_j[b + 1]; j++)
                    for (size_t p = m_off[j]; p < m_off[j + 1]; p++) {
                        const int k = m_k0[j] + (int)(p - m_off[j]);
                        const int64_t i = ent[p].i, dd = i - prev[k];
                        c_w[p] = ent[p].wh;
                        if (prev[k] < 0 || dd < 0 || dd >= 255) { c_d[p] = 255; x_row.push_back((int32_t)i); }
                        else c_d[p] = (uint8_t)dd;
                        prev[k] = i;
                    }
            }
            b_x[nb] = x_row.size();
        }
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
        Timed tm(T_FORWARD);
        if (windowed) { forward_windowed(f, out); return; }
        // each mass block its own part of the grid rows (masses in order, each by charge); then the parts are added
        // in block order
        const size_t nb = b_r0.size();
        pblocks(nb, [&](size_t b) {
            const size_t r0 = b_r0[b];
            double* P = part.data() + b_p[b];
            std::fill(P, P + (b_r1[b] - r0), 0.0);
            std::vector<int32_t> cur(b_row.begin() + b * K, b_row.begin() + (b + 1) * K);
            const int32_t* xr = x_row.data() + b_x[b];
            for (size_t j = b_j[b]; j < b_j[b + 1]; j++) {
                const double fj = f[j], tj = nt[j];
                const double* c0 = C.data() + (size_t)nq[j] * K;
                const float* w = c_w.data() + m_off[j];
                const uint8_t* st = c_d.data() + m_off[j];
                const int k0 = m_k0[j], ne = (int)(m_off[j + 1] - m_off[j]);
                int32_t* cu = cur.data() + k0;
                for (int c = 0; c < ne; c++) {
                    const int32_t row = st[c] == 255 ? *xr++ : cu[c] + st[c];
                    cu[c] = row;
                    const double wh = w[c];
                    if (wh < 0) continue;
                    const int k = k0 + c;
                    const double ce = c0[k] + tj * (c0[K + k] - c0[k]);
                    const size_t i = (size_t)row - r0;
                    P[i] += ce * (1.0 - wh) * fj;
                    P[i + 1] += ce * wh * fj;
                }
            }
        });
        parallel_range(n, [&](size_t a, size_t b) {
            double* t = tmp.data();
            std::fill(t + a, t + b, 0.0);
            for (size_t q = 0; q < nb; q++) {
                const size_t lo = std::max(a, b_r0[q]), hi = std::min(b, b_r1[q]);
                const double* P = part.data() + b_p[q] - b_r0[q];
                for (size_t i = lo; i < hi; i++) t[i] += P[i];
            }
        });
        blur.apply(tmp.data(), out, n);
    }

    void forward_windowed(const double* f, double* out) {
        const size_t nb = b_r0.size();
        pblocks(nb, [&](size_t b) {
            std::vector<double*> Pk(K);
            for (int k = 0; k < K; k++) {
                const size_t q = b * K + k;
                double* P = part_w.data() + w_off[q];
                std::fill(P, P + (w_r1[q] - w_r0[q]), 0.0);
                Pk[k] = P - w_r0[q];
            }
            std::vector<int32_t> cur(b_row.begin() + b * K, b_row.begin() + (b + 1) * K);
            const int32_t* xr = x_row.data() + b_x[b];
            for (size_t j = b_j[b]; j < b_j[b + 1]; j++) {
                const double fj = f[j], tj = nt[j];
                const double* c0 = C.data() + (size_t)nq[j] * K;
                const float* w = c_w.data() + m_off[j];
                const uint8_t* st = c_d.data() + m_off[j];
                const int k0 = m_k0[j], ne = (int)(m_off[j + 1] - m_off[j]);
                int32_t* cu = cur.data() + k0;
                for (int c = 0; c < ne; c++) {
                    const int32_t row = st[c] == 255 ? *xr++ : cu[c] + st[c];
                    cu[c] = row;
                    const double wh = w[c];
                    if (wh < 0) continue;
                    const int k = k0 + c;
                    const double ce = c0[k] + tj * (c0[K + k] - c0[k]);
                    double* P = Pk[k];
                    P[row] += ce * (1.0 - wh) * fj;
                    P[row + 1] += ce * wh * fj;
                }
            }
        });
        // the windows added up in a fixed order (block, then charge)
        const size_t nw = nb * K;
        parallel_range(n, [&](size_t a, size_t b) {
            double* t = tmp.data();
            std::fill(t + a, t + b, 0.0);
            for (size_t q = 0; q < nw; q++) {
                const size_t lo = std::max(a, w_r0[q]), hi = std::min(b, w_r1[q]);
                if (hi <= lo) continue;
                const double* P = part_w.data() + w_off[q] - w_r0[q];
                for (size_t i = lo; i < hi; i++) t[i] += P[i];
            }
        });
        blur.apply(tmp.data(), out, n);
    }

    // gradient with respect to the masses of a (weighted) residual r on the grid
    void adjoint_f(const double* r, double* out) { adjoint_f(r, out, [](size_t, double) {}); }
    // ... and post(j, out[j]) for every mass, in the same pass
    template <class Post> void adjoint_f(const double* r, double* out, Post&& post) {
        Timed tm(T_ADJOINT);
        blur.apply(r, tmp.data(), n);
        const double* g = tmp.data();
        // a mass adds its terms by charge
        pblocks(b_r0.size(), [&](size_t b) {
            std::vector<int32_t> cur(b_row.begin() + b * K, b_row.begin() + (b + 1) * K);
            const int32_t* xr = x_row.data() + b_x[b];
            for (size_t j = b_j[b]; j < b_j[b + 1]; j++) {
                const double tj = nt[j];
                const double* c0 = C.data() + (size_t)nq[j] * K;
                const float* w = c_w.data() + m_off[j];
                const uint8_t* st = c_d.data() + m_off[j];
                const int k0 = m_k0[j], ne = (int)(m_off[j + 1] - m_off[j]);
                int32_t* cu = cur.data() + k0;
                double s = 0, s2 = 0;   // two sums (lower and upper grid point): half the chain of additions
                for (int c = 0; c < ne; c++) {
                    const int32_t row = st[c] == 255 ? *xr++ : cu[c] + st[c];
                    cu[c] = row;
                    const double wh = w[c];
                    if (wh < 0) continue;
                    const int k = k0 + c;
                    const double ce = c0[k] + tj * (c0[K + k] - c0[k]);
                    const size_t i = (size_t)row;
                    s += ce * (1.0 - wh) * g[i];
                    s2 += ce * wh * g[i + 1];
                }
                s += s2;
                out[j] = s;
                post(j, s);
            }
        });
    }

    // unweighted m/z spectrum of each charge state: K x n (row major), for the output
    void columns(const double* f, std::vector<double>& out) {
        out.assign((size_t)K * n, 0.0);
        std::vector<double> raw((size_t)K * n, 0.0);
        // each charge its own row of raw; masses in order
        pblocks((size_t)K, [&](size_t k) {
            double* rk = raw.data() + k * n;
            for (size_t j = e_lo[k]; j < e_hi[k]; j++) {
                const Entry& e = entry(j, (int)k);
                const double wh = e.wh;
                if (wh < 0) continue;
                rk[e.i] += (1.0 - wh) * f[j];
                rk[e.i + 1] += wh * f[j];
            }
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
        Timed tm(T_NORMAL);
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
            // each charge and stretch of rows its own cells; masses in order
            const size_t NS = 8;
            pblocks((size_t)K * NS, [&](size_t task) {
                const int k = (int)(task / NS);
                const size_t s = task % NS;
                const long ra = (long)(ea + L * s / NS), rb = (long)(ea + L * (s + 1) / NS);
                if (rb <= ra) return;
                const auto [j0, j1] = rows_of(k, ra - 1, rb);
                for (size_t j = j0; j < j1; j++) {
                    const Entry& e = entry(j, k);
                    const double wh = e.wh;
                    if (wh < 0) continue;
                    const long i = e.i;
                    const double t = nt[j];
                    const size_t q = nq[j];
                    if (i >= ra) {
                        const double v = (1.0 - wh) * f[j];
                        const size_t ii = (size_t)i - ea;
                        raw[(q * K + k) * L + ii] += v * (1.0 - t);
                        if (t > 0) raw[((q + 1) * K + k) * L + ii] += v * t;
                    }
                    if (i + 1 < rb) {
                        const double v = wh * f[j];
                        const size_t ii = (size_t)i + 1 - ea;
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
    // EXPERIMENT (MS_ME_RATIO2): a charge state's two neighbours *together* carry at least `ratio` of it,
    // lam (max(0, r x_z - x_z-1 - x_z+1))^2: still no charge state alone, but a narrow native envelope
    // (BSA: 36 %, 100 %, 24 %, 0.2 %) is allowed; the pairwise form forces 1/3 on every neighbour
    static const bool ratio2 = getenv("MS_ME_RATIO2") != nullptr;
    const double step = ratio2 ? 1.0 / (1.05 * lam + 6.0 * (2.0 + ratio * ratio) * lpen)
                               : 1.0 / (1.05 * lam + 4.0 * (1.0 + ratio * ratio) * lpen);
    auto add_pen_grad = [&](const std::vector<double>& z, std::vector<double>& grad) {
        if (!use_pen) return;
        if (ratio2) {
            for (int q = 0; q < Q; q++)
                for (int k = 0; k < K; k++) {
                    const size_t i = (size_t)q * K + k;
                    const double l = pen[(size_t)q * K + std::min(k, K - 2)];
                    if (!(l > 0)) continue;
                    const double nb = (k > 0 ? z[i - 1] : 0.0) + (k + 1 < K ? z[i + 1] : 0.0);
                    const double d = std::max(0.0, ratio * z[i] - nb);
                    if (d <= 0) continue;
                    grad[i] += 2.0 * l * ratio * d;
                    if (k > 0) grad[i - 1] -= 2.0 * l * d;
                    if (k + 1 < K) grad[i + 1] -= 2.0 * l * d;
                }
            return;
        }
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
// EXPERIMENT (MS_ME_LAMLOCAL=w): the coupling of two neighbouring envelope nodes falls with their mass ratio,
// lam_node x exp(-(ln(m1/m0) / w)^2): nodes of one protein's cluster (variants, adducts: a few % apart) stay tied,
// unrelated proteins (29 against 148 kDa) get their own charge distributions. The caller sets the node masses.
// The engine of a run (since 4.36). noise_model 0 (automatic): the engine of 4.36, found on the experiment/memsys
// branch and frozen for the paper (commit 48e81f5): area output (each charge state divided by z, so species are not
// weighted by their charge), the noise measured once on the data and kept (no re-estimation from the residuals),
// chi2 over the points with signal (3 sd of their own noise), neighbour rule 0.2, charge envelopes coupled only
// within a few % of mass, envelope nodes seeded where signal is unexplained, final pass by the Skilling-Bryan path
// (hybrid solver). noise_model 1 or 4: the engine as in 4.35 (every one of these off); 2 and 3: the comparison noise
// models of the 4.35 engine. The MS_* environment variables still override each part (development, benchmarks).
struct MaxentOptions {
    bool legacy = true;
    int area = 0;              // 0 heights (4.35), 1 divided by z (signal), 2 by z^2 (ions on image-current detectors)
    bool memsys = false, hybrid = false;
    bool noise2 = false, fixnoise = false, snrmask = false, seednodes = false;
    bool seedharm = false;     // lighter nodes not on charge states of a stronger species (harmonic ghosts)
    double ratio = -1;         // minimum neighbour ratio; < 0: the default of envelope_fit (0.33)
    double lamlocal = 0, lamlocal0 = 0, lamend = 0;
};
static bool opt_flag(const char* name, bool dflt) {
    const char* v = getenv(name);
    return v ? std::strcmp(v, "0") != 0 : dflt;
}
static double opt_num(const char* name, double dflt) {
    const char* v = getenv(name);
    return v ? atof(v) : dflt;
}
static MaxentOptions maxent_options(bool legacy) {
    MaxentOptions o;
    o.legacy = legacy;
    const bool n = !legacy;
    o.area = (int)opt_num("MS_ME_AREA", n ? 1 : 0);
    o.memsys = opt_flag("MS_MEMSYS", n);
    o.hybrid = o.memsys && opt_flag("MS_MEMSYS_SAFE", n);
    o.noise2 = opt_flag("MS_ME_NOISE2", n);
    o.fixnoise = opt_flag("MS_ME_FIXNOISE", n);
    o.snrmask = opt_flag("MS_ME_SNRMASK", n);
    o.seednodes = opt_flag("MS_ME_SEEDNODES", n);
    o.seedharm = opt_flag("MS_ME_SEEDHARM", n);
    o.ratio = opt_num("MS_ME_RATIO", n ? 0.2 : -1.0);
    o.lamlocal = opt_num("MS_ME_LAMLOCAL", n ? 0.03 : 0.0);
    o.lamlocal0 = opt_num("MS_ME_LAMLOCAL0", n ? 0.03 : 0.0);
    o.lamend = opt_num("MS_ME_LAMEND", n ? 1.0 : 0.0);
    return o;
}
static thread_local const MaxentOptions* g_opt = nullptr;   // set by run_maxent for the functions it calls
static const MaxentOptions& opt() {
    static const MaxentOptions legacy_default;   // (a call outside a run: as 4.35)
    return g_opt ? *g_opt : legacy_default;
}
static thread_local const std::vector<double>* g_env_node_mass = nullptr;
std::vector<double> envelope_fit(std::vector<double> G, const std::vector<double>& b, int Q, int K, int zfirst,
                                 const std::vector<double>& warm, int smooth = 15, double lam_node = 0.1,
                                 double ratio = 0.33, double lam_ratio = 3.0) {
    Timed tm(T_ENVFIT);
    if (const char* ln = getenv("MS_ME_LAMNODE")) lam_node = atof(ln);   // EXPERIMENT: coupling of neighbouring nodes
    if (opt().ratio >= 0) ratio = opt().ratio;   // minimum neighbour ratio (4.36: 0.2)
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
        const double wloc = opt().lamlocal;
        for (int q = 0; q + 1 < Q; q++) {
            double lq = lam_node;
            if (wloc > 0 && g_env_node_mass && (int)g_env_node_mass->size() == Q) {
                const double m0 = (*g_env_node_mass)[q], m1 = (*g_env_node_mass)[q + 1];
                // MS_ME_LAMLOCAL0=w0: full coupling within w0 (one protein's family and its charge off-by-one
                // ghosts, M (z-1)/z: a few %), falling over w beyond it
                const double w0 = opt().lamlocal0;
                if (m0 > 0 && m1 > 0) { const double t = std::max(0.0, std::fabs(std::log(m1 / m0)) - w0) / wloc; lq = lam_node * std::exp(-t * t); }
                // MS_ME_LAMEND=f: the two ends of the mass range carry no species; their nodes keep at least f of
                // the coupling (a free envelope there fits the noise as a broad floor)
                const double fend = opt().lamend;
                if (q == 0 || q + 2 == Q) lq = std::max(lq, fend * lam_node);
            }
            for (int k = 0; k < K; k++) {
                const size_t i0 = (size_t)q * K + k, i1 = i0 + K;
                G[i0 * P + i0] += lq * scale; G[i1 * P + i1] += lq * scale;
                G[i0 * P + i1] -= lq * scale; G[i1 * P + i0] -= lq * scale;
            }
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
        Timed tm(T_OBJECTIVE);
        size_t na = u.size(), n = y->size();
        f.resize(na); r.resize(n); adj.resize(na); fit.resize(n);
        // (on every core; the sums in fixed blocks: psum)
        double* ff = f.data(); const double* uu = u.data();
        const double dl = dflt, ld = log_dflt;
        double ent;
        {   // f = exp(u) and the entropy in one pass
            Timed te(T_EXP);
            ent = psum(na, [&](size_t j) { const double fj = fexp(uu[j]); ff[j] = fj; return fj - dl - fj * (uu[j] - ld); });
        }
        model->forward(ff, fit.data());
        const double* yy = y->data(); const double* ww = w->data();
        double* rr = r.data();
        const double* ft = fit.data();
        double chi2 = psum(n, [&](size_t i) { double d = ft[i] - yy[i]; rr[i] = ww[i] * d; return ww[i] * d * d; });
        chi2 *= 0.5;
        grad.resize(na);
        double* gg = grad.data(); const double al = alpha;
        // the gradient as the adjoint gives each mass its sum
        model->adjoint_f(r.data(), adj.data(), [&](size_t j, double s) { gg[j] = ff[j] * (s + al * (uu[j] - ld)); });
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

// L-BFGS with bounds: compact form (m pairs), projected search direction,
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
    Timed tm(T_LBFGS);
    // every loop over the n variables runs on every core (sums: psum / pdot, fixed blocks)
    size_t n = x.size();
    auto clip = [&](double v) { return std::min(hi, std::max(lo, v)); };
    pfor(n, [&](size_t a, size_t b) { for (size_t j = a; j < b; j++) x[j] = clip(x[j]); });
    std::vector<double> g, gn, xn(n), d(n), hd(n);
    double fx = obj(x, g);
    std::vector<std::vector<float>> S, Y;    // the pairs, oldest first (single precision: half the memory to read)
    // per pair i, j: SY = s_i'y_j (used for i <= j) and W = y_i' hd y_j (m x m, row major)
    std::vector<double> SY((size_t)m * m), W((size_t)m * m);
    // vectors kept from iteration to iteration (a new one of n doubles each time cost fresh memory pages)
    std::vector<float> s_buf, y_buf;
    int restarts = 0;
    auto set_diag = [&]() {
        pfor(n, [&](size_t a, size_t b) { for (size_t j = a; j < b; j++) hd[j] = 1.0 / std::max(obj.f[j], 1e-300); });
    };
    auto first_direction = [&]() {
        double gmax = pmax(n, 0.0, [&](size_t j) { return std::fabs(hd[j] * g[j]); });
        double sc = gmax > 0 ? 1.0 / gmax : 1.0;
        pfor(n, [&](size_t a, size_t b) { for (size_t j = a; j < b; j++) d[j] = clip(x[j] - sc * hd[j] * g[j]) - x[j]; });
    };
    // Kept from the pass that stores a pair (it reads every vector the next iteration needs first): the
    // projected gradient's largest value and S'g, Y' hd g at the new point
    bool have_pg = false, have_ab = false;
    double pg_next = 0;
    std::vector<double> ab_next;
    const size_t nblk = std::max<size_t>(1, (n + PBLOCK - 1) / PBLOCK);
    std::vector<double> bmax(nblk);
    for (int it = 0; it < maxiter; it++) {
        double pgmax = have_pg ? pg_next : pmax(n, 0.0, [&](size_t i) {
            double gi = g[i];
            if ((x[i] <= lo && gi > 0) || (x[i] >= hi && gi < 0)) gi = 0;
            return std::fabs(gi);
        });
        have_pg = false;
        if (pgmax <= pgtol) { if (getenv("MS_DEBUG")) fprintf(stderr, "  lbfgs pgtol stop it %d f %.6g\n", it, fx); break; }
        if (it == 0) set_diag();   // one metric per call (a fixed variable scaling keeps the pairs consistent)
        const auto tdir0 = std::chrono::steady_clock::now();
        int h = (int)S.size();
        double dg = 0, stpmx_neg = -1e10;
        bool have_dg = false;
        if (h == 0) {
            first_direction();
        } else {
            // H g = gamma hd (g - Y u) + S R^-T ((D + gamma W) u - gamma Y' hd g), u = R^-1 S'g, with R the upper
            // triangle of S'Y and D its diagonal; gamma = s'y / y' hd y of the newest pair
            std::vector<double> ab(2 * h);
            if (have_ab && ab_next.size() == ab.size()) ab = ab_next;
            else psums(n, 2 * h, ab.data(), [&](size_t i0, size_t i1, double* acc) {
                for (int i = 0; i < h; i++) {
                    const float* si = S[i].data();
                    double v = 0;
                    for (size_t j = i0; j < i1; j++) v += si[j] * g[j];
                    acc[i] = v;
                }
                for (int i = 0; i < h; i++) {
                    const float* yi = Y[i].data();
                    double v = 0;
                    for (size_t j = i0; j < i1; j++) v += yi[j] * (hd[j] * g[j]);
                    acc[h + i] = v;
                }
            });
            const size_t M = (size_t)m;
            const double gamma = SY[(h - 1) * M + h - 1] / std::max(W[(h - 1) * M + h - 1], 1e-300);
            std::vector<double> u(h), p(h);
            for (int i = h - 1; i >= 0; i--) {
                double v = ab[i];
                for (int j = i + 1; j < h; j++) v -= SY[i * M + j] * u[j];
                u[i] = v / SY[i * M + i];
            }
            for (int i = 0; i < h; i++) {
                double w = 0;
                for (int j = 0; j < h; j++) w += W[i * M + j] * u[j];
                double v = SY[i * M + i] * u[i] + gamma * w - gamma * ab[h + i];
                for (int j = 0; j < i; j++) v -= SY[j * M + i] * p[j];
                p[i] = v / SY[i * M + i];
            }
            // d = clip(x - H g) - x, in short pieces that stay in the cache while every pair is added; in the
            // same pass d'g (the blocks of psum, each summed in order) and the largest step inside the box
            psums(n, 1, &dg, [&](size_t i0, size_t i1, double* acc) {
                double r[1024];
                double v = 0, mx = -1e10;
                for (size_t c0 = i0; c0 < i1; c0 += 1024) {
                    const size_t L = std::min(i1, c0 + 1024) - c0;
                    const double *gc = g.data() + c0, *hc = hd.data() + c0, *xc = x.data() + c0;
                    double* dc = d.data() + c0;
                    for (size_t j = 0; j < L; j++) r[j] = gc[j];
                    for (int i = 0; i < h; i++) {
                        const float* yi = Y[i].data() + c0; const double ui = u[i];
                        for (size_t j = 0; j < L; j++) r[j] -= ui * yi[j];
                    }
                    for (size_t j = 0; j < L; j++) r[j] *= gamma * hc[j];
                    for (int i = 0; i < h; i++) {
                        const float* si = S[i].data() + c0; const double pi = p[i];
                        for (size_t j = 0; j < L; j++) r[j] += pi * si[j];
                    }
                    for (size_t j = 0; j < L; j++) {
                        const double dj = clip(xc[j] - r[j]) - xc[j];
                        dc[j] = dj;
                        v += dj * gc[j];
                        mx = std::max(mx, dj > 0 ? -((hi - xc[j]) / dj) : dj < 0 ? -((lo - xc[j]) / dj) : -1e10);
                    }
                }
                acc[0] = v;
                bmax[i0 / PBLOCK] = mx;
            });
            for (size_t b = 0; b < nblk; b++) stpmx_neg = std::max(stpmx_neg, bmax[b]);
            have_dg = true;
        }
        have_ab = false;
        prof_t[T_DIR] += std::chrono::duration<double>(std::chrono::steady_clock::now() - tdir0).count(); prof_n[T_DIR]++;
        if (!have_dg) dg = pdot(d.data(), g.data(), n);
        if (!(dg < 0)) {
            S.clear(); Y.clear();
            first_direction();
            dg = pdot(d.data(), g.data(), n);
            have_dg = false;
            if (!(dg < 0)) break;
        }
        // largest step inside the box (the smallest ratio: any order gives the same minimum)
        double stpmx = -(have_dg ? stpmx_neg : pmax(n, -1e10, [&](size_t j) {
            if (d[j] > 0) return -((hi - x[j]) / d[j]);
            if (d[j] < 0) return -((lo - x[j]) / d[j]);
            return -1e10;
        }));
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
                pfor(n, [&](size_t a, size_t b) { for (size_t j = a; j < b; j++) xn[j] = clip(x[j] + stp * d[j]); });
                fn = obj(xn, gn);
                nev++;
                double stp_eval = stp;
                double gd = pdot(gn.data(), d.data(), n);
                if (!std::isfinite(fn)) { fn = 1e300; gd = 1e300; }
                task = ls.iterate(stp, fn, gd);
                if (task == Dcsrch::CONV || task == Dcsrch::WARN) {
                    if (stp != stp_eval) {   // the best point of the search: evaluate it
                        pfor(n, [&](size_t a, size_t b) { for (size_t j = a; j < b; j++) xn[j] = clip(x[j] + stp * d[j]); });
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
            S.clear(); Y.clear();
            restarts++;
            continue;
        }
        Timed tpair(T_PAIR);
        s_buf.resize(n); y_buf.resize(n);
        float* s = s_buf.data(); float* yk = y_buf.data();
        // in one pass: s and y; s_i'y and y_i' hd y of the pairs kept, then s'y, y'y and y' hd y; for the next
        // iteration S'g and Y' hd g at the new point (old pairs, then the new one) and the projected gradient
        const int hp = (int)S.size();
        const int na2 = 4 * hp + 5;
        std::vector<double> nw(na2);
        psums(n, na2, nw.data(), [&](size_t i0, size_t i1, double* acc) {
            double pg = 0;
            for (size_t j = i0; j < i1; j++) {
                s[j] = (float)(xn[j] - x[j]); yk[j] = (float)(gn[j] - g[j]);
                double gi = gn[j];
                if ((xn[j] <= lo && gi > 0) || (xn[j] >= hi && gi < 0)) gi = 0;
                pg = std::max(pg, std::fabs(gi));
            }
            bmax[i0 / PBLOCK] = pg;
            for (int i = 0; i < hp; i++) {
                const float* si = S[i].data();
                double v = 0, w = 0;
                for (size_t j = i0; j < i1; j++) { v += (double)si[j] * yk[j]; w += si[j] * gn[j]; }
                acc[i] = v; acc[2 * hp + 3 + i] = w;
            }
            for (int i = 0; i < hp; i++) {
                const float* yi = Y[i].data();
                double v = 0, w = 0;
                for (size_t j = i0; j < i1; j++) { v += yi[j] * (hd[j] * yk[j]); w += yi[j] * (hd[j] * gn[j]); }
                acc[hp + i] = v; acc[3 * hp + 4 + i] = w;
            }
            double sy = 0, yy = 0, yhy = 0, sg = 0, yg = 0;
            for (size_t j = i0; j < i1; j++) {
                sy += (double)s[j] * yk[j]; yy += (double)yk[j] * yk[j]; yhy += yk[j] * (hd[j] * yk[j]);
                sg += s[j] * gn[j]; yg += yk[j] * (hd[j] * gn[j]);
            }
            acc[2 * hp] = sy; acc[2 * hp + 1] = yy; acc[2 * hp + 2] = yhy;
            acc[3 * hp + 3] = sg; acc[4 * hp + 4] = yg;
        });
        pg_next = 0;
        for (size_t b = 0; b < nblk; b++) pg_next = std::max(pg_next, bmax[b]);
        have_pg = true;
        const double sy = nw[2 * hp];
        const bool keep = sy > 2.220446049250313e-16 * nw[2 * hp + 1];
        {   // S'g and Y' hd g of the pairs of the next iteration
            const int first = keep && hp >= m ? 1 : 0, last = keep ? hp + 1 : hp, hn = last - first;
            ab_next.assign(2 * hn, 0.0);
            for (int i = first; i < last; i++) {
                ab_next[i - first] = nw[2 * hp + 3 + i];
                ab_next[hn + i - first] = nw[3 * hp + 4 + i];
            }
            have_ab = true;
        }
        if (keep) {
            const size_t M = (size_t)m;
            std::vector<float> old_s, old_y;    // the oldest pair: its vectors take the next pair
            int first = 0;                      // the first pair kept
            if (hp >= m) {
                old_s = std::move(S.front()); old_y = std::move(Y.front());
                S.erase(S.begin()); Y.erase(Y.begin());
                for (int i = 1; i < hp; i++)
                    for (int j = 1; j < hp; j++) { SY[(i - 1) * M + j - 1] = SY[i * M + j]; W[(i - 1) * M + j - 1] = W[i * M + j]; }
                first = 1;
            }
            const size_t k = S.size();   // the new pair
            for (int i = first; i < hp; i++) {
                const size_t r = (size_t)(i - first);
                SY[r * M + k] = nw[i];
                W[r * M + k] = W[k * M + r] = nw[hp + i];
            }
            SY[k * M + k] = sy; W[k * M + k] = nw[2 * hp + 2];
            S.push_back(std::move(s_buf)); Y.push_back(std::move(y_buf));
            s_buf = std::move(old_s); y_buf = std::move(old_y);
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
    // EXPERIMENT (MS_ME_NOISE2): windows may hold zeros (up to half: counting data of a single quadrupole has them
    // between the ions), and the median based scatter is scaled by the tail factor of the quiet half of the spectrum
    // (mean square against median based scatter: 1 for Gaussian noise, 5 to 100 for spiky data, whose variance the
    // median does not see)
    const bool noise2 = opt().noise2;
    // windows with zeros only when too few windows are free of them (counting data); where zeros come from a
    // threshold (a clipped baseline) the windows without zeros measure the noise and the others would bias it
    std::vector<double> t_mean, t_ms, t_md;
    bool used_zeros = false;
    for (int zpass = 0; zpass < (noise2 ? 2 : 1); zpass++) {
    used_zeros = zpass == 1;
    if (zpass == 1) { if (mw.size() >= 40) break; mw.clear(); vw.clear(); t_mean.clear(); t_ms.clear(); t_md.clear(); }
    const bool allow_zeros = zpass == 1;
    for (size_t a = 0; a + W <= n; a += W) {
        bool ok = true; double mean = 0; size_t zeros = 0;
        for (size_t i = a; i < a + W; i++) {
            if (!keep.empty() && !keep[i]) { ok = false; break; }
            if (!(s.y[i] > 0)) { if (!allow_zeros) { ok = false; break; } zeros++; }
            mean += s.y[i];
        }
        if (!ok || (allow_zeros && (2 * zeros > W || !(mean > 0)))) continue;
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
        if (noise2 && sd > 0) { t_mean.push_back(mean); t_ms.push_back(rms * rms); t_md.push_back(sd * sd); }
        if (!(sd > 0) || rms > 2.5 * sd) continue;   // a peak inside the window
        mw.push_back(mean); vw.push_back(sd * sd);
    }
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
    if (noise2 && t_mean.size() >= 40) {
        std::vector<double> sm = t_mean;
        const double mmed = ms::e2::median(sm);
        std::vector<double> a1, b1;
        for (size_t i = 0; i < t_mean.size(); i++) if (t_mean[i] <= mmed) { a1.push_back(t_ms[i]); b1.push_back(t_md[i]); }
        const double T = std::min(1000.0, std::max(1.0, ms::e2::median(a1) / std::max(ms::e2::median(b1), 1e-300)));
        ab[0] *= T; ab[1] *= T;
        if (getenv("MS_DEBUG")) fprintf(stderr, "  noise from data: tail factor %.3g, floor %.3g, per intensity %.3g (raw units)\n", T, std::sqrt(ab[0]), ab[1]);
    }
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

// EXPERIMENT (MS_MEMSYS): eigenvalues and vectors of a small symmetric matrix (Jacobi rotations; n <= 8):
// A (n x n, row major) -> ev (ascending order not kept), V (columns: the vectors, row major n x n)
static void small_eigen(std::vector<double> A, int n, std::vector<double>& ev, std::vector<double>& V) {
    V.assign((size_t)n * n, 0.0);
    for (int i = 0; i < n; i++) V[(size_t)i * n + i] = 1.0;
    for (int sweep = 0; sweep < 100; sweep++) {
        double off = 0;
        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++) off += A[(size_t)i * n + j] * A[(size_t)i * n + j];
        if (off < 1e-30) break;
        for (int p = 0; p < n; p++)
            for (int q = p + 1; q < n; q++) {
                double apq = A[(size_t)p * n + q];
                if (std::fabs(apq) < 1e-300) continue;
                double app = A[(size_t)p * n + p], aqq = A[(size_t)q * n + q];
                double th = 0.5 * (aqq - app) / apq;
                double t = (th >= 0 ? 1.0 : -1.0) / (std::fabs(th) + std::sqrt(th * th + 1.0));
                double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
                for (int k = 0; k < n; k++) {   // A <- J' A J
                    double akp = A[(size_t)k * n + p], akq = A[(size_t)k * n + q];
                    A[(size_t)k * n + p] = c * akp - s * akq;
                    A[(size_t)k * n + q] = s * akp + c * akq;
                }
                for (int k = 0; k < n; k++) {
                    double apk = A[(size_t)p * n + k], aqk = A[(size_t)q * n + k];
                    A[(size_t)p * n + k] = c * apk - s * aqk;
                    A[(size_t)q * n + k] = s * apk + c * aqk;
                }
                for (int k = 0; k < n; k++) {
                    double vkp = V[(size_t)k * n + p], vkq = V[(size_t)k * n + q];
                    V[(size_t)k * n + p] = c * vkp - s * vkq;
                    V[(size_t)k * n + q] = s * vkp + c * vkq;
                }
            }
    }
    ev.resize(n);
    for (int i = 0; i < n; i++) ev[i] = A[(size_t)i * n + i];
}

// ============================================================ driver
void run_maxent(const double* mz, const double* it, long n0, const ms_maxent_params& p, ms_maxent_result& out,
                ms_progress_fn cb, void* user) {
    using namespace ms::e2;
    const bool debug = getenv("MS_DEBUG") != nullptr;
    struct PrintTimes {  // (made first: printed after the total below is added)
        PrintTimes() { for (double& t : prof_t) t = 0; for (long& c : prof_n) c = 0; }
        ~PrintTimes() {
            if (!getenv("MS_PROFILE")) return;
            for (int k = 0; k < T_N; k++) fprintf(stderr, "  %-18s %8.2f s  %6ld calls\n", prof_name[k], prof_t[k], prof_n[k]);
        }
    } print_times;
    Timed total_time(T_TOTAL);
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
    // noise_model 0: the 4.36 engine; 1 or 4: as in 4.35 (floor + counting, re-estimated); 2, 3: comparison models
    const bool legacy = p.noise_model != 0;
    const int noise_kind = p.noise_model >= 0 && p.noise_model <= 3 ? p.noise_model : 0;
    const MaxentOptions run_opt = maxent_options(legacy);
    struct OptScope {
        const MaxentOptions* prev;
        explicit OptScope(const MaxentOptions* o) : prev(g_opt) { g_opt = o; }
        ~OptScope() { g_opt = prev; }
    } opt_scope(&run_opt);
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
    if (getenv("MS_PROFILE"))
        fprintf(stderr, "  model: %zu grid points, %zu masses, %d charges, %zu entries, %d threads\n", model.n, model.na,
                model.K, model.nnz, omp_threads());
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
        // EXPERIMENT (MS_ME_NOISE="a,b"): the noise measured outside, from the scan to scan variance of the averaged
        // scans (independent of any fit): variance of one averaged raw point = a + b y (raw intensity units),
        // converted to the grid as noise_from_data converts its own estimate
        if (const char* ev = getenv("MS_ME_NOISE")) {
            double va = 0, vb = 0;
            if (sscanf(ev, "%lf,%lf", &va, &vb) == 2 && va >= 0 && vb >= 0) {
                const double F = (d_mz > 0 && dx_raw > 0) ? std::min(2.0 / 3.0, dx_raw / d_mz) : 2.0 / 3.0;
                noise.s0 = std::sqrt(F * va) / scale; noise.ginv = F * vb / scale; noise.eps = 0; noise.F = F;
                noise.from_data = true;
                if (debug) fprintf(stderr, "  noise from outside: s0 %.4g, 1/g %.4g (grid units, F %.3g)\n", noise.s0, noise.ginv, F);
            }
        }
        if (noise_kind == 2) noise.s0 = std::max(noise_old, 2e-3);
        else noise.s0 = std::max(noise.s0, 1e-4);
    }
    const NoiseModel noise_data = noise;   // lower bounds for the re-estimation
    // EXPERIMENT (MS_ME_KAPPA): instead of the re-estimation from the residuals binned by intensity, the noise measured
    // on the data (floor, counting term) plus a model-error term (kappa y)^2, kappa = 1.4826 x the median of
    // |residual| / intensity over the points above a tenth of the largest intensity (the peak shape, adducts and
    // isotopes the model does not have: its error grows with the signal, it is not noise)
    const bool kappa_mode = getenv("MS_ME_KAPPA") != nullptr;
    auto kappa_noise = [&](const std::vector<double>& yb, const std::vector<double>& r) {
        NoiseModel nm = noise_data;
        double ymax = 0;
        for (size_t i = 0; i < yb.size(); i++) ymax = std::max(ymax, yb[i]);
        std::vector<double> rel;
        for (size_t i = 0; i < yb.size(); i++) if (yb[i] >= 0.1 * ymax && yb[i] > 0) rel.push_back(std::fabs(r[i]) / yb[i]);
        if (rel.size() >= 10) nm.eps = 1.4826 * median(rel);
        if (getenv("MS_ME_KAPPA_SCALE")) nm.eps *= atof(getenv("MS_ME_KAPPA_SCALE"));
        return nm;
    };
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
    // EXPERIMENT (MS_ME_SNRMASK): a point has signal when its intensity is 3 sd of its own noise above zero (floor
    // and counting term): with counting noise most of a spectrum is weak baseline that an empty result fits within
    // its noise, and the average over those points would stop the search before any peak is fitted
    const bool snrmask = opt().snrmask;
    for (size_t i = 0; i < ng; i++) signal[i] = snrmask ? ybar[i] > 3.0 * std::sqrt(noise.var(ybar[i])) : ybar[i] > 3.0 * noise.s0;
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
        pfor(na, [&](size_t a0, size_t b0) { for (size_t j = a0; j < b0; j++) f[j] = fexp(u[j]); });
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
    std::vector<double> seed_used;   // EXPERIMENT (MS_ME_SEEDNODES): nodes placed at seeds by the last species_nodes()
    auto refit_envelope = [&]() {
        model.node_normal(f.data(), y.data(), w.data(), G, bvec);
        g_env_node_mass = &model.node_mass;
        const std::vector<double> prev = model.envelope();
        std::vector<double> cn = envelope_fit(G, bvec, model.Q, K, charges.front(), prev);
        g_env_node_mass = nullptr;
        // a node whose masses carry (almost) no intensity yet has nothing to fit its envelope to; it keeps the
        // envelope it had (a seed node's start from the data), rescaled with the others
        const bool keep_unsupported = opt().seednodes;
        if (keep_unsupported && model.Q > 1) {
            const size_t P = (size_t)model.Q * K;
            std::vector<double> sup(model.Q, 0.0);
            double smax = 0;
            for (int q = 0; q < model.Q; q++) {
                for (int k = 0; k < K; k++) sup[q] += G[((size_t)q * K + k) * P + (size_t)q * K + k];
                smax = std::max(smax, sup[q]);
            }
            double an = 0, ap = 0;
            for (int q = 0; q < model.Q; q++)
                if (sup[q] >= 1e-4 * smax) for (int k = 0; k < K; k++) { an += cn[(size_t)q * K + k]; ap += prev[(size_t)q * K + k]; }
            const double rs = ap > 0 ? an / ap : 1.0;
            for (int q = 0; q < model.Q; q++)
                if (sup[q] < 1e-4 * smax) for (int k = 0; k < K; k++) cn[(size_t)q * K + k] = rs * prev[(size_t)q * K + k];
        }
        // a charge state of a lighter node that falls on a charge state of a stronger, heavier species (both within a
        // species' width in mass) is that species' signal: the node's envelope is zero there. With envelopes free between distant
        // masses, a node at M/2 put its whole envelope on the charges z where M sits at 2z: a 12229 Da protein at
        // m/z 599-1551 (mass range 5000-20000 Da) gave a 6115 Da half-mass ghost at charges 8-10, whose other
        // charges fall below the window (4-8 % of the signal, 4.36). A real species there keeps its charge states
        // where the stronger one has none (a half antibody at z, the antibody's envelope far from 2z).
        if (opt().seedharm && model.Q > 1) {
            const size_t P = (size_t)model.Q * K;
            std::vector<double> sup(model.Q, 0.0);
            double smax = 0;
            for (int q = 0; q < model.Q; q++) {
                for (int k = 0; k < K; k++) sup[q] += G[((size_t)q * K + k) * P + (size_t)q * K + k] * cn[(size_t)q * K + k];
                smax = std::max(smax, sup[q]);
            }
            const double hw = 0.0005 * 0.5 * (masses.front() + masses.back());
            const double mz_win_lo = use_log ? std::exp(coord.front()) : coord.front();
            const double mz_win_hi = use_log ? std::exp(coord.back()) : coord.back();
            for (int j = 0; j < model.Q; j++) {
                if (!(sup[j] >= 0.1 * smax) || smax <= 0) continue;
                const double Mj = model.node_mass[j];
                for (int q = 0; q < model.Q; q++) {
                    const double Mq = model.node_mass[q];
                    if (q == j || Mq >= Mj || sup[q] >= sup[j] || std::fabs(std::log(Mj / Mq)) < 0.02) continue;   // lighter only
                    int m = 0;   // Mq = Mj / m (m = 2, 3, 4): a harmonic; only then do its charge states z sit on Mj's m z
                    for (int mm = 2; mm <= 4 && !m; mm++) if (std::fabs(Mq - Mj / mm) < hw) m = mm;
                    if (!m) continue;
                    // only when the data cannot tell: every charge z of the node inside the m/z window sits on a charge
                    // m z of Mj within the charge range (Mj alone can explain all of it); a real species at Mj / m (the
                    // haemoglobin dimer next to the tetramer) has charges in the window whose m z lies beyond the
                    // charge range, which the fit has to answer to
                    bool falsifiable = false;
                    std::vector<int> on;
                    for (int k = 0; k < K && !falsifiable; k++) {
                        const double pos = (Mq + charges[k] * a) / charges[k];
                        if (pos < mz_win_lo || pos > mz_win_hi) continue;
                        const int kk = charges[k] * m - charges.front();   // index of charge m z
                        if (kk >= 0 && kk < K) on.push_back(k);   // (Mj may carry it there: within the charge range)
                        else falsifiable = true;
                    }
                    if (falsifiable || on.empty()) continue;
                    for (int k : on) {
                        if (getenv("MS_DEBUG") && cn[(size_t)q * K + k] > 0)
                            fprintf(stderr, "  harmonic: node %.1f Da = %.1f / %d, charge %d off (%.3g)\n", Mq, Mj, m, charges[k], cn[(size_t)q * K + k]);
                        cn[(size_t)q * K + k] = 0;
                    }
                }
            }
        }
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
        seed_used.clear();
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
        int qmax = std::max(2, std::min(12, qcost));
        // EXPERIMENT (MS_ME_SEEDNODES): a species the first pass did not fit has no node, so the envelope there is
        // interpolated from unrelated proteins and the final pass cannot fit it either (measured: GroEL monomer next
        // to the 13-mer, IgG4 half-molecule next to the antibody, carbonic anhydrase in a protein mix). Seeds: the
        // signal the fit leaves, projected back onto mass with every charge weighted alike, local maxima above 2 % of
        // the same projection of the data; they come after the species found and get at least 4 nodes in all
        const bool seeds = opt().seednodes;
        std::vector<std::pair<double, double>> seedc;
        if (seeds) {
            // per charge the largest value within a peak width of the mass's m/z (v: the data or what the fit leaves)
            const long rad = std::max(1L, (long)(2 * model.sig));
            auto sample = [&](double M, const std::vector<double>& v, std::vector<double>& sk) {
                for (int k = 0; k < K; k++) {
                    const double z = charges[k];
                    double pos = (M + z * a) / z;
                    if (use_log) pos = std::log(std::max(pos, 1e-9));
                    long ic = (long)(std::upper_bound(coord.begin(), coord.end(), pos) - coord.begin()) - 1;
                    double best = 0;
                    if (ic >= 0 && ic < (long)ng)
                        for (long i = std::max(0L, ic - rad); i <= std::min((long)ng - 1, ic + rad); i++) best = std::max(best, v[i]);
                    sk[k] = best;
                }
            };
            // a species shows at neighbouring charges: the score sums min(s_z, s_z+1), which a mass explained by
            // scattered residuals of other species (every other charge, random coincidences) does not reach
            std::vector<double> rp(ng);
            for (size_t i = 0; i < ng; i++) rp[i] = std::max(0.0, y[i] - fitv[i]);
            const double dl = 0.0005;   // coarse mass grid, 0.05 % steps
            const double l0 = std::log(masses.front()), l1 = std::log(masses.back());
            const int nc = std::max(3, (int)((l1 - l0) / dl));
            std::vector<double> sr(nc), sy(nc), skr(K), sky(K);
            for (int c = 0; c < nc; c++) {
                const double M = std::exp(l0 + (l1 - l0) * c / (nc - 1));
                sample(M, rp, skr); sample(M, y, sky);
                double a1 = 0, a2 = 0;
                for (int k = 0; k + 1 < K; k++) { a1 += std::min(skr[k], skr[k + 1]); a2 += std::min(sky[k], sky[k + 1]); }
                sr[c] = a1; sy[c] = a2;
            }
            const double symax = *std::max_element(sy.begin(), sy.end());
            // significance: the score of noise alone is the typical (median) score of all masses; a seed must stand
            // 3 times above it as well as at 2 % of the data's best (a noisy single scan otherwise seeds on noise)
            std::vector<double> srs = sr;
            const double sbg = ms::e2::median(srs);
            const double seed_min = opt_num("MS_ME_SEEDMIN", 0.02), seed_bg = opt_num("MS_ME_SEEDBG", 3.0);   // (development)
            for (int c = 1; c + 1 < nc; c++)
                if (sr[c] > sr[c - 1] && sr[c] >= sr[c + 1] && sr[c] >= seed_min * symax && sr[c] >= seed_bg * sbg)
                    seedc.push_back({sr[c], std::exp(l0 + (l1 - l0) * c / (nc - 1))});
            std::sort(seedc.begin(), seedc.end(), [](const std::pair<double, double>& p1, const std::pair<double, double>& p2) { return p1.first > p2.first; });
            if (getenv("MS_DEBUG"))
                for (size_t s2 = 0; s2 < std::min<size_t>(seedc.size(), 6); s2++)
                    fprintf(stderr, "  seed %.1f Da (%.3g of the data score)\n", seedc[s2].second, seedc[s2].first / std::max(symax, 1e-300));
        }
        auto add_nodes = [&](const std::vector<std::pair<double, double>>& cs, bool is_seed) {
            for (auto& c : cs) {
                if ((int)nodes.size() >= qmax) break;
                bool near = false;
                for (double m : nodes) if (std::fabs(std::log(c.second / m)) < 0.02) { near = true; break; }
                if (!near) { nodes.push_back(c.second); if (is_seed) seed_used.push_back(c.second); }
            }
        };
        add_nodes(cand, false);
        if (!seedc.empty()) qmax = std::max(qmax, (int)nodes.size() + std::min<int>(2, (int)seedc.size()));   // room for up to two seeds
        add_nodes(seedc, true);
        std::sort(nodes.begin(), nodes.end());
        return nodes;
    };
    // the envelope of a node placed at a seed starts from the data, not from its neighbours: per charge the largest
    // signal the fit leaves within a peak width of the seed's m/z, smoothed over neighbouring charges (a copy of
    // the neighbours' charges would put the seed's intensity where its peaks are not, and the refit, which only
    // sees the intensity the masses already carry, could not move it)
    auto seed_envelopes = [&]() {
        if (seed_used.empty()) return;
        std::vector<double> E = model.envelope();
        const long rad = std::max(1L, (long)(2 * model.sig));
        for (int q = 0; q < model.Q; q++) {
            const double M = model.node_mass[q];
            if (std::none_of(seed_used.begin(), seed_used.end(), [&](double sm) { return std::fabs(sm - M) < 1e-6 * M; })) continue;
            std::vector<double> sc(K, 0.0);
            for (int k = 0; k < K; k++) {
                const double z = charges[k];
                double pos = (M + z * a) / z;
                if (use_log) pos = std::log(std::max(pos, 1e-9));
                long i0 = (long)(std::upper_bound(coord.begin(), coord.end(), pos) - coord.begin()) - 1;
                if (i0 < 0 || i0 >= (long)ng) continue;
                for (long i = std::max(0L, i0 - rad); i <= std::min((long)ng - 1, i0 + rad); i++) sc[k] = std::max(sc[k], y[i] - fitv[i]);
            }
            double tot = 0;
            std::vector<double> sm(K, 0.0);
            for (int k = 0; k < K; k++) {
                sm[k] = 0.5 * std::max(0.0, sc[k]) + 0.25 * std::max(0.0, k > 0 ? sc[k - 1] : 0.0) + 0.25 * std::max(0.0, k + 1 < K ? sc[k + 1] : 0.0);
                tot += sm[k];
            }
            if (!(tot > 0)) continue;
            double rowsum = 0;
            for (int k = 0; k < K; k++) rowsum += E[(size_t)q * K + k];
            for (int k = 0; k < K; k++) E[(size_t)q * K + k] = (rowsum > 0 ? rowsum : 1.0) * sm[k] / tot;
            if (getenv("MS_DEBUG")) {
                int kb = (int)(std::max_element(sm.begin(), sm.end()) - sm.begin());
                fprintf(stderr, "  seed node %.1f Da: envelope from the data, peak at charge %d\n", M, (int)charges[kb]);
            }
        }
        model.set_envelope(E);
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
    // EXPERIMENT (MS_MEMSYS): the alpha search of Skilling and Bryan (1984, the "Cambridge algorithm" of
    // MemSys): one path instead of rounds of L-BFGS runs. Each step builds 4 search directions in the entropy
    // metric (f grad S, f grad C, and the curvature of C applied to each: f A'WA of them), models S and C as
    // quadratics in their span and takes the step (and alpha) that brings C two thirds of the way down to the
    // least C within a trust region (distance^2 <= r * sum f in the metric 1/f). The pass ends when C stops
    // falling (the best fit) or reaches the target; as in the rounds, the noise is then estimated again and the
    // envelope nodes put at the species. C here: sum w (A f - y)^2 over the grid.
    const bool memsys = opt().memsys;
    // safe (MS_MEMSYS_SAFE): the first pass as in 4.33 (the rounds, the same noise re-estimation and nodes), the
    // path for the final pass only; the noise is then decided exactly as 4.33 decides it
    const bool hybrid = opt().hybrid;
    int memsys_start = 0;
    if (memsys && !hybrid && !getenv("MS_MEMSYS_NOWIN")) model.make_windows();
    auto memsys_run = [&]() {
        Timed tmm(T_LBFGS);   // (profile: the whole path)
        const int ND = getenv("MS_MEMSYS_ND") ? std::max(2, std::min(8, atoi(getenv("MS_MEMSYS_ND")))) : 5;   // directions (5: the masses of 4.33 within 0.5 mDa more often than with 4)
        const double rfac = getenv("MS_MEMSYS_R") ? atof(getenv("MS_MEMSYS_R")) : 10.0;
        const int env_every = getenv("MS_MEMSYS_ENV") ? std::max(1, atoi(getenv("MS_MEMSYS_ENV"))) : 5;
        const int max_it = getenv("MS_MEMSYS_IT") ? atoi(getenv("MS_MEMSYS_IT")) : 300;   // (a cap: passes end when settled)
        const double flat_tol = getenv("MS_MEMSYS_FLAT") ? atof(getenv("MS_MEMSYS_FLAT")) : 2e-3;
        const double flat_tol1 = getenv("MS_MEMSYS_FLAT1") ? atof(getenv("MS_MEMSYS_FLAT1")) : 1e-2;   // first pass
        const double test_tol = getenv("MS_MEMSYS_TEST") ? atof(getenv("MS_MEMSYS_TEST")) : 1.0;
        const bool adapt = getenv("MS_MEMSYS_ADAPT") != nullptr;   // the trust region grows when steps go as predicted
        // inner / outer: steps at a fixed envelope and fixed weights until C settles, then one envelope refit;
        // the pass ends when a refit no longer moves the settled C (instead of a refit every env_every steps)
        const bool inner = getenv("MS_MEMSYS_NOINNER") == nullptr;
        const int inner_max = getenv("MS_MEMSYS_INNER_MAX") ? atoi(getenv("MS_MEMSYS_INNER_MAX")) : 40;
        const double outer_tol = getenv("MS_MEMSYS_OUTER") ? atof(getenv("MS_MEMSYS_OUTER")) : 5e-3;
        // steps that would lower a mass: multiplicative (f exp(d / f)), so that no mass is clipped at 0.1 f
        const bool logstep = getenv("MS_MEMSYS_LOGSTEP") != nullptr;
        // how far each step goes: C aims at Cmin + frac (C0 - Cmin) (frac 1/3: two thirds of the way down), in the
        // first pass and in the final pass; push: steps at frac 0 (the best fit within the trust region) at the end
        // of the first pass, so that the noise is estimated from residuals of a fit pushed as hard as 4.33's rounds
        const double aim1 = getenv("MS_MEMSYS_AIM1") ? atof(getenv("MS_MEMSYS_AIM1")) : 1.0 / 3.0;
        const double aimF = getenv("MS_MEMSYS_AIM") ? atof(getenv("MS_MEMSYS_AIM")) : 1.0 / 3.0;
        const int push_n = getenv("MS_MEMSYS_PUSH") ? atoi(getenv("MS_MEMSYS_PUSH")) : 0;
        // a step that raised C is taken back and the trust region made 4 times smaller (it grows back by 2 after
        // each step that lowered C)
        const bool reject = getenv("MS_MEMSYS_NOREJECT") == nullptr;
        std::vector<double> f_keep, fit_keep;
        bool have_keep = false;
        double C_keep = 0;
        int n_rej = 0;
        double rcur = rfac, pred_prev = -1;
        bool limited_prev = false;
        const double flo = std::exp(lo_u), fhi = std::exp(hi_u);
        std::vector<std::vector<double>> e(ND, std::vector<double>(na)), Ae(ND, std::vector<double>(ng));
        std::vector<double> gS(na), gC(na), wr(ng), back(ng);
        int total_it = 0, n_clip = 0;
        double C_last = 0;
        auto sync_u = [&]() { pfor(na, [&](size_t a0, size_t b0) { for (size_t j = a0; j < b0; j++) u[j] = std::log(f[j]); }); };
        auto sync_f = [&]() { pfor(na, [&](size_t a0, size_t b0) { for (size_t j = a0; j < b0; j++) f[j] = fexp(u[j]); }); };
        sync_f();
        for (pass = memsys_start; pass < max_passes; pass++) {
            const bool final_pass = pass + 1 >= max_passes;
            double C_prev = -1;
            int flat_n = 0;
            bool fit_ok = false;   // fitv is A f already
            std::string why = "steps";
            bool do_refit = false;
            int since_refit = 0;
            double C_settled = -1;
            int push_left = -1;
            for (int it = 0; it < max_it; it++, total_it++, since_refit++) {
                if (inner ? do_refit : (it > 0 && it % env_every == 0)) {   // the envelope and the weights follow the fit
                    do_refit = false;
                    since_refit = 0;
                    flat_n = 0;
                    have_keep = false;
                    sync_u();
                    refit_envelope();
                    sync_f();
                    set_weights();
                    C_prev = -1;   // C changed with them
                    fit_ok = false;
                }
                // C and its gradient at f
                if (!fit_ok) model.forward(f.data(), fitv.data());
                have_fit = true;
                const double C0 = psum(ng, [&](size_t i) { double d = fitv[i] - y[i]; wr[i] = 2.0 * w[i] * d; return w[i] * d * d; });
                if (reject && have_keep) {
                    have_keep = false;
                    if (C0 > C_keep) {   // worse: back to the point before the step, a smaller trust region
                        f = f_keep;
                        sync_u();
                        fitv = fit_keep;
                        fit_ok = true;
                        rcur = std::max(rcur * 0.25, 1e-4);
                        n_rej++;
                        if (debug) fprintf(stderr, "memsys step rejected: C %.6g > %.6g, trust region now %.3g\n", C0, C_keep, rcur);
                        continue;
                    }
                    rcur = std::min(rcur * 2.0, rfac);
                }
                const double chi2n = chi2_signal();
                model.adjoint_f(wr.data(), gC.data());
                // grad S = -ln(f/m) = -(u - ln m); the first two directions, f grad S and f grad C
                double raw[3];   // sum f gS^2, sum f gC^2, sum f gS gC (the metric lengths and their angle)
                psums(na, 3, raw, [&](size_t i0, size_t i1, double* acc) {
                    double a = 0, b = 0, c = 0;
                    for (size_t j = i0; j < i1; j++) {
                        const double s = -(u[j] - log_dflt), g = gC[j];
                        gS[j] = s;
                        e[0][j] = f[j] * s;
                        e[1][j] = f[j] * g;
                        a += f[j] * s * s; b += f[j] * g * g; c += f[j] * s * g;
                    }
                    acc[0] = a; acc[1] = b; acc[2] = c;
                });
                const double test = 1.0 - raw[2] / std::sqrt(std::max(raw[0] * raw[1], 1e-300));   // 0 at the solution for the alpha
                const double alpha_eff = 0.5 * std::sqrt(raw[1] / std::max(raw[0], 1e-300));
                // (the directions keep their own lengths: the small matrices are scaled instead, below)
                for (int k = 0; k < 2; k++) model.forward(e[k].data(), Ae[k].data());
                // the curvature of C on them: f A'W A e
                for (int k = 2; k < ND; k++) {   // direction k: the curvature on direction k - 2
                    pfor(ng, [&](size_t a0, size_t b0) { for (size_t i = a0; i < b0; i++) back[i] = 2.0 * w[i] * Ae[k - 2][i]; });
                    double* ek = e[k].data();
                    model.adjoint_f(back.data(), ek, [&](size_t j, double v) { ek[j] = f[j] * v; });
                    model.forward(ek, Ae[k].data());
                }
                // the quadratic models in the span: G (metric), s (grad S), c (grad C), M (curvature of C)
                const int NG = ND * (ND + 1) / 2;
                std::vector<double> sums(NG + 2 * ND + 1), msum(NG);
                psums(na, NG + 2 * ND + 1, sums.data(), [&](size_t i0, size_t i1, double* acc) {
                    int q = 0;
                    for (int a = 0; a < ND; a++)
                        for (int b = a; b < ND; b++) {
                            double v = 0;
                            for (size_t j = i0; j < i1; j++) v += e[a][j] * e[b][j] / f[j];
                            acc[q++] = v;
                        }
                    for (int a = 0; a < ND; a++) {
                        double vs = 0, vc = 0;
                        for (size_t j = i0; j < i1; j++) { vs += e[a][j] * gS[j]; vc += e[a][j] * gC[j]; }
                        acc[NG + a] = vs; acc[NG + ND + a] = vc;
                    }
                    double sf = 0;
                    for (size_t j = i0; j < i1; j++) sf += f[j];
                    acc[NG + 2 * ND] = sf;
                });
                psums(ng, NG, msum.data(), [&](size_t i0, size_t i1, double* acc) {
                    int q = 0;
                    for (int a = 0; a < ND; a++)
                        for (int b = a; b < ND; b++) {
                            double v = 0;
                            for (size_t i = i0; i < i1; i++) v += w[i] * Ae[a][i] * Ae[b][i];
                            acc[q++] = 2.0 * v;
                        }
                });
                std::vector<double> G((size_t)ND * ND), M((size_t)ND * ND), s(ND), c(ND);
                {
                    int q = 0;
                    for (int a = 0; a < ND; a++)
                        for (int b = a; b < ND; b++, q++) {
                            G[(size_t)a * ND + b] = G[(size_t)b * ND + a] = sums[q];
                            M[(size_t)a * ND + b] = M[(size_t)b * ND + a] = msum[q];
                        }
                    for (int a = 0; a < ND; a++) { s[a] = sums[NG + a]; c[a] = sums[NG + ND + a]; }
                }
                const double sum_f = sums[NG + 2 * ND];
                // the directions at unit length in the metric: the matrices scaled (the steps scaled back below)
                std::vector<double> dsc(ND);
                for (int a = 0; a < ND; a++) dsc[a] = 1.0 / std::sqrt(std::max(G[(size_t)a * ND + a], 1e-300));
                for (int a = 0; a < ND; a++) {
                    s[a] *= dsc[a];
                    c[a] *= dsc[a];
                    for (int b = 0; b < ND; b++) {
                        G[(size_t)a * ND + b] *= dsc[a] * dsc[b];
                        M[(size_t)a * ND + b] *= dsc[a] * dsc[b];
                    }
                }
                // a basis orthonormal in G (directions that add nothing are left out), then M diagonal in it
                std::vector<double> gev, gV;
                small_eigen(G, ND, gev, gV);
                double gmax = 0;
                for (double v : gev) gmax = std::max(gmax, v);
                std::vector<int> keep;
                for (int k = 0; k < ND; k++) if (gev[k] > 1e-8 * gmax) keep.push_back(k);
                const int dim = (int)keep.size();
                std::vector<double> B((size_t)ND * dim);   // columns: G^-1/2 directions
                for (int col = 0; col < dim; col++)
                    for (int a = 0; a < ND; a++) B[(size_t)a * dim + col] = gV[(size_t)a * ND + keep[col]] / std::sqrt(gev[keep[col]]);
                std::vector<double> Mp((size_t)dim * dim, 0.0), sp(dim, 0.0), cp(dim, 0.0);
                for (int i = 0; i < dim; i++) {
                    for (int a = 0; a < ND; a++) { sp[i] += B[(size_t)a * dim + i] * s[a]; cp[i] += B[(size_t)a * dim + i] * c[a]; }
                    for (int j = 0; j < dim; j++) {
                        double v = 0;
                        for (int a = 0; a < ND; a++)
                            for (int b = 0; b < ND; b++) v += B[(size_t)a * dim + i] * M[(size_t)a * ND + b] * B[(size_t)b * dim + j];
                        Mp[(size_t)i * dim + j] = v;
                    }
                }
                std::vector<double> mu, U;
                small_eigen(Mp, dim, mu, U);
                std::vector<double> s2(dim, 0.0), c2(dim, 0.0);
                for (int k = 0; k < dim; k++) {
                    mu[k] = std::max(mu[k], 0.0);
                    for (int i = 0; i < dim; i++) { s2[k] += U[(size_t)i * dim + k] * sp[i]; c2[k] += U[(size_t)i * dim + k] * cp[i]; }
                }
                // the step for alpha (and the trust region): z_k = (alpha s_k - c_k) / (alpha + mu_k + beta)
                if (adapt && C_prev > 0 && pred_prev > 0 && C_prev > pred_prev) {
                    const double ratio = (C_prev - C0) / (C_prev - pred_prev);
                    if (ratio > 0.75 && limited_prev) rcur = std::min(rcur * 2, 100.0);
                    else if (ratio < 0.25) rcur = std::max(rcur / 2, 0.01);
                }
                const double r2 = rcur * sum_f;
                std::vector<double> z(dim);
                auto zfor = [&](double al, double be) {
                    double d2 = 0;
                    for (int k = 0; k < dim; k++) { z[k] = (al * s2[k] - c2[k]) / (al + mu[k] + be); d2 += z[k] * z[k]; }
                    return d2;
                };
                auto step_for = [&](double al) {   // z at alpha, inside the trust region
                    if (zfor(al, 0.0) <= r2) return;
                    double blo = 0, bhi = std::max(1e-300, al);
                    for (int t = 0; t < 400 && zfor(al, bhi) > r2; t++) bhi *= 2;
                    for (int t = 0; t < 100; t++) {
                        double bm = 0.5 * (blo + bhi);
                        if (zfor(al, bm) > r2) blo = bm; else bhi = bm;
                    }
                    zfor(al, bhi);
                };
                auto predC = [&]() {
                    double v = C0;
                    for (int k = 0; k < dim; k++) v += c2[k] * z[k] + 0.5 * mu[k] * z[k] * z[k];
                    return v;
                };
                double ns = 0, nc = 0;
                for (int k = 0; k < dim; k++) { ns += s2[k] * s2[k]; nc += c2[k] * c2[k]; }
                // the scale of alpha (at the start grad S is 0: f is the default everywhere)
                const double scale = ns > 1e-30 * nc && ns > 0 ? std::sqrt(nc / ns) : std::sqrt(nc) + 1.0;
                const double alo = scale * 1e-12, ahi = scale * 1e12;
                step_for(alo);
                const double Cmin = predC();
                const bool on_target = chi2n <= target_cur;
                const double frac = push_left > 0 ? 0.0 : final_pass ? aimF : aim1;
                const double aim = on_target ? C0 : Cmin + (C0 - Cmin) * frac;
                double llo = std::log(alo), lhi = std::log(ahi);
                for (int t = 0; t < 120; t++) {   // C rises with alpha
                    double lm = 0.5 * (llo + lhi);
                    step_for(std::exp(lm));
                    if (predC() > aim) lhi = lm; else llo = lm;
                }
                const double al = std::exp(llo);
                step_for(al);
                {
                    double d2 = 0;
                    for (int k = 0; k < dim; k++) d2 += z[k] * z[k];
                    limited_prev = d2 >= 0.99 * r2;
                    pred_prev = predC();
                }
                // the step in f: coefficients on the 4 directions
                std::vector<double> x(ND, 0.0);
                for (int a = 0; a < ND; a++)
                    for (int i = 0; i < dim; i++)
                        for (int k = 0; k < dim; k++) x[a] += B[(size_t)a * dim + i] * U[(size_t)i * dim + k] * z[k];
                for (int a = 0; a < ND; a++) x[a] *= dsc[a];   // on the directions as they are in memory
                for (double& v : x) if (!std::isfinite(v)) v = 0;   // never a step that is not a number
                if (reject) { f_keep = f; fit_keep = fitv; C_keep = C0; have_keep = true; }
                double clipped = 0;   // masses held at a limit (then the fit is not the old fit plus the step's)
                psums(na, 1, &clipped, [&](size_t i0, size_t i1, double* acc) {
                    double cnt = 0;
                    for (size_t j = i0; j < i1; j++) {
                        double v = f[j];
                        for (int a = 0; a < ND; a++) v += x[a] * e[a][j];
                        const double lin = v;
                        if (logstep && v < f[j]) v = f[j] * std::exp(std::max(-30.0, (v - f[j]) / f[j]));
                        v = std::max(v, 0.1 * f[j]);   // stays positive
                        f[j] = std::min(fhi, std::max(flo, v));
                        u[j] = std::log(f[j]);
                        cnt += f[j] != lin;
                    }
                    acc[0] = cnt;
                });
                if (clipped == 0 && !getenv("MS_MEMSYS_NOLIN")) {   // the new fit: the old one plus the step's (A is linear)
                    pfor(ng, [&](size_t a0, size_t b0) {
                        for (size_t i = a0; i < b0; i++) {
                            double v = fitv[i];
                            for (int a = 0; a < ND; a++) v += x[a] * Ae[a][i];
                            fitv[i] = v;
                        }
                    });
                    fit_ok = true;
                } else {
                    fit_ok = false;
                }
                n_clip += clipped > 0;
                history.push_back({alpha_eff, chi2n});
                alpha = alpha_eff;
                if (debug)
                    fprintf(stderr, "memsys pass %d step %d: C %.6g chi2n %.4f test %.3g alpha %.3g dim %d aim %.6g Cmin %.6g\n",
                            pass + 1, it + 1, C0, chi2n, test, alpha_eff, dim, aim, Cmin);
                ms::tick(cb, user, std::min(rounds - 1, total_it / 10), rounds);
                // the end of the pass: C stopped falling (best fit), or on target, with the step settled
                if (C_prev > 0 && (inner ? since_refit >= 2 : it >= 5)) flat_n = (C_prev - C0) / C0 < (final_pass ? flat_tol : flat_tol1) ? flat_n + 1 : 0;
                C_prev = C0;
                C_last = C0;
                if (push_left > 0) {
                    if (--push_left == 0) break;
                    continue;
                }
                if (inner) {
                    const bool settled = ((on_target || flat_n >= 3) && test < test_tol) || since_refit + 1 >= inner_max;
                    if (settled) {
                        const bool same = C_settled > 0 && std::fabs(C_settled - C0) / C0 < (final_pass ? outer_tol : 5 * outer_tol);
                        if (debug) fprintf(stderr, "memsys settled after %d steps: C %.6g (before the last refit %.6g)\n", since_refit + 1, C0, C_settled);
                        if (same) {
                            why = on_target ? "target" : "best fit";
                            if (!final_pass && push_n > 0) { push_left = push_n; continue; }
                            break;
                        }
                        C_settled = C0;
                        do_refit = true;
                    }
                    continue;
                }
                if (on_target && test < test_tol) { why = "target"; break; }
                if (flat_n >= 3 && test < test_tol) { why = "best fit"; break; }
            }
            sync_f();
            // polish: L-BFGS (as in 4.33's rounds) at the alpha the path found, with an envelope refit after each
            // sweep: the fine, weak structure that the few directions of the path converge on slowly
            {
                const int polish1 = getenv("MS_MEMSYS_POLISH") ? atoi(getenv("MS_MEMSYS_POLISH")) : 50;
                // (the final pass: none by default; the polish matters for the noise estimated after the first pass)
                const int polish = !final_pass ? polish1 : getenv("MS_MEMSYS_POLISH_FINAL") ? atoi(getenv("MS_MEMSYS_POLISH_FINAL")) : 0;
                const int sweeps = getenv("MS_MEMSYS_POLISH_SWEEPS") ? atoi(getenv("MS_MEMSYS_POLISH_SWEEPS")) : 1;
                // first pass: at a lower alpha (the fit pushed as far as 4.33's rounds push it before the noise is
                // estimated from its residuals)
                const double pa1 = getenv("MS_MEMSYS_PA1") ? atof(getenv("MS_MEMSYS_PA1")) : 1.0;
                if (polish > 0) {
                    obj.alpha = final_pass ? alpha : alpha * pa1;
                    sync_u();
                    for (int sw = 0; sw < sweeps; sw++) {
                        lbfgs_b(obj, u, lo_u, hi_u, polish);
                        pfor(na, [&](size_t a0, size_t b0) { for (size_t j = a0; j < b0; j++) { u[j] = std::min(hi_u, std::max(lo_u, u[j])); f[j] = fexp(u[j]); } });
                        refit_envelope();
                    }
                    sync_f();
                    if (debug) fprintf(stderr, "memsys pass %d polished at alpha %.4g (%zu evaluations so far)\n", pass + 1, alpha, obj.evals);
                }
            }
            model.forward(f.data(), fitv.data());
            if (debug) fprintf(stderr, "memsys pass %d ends (%s) after %d steps in all (%d with masses at a limit, %d rejected), chi2n %.4f\n", pass + 1, why.c_str(), total_it, n_clip, n_rej, chi2_signal());
            if (final_pass) {
                reached = why == "target";
                stop_reason = why == "target" ? "target" : why == "best fit" ? "unattainable" : "rounds";
                if (why == "best fit") plateau = chi2_signal();
                break;
            }
            // between the passes: the noise from the residuals, the envelope nodes at the species
            std::vector<double> r(ng);
            for (size_t i = 0; i < ng; i++) r[i] = fitv[i] - y[i];
            std::string dbg;
            {
                NoiseModel nm = noise_from_residuals(yvar, r, signal, noise_data, noise_kind, debug ? &dbg : nullptr);
            if (kappa_mode) nm = kappa_noise(yvar, r);
            if (opt().fixnoise) nm = noise;   // the noise measured on the data kept (no re-estimation; 4.36)
                if (debug) fprintf(stderr, "  noise re-estimation: s0 %.3g -> %.3g, 1/g %.3g -> %.3g, eps %.3g -> %.3g\n%s",
                                   noise.s0, nm.s0, noise.ginv, nm.ginv, noise.eps, nm.eps, dbg.c_str());
                noise = nm;
            }
            if (p.envelope_nodes <= 0) {
                std::vector<double> nodes = species_nodes();
                if (nodes != model.node_mass && nodes.size() >= 2) { model.set_nodes(nodes, masses); seed_envelopes(); }
            }
            set_weights();
        }
        (void)C_last;
    };
    // EXPERIMENT (MS_ME_CLASSIC): the final pass stops at the alpha of "Classic" maximum entropy with automatic noise
    // scaling (Gull and Skilling, MemSys5 users' manual, sections 2.1-2.6): the evidence for alpha and for a common
    // scale c of the noise is largest where -2 alpha S = G c^2 and c^2 = chi2 / (N - G), i.e. where
    // Omega = -2 alpha S (N - G) / (G chi2) = 1. G = sum lambda / (lambda + alpha) is the number of good measurements
    // (lambda: the eigenvalues of the curvature of chi2 / 2 in the entropy metric, B = f^1/2 A'WA f^1/2): estimated as
    // na - alpha r'(alpha I + B)^-1 r with one fixed random +-1 vector r and conjugate gradients (each step one forward
    // and one adjoint product). A model that cannot reach chi2 = N (peak shape, adducts) raises c instead of
    // driving alpha to the best-fit plateau. The noise model's shape stays the one measured on the data.
    const bool classic = getenv("MS_ME_CLASSIC") != nullptr;
    const int classic_cg = getenv("MS_ME_CLASSIC_CG") ? atoi(getenv("MS_ME_CLASSIC_CG")) : 100;
    double last_G = 0, last_c2 = 0;
    auto evidence_omega = [&]() {
        const double S = psum(na, [&](size_t j) { const double fj = f[j]; return fj - dflt - fj * (std::log(std::max(fj, 1e-300)) - log_dflt); });
        double chi2 = 0, Nd = 0;
        {
            double acc[2];
            psums(ng, 2, acc, [&](size_t i0, size_t i1, double* a) {
                double c = 0, n = 0;
                for (size_t i = i0; i < i1; i++) if (w[i] > 0) { const double d = fitv[i] - y[i]; c += w[i] * d * d; n += 1; }
                a[0] = c; a[1] = n;
            });
            chi2 = acc[0]; Nd = acc[1];
        }
        // conjugate gradients on (alpha I + B) x = r
        std::vector<double> rv(na), x(na, 0.0), res(na), pv(na), q(na), sq(na), t(na), tg(ng);
        for (size_t j = 0; j < na; j++) {
            uint64_t hsh = (uint64_t)j * 0x9E3779B97F4A7C15ull;
            hsh ^= hsh >> 29; hsh *= 0xBF58476D1CE4E5B9ull; hsh ^= hsh >> 32;
            rv[j] = (hsh & 1) ? 1.0 : -1.0;
            sq[j] = std::sqrt(std::max(f[j], 0.0));
        }
        auto applyA = [&](const std::vector<double>& v, std::vector<double>& out) {
            pfor(na, [&](size_t a0, size_t b0) { for (size_t j = a0; j < b0; j++) t[j] = sq[j] * v[j]; });
            model.forward(t.data(), tg.data());
            pfor(ng, [&](size_t a0, size_t b0) { for (size_t i = a0; i < b0; i++) tg[i] *= w[i]; });
            model.adjoint_f(tg.data(), out.data());
            pfor(na, [&](size_t a0, size_t b0) { for (size_t j = a0; j < b0; j++) out[j] = alpha * v[j] + sq[j] * out[j]; });
        };
        res = rv; pv = rv;
        double rr = psum(na, [&](size_t j) { return res[j] * res[j]; });
        const double rr0 = rr;
        int it = 0;
        for (; it < classic_cg && rr > 1e-8 * rr0; it++) {
            applyA(pv, q);
            const double pq = psum(na, [&](size_t j) { return pv[j] * q[j]; });
            if (!(pq > 0)) break;
            const double a_ = rr / pq;
            pfor(na, [&](size_t a0, size_t b0) { for (size_t j = a0; j < b0; j++) { x[j] += a_ * pv[j]; res[j] -= a_ * q[j]; } });
            const double rr1 = psum(na, [&](size_t j) { return res[j] * res[j]; });
            const double b_ = rr1 / rr;
            rr = rr1;
            pfor(na, [&](size_t a0, size_t b0) { for (size_t j = a0; j < b0; j++) pv[j] = res[j] + b_ * pv[j]; });
        }
        const double rx = psum(na, [&](size_t j) { return rv[j] * x[j]; });
        const double G = std::max(1e-6, std::min((double)na - alpha * rx, Nd - 1.0));
        const double c2 = chi2 / std::max(Nd - G, 1.0);
        const double om = -2.0 * alpha * S * (Nd - G) / std::max(G * chi2, 1e-300);
        last_G = G; last_c2 = c2;
        if (debug) fprintf(stderr, "  classic: S %.6g chi2 %.6g N %.0f G %.1f (cg %d, residual %.1e) c^2 %.4g Omega %.4f\n",
                           S, chi2, Nd, G, it, std::sqrt(rr / rr0), c2, om);
        return om;
    };
    ms::tick(cb, user, 0, rounds);
    if (memsys && !hybrid) memsys_run();
    else for (int k = 0; k < rounds; k++) {
        alpha = std::exp(la);
        obj.alpha = alpha;
        if (have_fit) set_weights();
        for (int sweep = 0; sweep < 2; sweep++) {
            lbfgs_b(obj, u, lo_u, hi_u, iterations);
            pfor(na, [&](size_t a0, size_t b0) { for (size_t j = a0; j < b0; j++) { u[j] = std::min(hi_u, std::max(lo_u, u[j])); f[j] = fexp(u[j]); } });
            refit_envelope();
        }
        pfor(na, [&](size_t a0, size_t b0) { for (size_t j = a0; j < b0; j++) f[j] = fexp(u[j]); });
        model.forward(f.data(), fitv.data());
        have_fit = true;
        double chi2n = chi2_signal();
        double lc = std::log(std::max(chi2n, 1e-300) / target_cur);
        if (classic && pass + 1 >= max_passes) lc = std::log(std::max(evidence_omega(), 1e-300));   // the evidence target
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
            if (kappa_mode) nm = kappa_noise(yvar, r);
            if (opt().fixnoise) nm = noise;   // the noise measured on the data kept (no re-estimation; 4.36)
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
                    seed_envelopes();
                    renoded = true;
                    if (debug) { fprintf(stderr, "  envelope nodes now at:"); for (double m : nodes) fprintf(stderr, " %.1f", m); fprintf(stderr, "\n"); }
                }
            }
            noise = nm;
            set_weights();
            search.reset(); pass_pts.clear();
            if (hybrid) {   // the final pass by the path
                if (!getenv("MS_MEMSYS_NOWIN")) model.make_windows();
                memsys_start = pass;
                memsys_run();
                break;
            }
            if (renoded || classic) continue;   // with new nodes the fit changes a lot at the same alpha: measure it first
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
        if (classic) {   // Omega = 1 (lc = ln Omega): no plateau, no cap
            if ((on_target && enough) || last_round || (stalled && enough)) {
                reached = on_target;
                stop_reason = on_target ? "evidence" : stalled ? "stalled" : "rounds";
                if (debug) fprintf(stderr, "  classic stop (%s): alpha %.6g G %.1f c^2 %.4g\n", stop_reason.c_str(), alpha, last_G, last_c2);
                break;
            }
            la = search.next();
            continue;
        }
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
    Timed tout(T_OUTPUT);
    pfor(na, [&](size_t a0, size_t b0) { for (size_t j = a0; j < b0; j++) f[j] = fexp(u[j]); });
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
    // the place of each output mass between the envelope nodes, once (env_at would find it again for every charge)
    std::vector<int> o_q(nmg);
    std::vector<double> o_t(nmg);
    parallel_range(nmg, [&](size_t a0, size_t b0) { for (size_t i = a0; i < b0; i++) model.node_of(mgrid[i], o_q[i], o_t[i]); });
    // EXPERIMENT (MS_ME_AREA): each charge's profile divided by its charge on the way to the mass axis. A peak of
    // width w in m/z is z w wide in mass, so heights copied over make a species' area the sum over its charges of
    // z x its m/z area (the more charges, the larger: +20 % on two proteins at mean charges 14 and 17, +70 % at 12
    // and 20 against the m/z signal they have). Divided by z the area is the m/z signal summed over the charges;
    // MS_ME_AREA=2 divides by z^2 (the number of ions where the signal grows with the charge, as in an Orbitrap).
    const int area_mode = opt().area;
    for (int k = 0; k < K; k++) {
        double z = charges[k];
        const double zdiv = area_mode == 1 ? z : area_mode == 2 ? z * z : 1.0;
        const double* col = cols.data() + (size_t)k * ng;
        parallel_range(nmg, [&](size_t a0, size_t b0) {
            for (size_t i = a0; i < b0; i++) {
                double pos = (mgrid[i] + z * a) / z;
                double xq = use_log ? std::log(std::max(pos, 1e-9)) : pos;
                const double* c0 = model.C.data() + (size_t)o_q[i] * K;   // as env_at(mgrid[i], k)
                h[i] = interp0(xq, coord.data(), col, ng) * (c0[k] + o_t[i] * (c0[K + k] - c0[k])) * scale / zdiv;
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

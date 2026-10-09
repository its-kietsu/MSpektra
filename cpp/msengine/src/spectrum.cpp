// Spectrum tools: centroiding, peak apex, resolving power, log m/z bins.
// Mirrors hrms_data.centroid / _LogBins and ms_deconv.refine_peak_mass /
// estimate_resolution of the Python reference (same arithmetic, same order
// of the floating point operations where it matters).
#include "common.h"
#include <algorithm>
#include <numeric>
#include <limits>

namespace ms {

// ------------------------------------------------------------- centroid
void centroid(const double* mz, const double* it, long n, double rel,
              std::vector<double>& cmz, std::vector<double>& cit) {
    cmz.clear();
    cit.clear();
    if (n < 3) return;
    double top = it[0];
    for (long i = 1; i < n; i++) if (it[i] > top) top = it[i];
    if (!(top > 0)) return;
    const double thr = rel * top;
    for (long k = 1; k < n - 1; k++) {
        const double a = it[k - 1], b = it[k], c = it[k + 1];
        if (!(b > a && b >= c && b > thr)) continue;
        double off = 0.0;
        if (a > 0 && c > 0) {
            const double la = std::log(std::max(a, 1e-12)), lb = std::log(b), lc = std::log(std::max(c, 1e-12));
            const double den = la - 2 * lb + lc;
            if (den < 0) off = 0.5 * (la - lc) / den;
        }
        if (off < -0.5) off = -0.5;
        if (off > 0.5) off = 0.5;
        const double step = off >= 0 ? mz[k + 1] - mz[k] : mz[k] - mz[k - 1];
        cmz.push_back(mz[k] + off * step);
        cit.push_back(b);
    }
}

// ----------------------------------------------------------------- apex
// weighted least squares fit of c2 x^2 + c1 x + c0 to log y (np.polyfit with
// w: columns scaled by their norms, solved by Householder QR). false when
// the system is degenerate.
static bool polyfit2(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& w,
                     double& c2, double& c1) {
    const size_t m = x.size();
    if (m < 3) return false;
    std::vector<double> A(3 * m), b(m);   // column major: A[j*m + i]
    for (size_t i = 0; i < m; i++) {
        A[0 * m + i] = x[i] * x[i] * w[i];
        A[1 * m + i] = x[i] * w[i];
        A[2 * m + i] = w[i];
        b[i] = y[i] * w[i];
    }
    double scale[3];
    for (int j = 0; j < 3; j++) {
        double s = 0;
        for (size_t i = 0; i < m; i++) s += A[j * m + i] * A[j * m + i];
        scale[j] = std::sqrt(s);
        if (!(scale[j] > 0)) return false;
        for (size_t i = 0; i < m; i++) A[j * m + i] /= scale[j];
    }
    // Householder QR
    double R[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (int j = 0; j < 3; j++) {
        double norm = 0;
        for (size_t i = j; i < m; i++) norm += A[j * m + i] * A[j * m + i];
        norm = std::sqrt(norm);
        if (!(norm > 0)) return false;
        const double alpha = A[j * m + j] > 0 ? -norm : norm;
        std::vector<double> v(m, 0.0);
        for (size_t i = j; i < m; i++) v[i] = A[j * m + i];
        v[j] -= alpha;
        double vv = 0;
        for (size_t i = j; i < m; i++) vv += v[i] * v[i];
        if (vv > 0) {
            for (int k = j; k < 3; k++) {
                double d = 0;
                for (size_t i = j; i < m; i++) d += v[i] * A[k * m + i];
                d = 2 * d / vv;
                for (size_t i = j; i < m; i++) A[k * m + i] -= d * v[i];
            }
            double d = 0;
            for (size_t i = j; i < m; i++) d += v[i] * b[i];
            d = 2 * d / vv;
            for (size_t i = j; i < m; i++) b[i] -= d * v[i];
        }
        for (int k = j; k < 3; k++) R[j][k] = A[k * m + j];
    }
    double c[3];
    for (int j = 2; j >= 0; j--) {
        double s = b[j];
        for (int k = j + 1; k < 3; k++) s -= R[j][k] * c[k];
        if (!(std::fabs(R[j][j]) > 0)) return false;
        c[j] = s / R[j][j];
    }
    c2 = c[0] / scale[0];
    c1 = c[1] / scale[1];
    return std::isfinite(c2) && std::isfinite(c1);
}

double apex(const double* x, const double* y, long n, double m0) {
    if (n < 3) return m0;
    long k = 0;
    double best = std::fabs(x[0] - m0);
    for (long i = 1; i < n; i++) {
        const double d = std::fabs(x[i] - m0);
        if (d < best) { best = d; k = i; }
    }
    const double top = y[k];
    if (top <= 0 || k == 0 || k == n - 1 || std::max(y[k - 1], y[k + 1]) > top) return m0;
    long a = k, b = k;
    while (a > 0 && y[a - 1] >= 0.5 * top && y[a - 1] <= y[a]) a--;
    while (b < n - 1 && y[b + 1] >= 0.5 * top && y[b + 1] <= y[b]) b++;
    if (b - a >= 2) {
        std::vector<double> xx, ly, w;
        for (long i = a; i <= b; i++) {
            xx.push_back(x[i] - x[k]);
            ly.push_back(std::log(std::max(y[i], 1e-300)));
            w.push_back(y[i] / top);
        }
        double c2 = 0, c1 = 0;
        if (polyfit2(xx, ly, w, c2, c1) && c2 < 0 && std::fabs(c1 / (2 * c2)) <= (x[b] - x[a]))
            return x[k] - c1 / (2 * c2);
    }
    if (y[k - 1] > 0 && y[k + 1] > 0) {
        const double la = std::log(y[k - 1]), lb = std::log(top), lc = std::log(y[k + 1]);
        const double den = la - 2 * lb + lc;
        if (den < 0) {
            const double off = std::max(-0.5, std::min(0.5, 0.5 * (la - lc) / den));
            const double step = off >= 0 ? x[k + 1] - x[k] : x[k] - x[k - 1];
            return x[k] + off * step;
        }
    }
    return x[k];
}

// ----------------------------------------------------------- resolution
// FWHM of the peak at i, half height crossings interpolated; NaN when the
// peak is not enclosed by points below half height.
static double fwhm_at(const double* x, const double* y, long n, long i) {
    const double half = y[i] / 2.0;
    long lo = i, hi = i;
    while (lo > 0 && y[lo] > half) lo--;
    while (hi < n - 1 && y[hi] > half) hi++;
    if (y[lo] > half || y[hi] > half || lo == i || hi == i) return std::numeric_limits<double>::quiet_NaN();
    const double xl = x[lo] + (half - y[lo]) * (x[lo + 1] - x[lo]) / std::max(y[lo + 1] - y[lo], 1e-30);
    const double xr = x[hi - 1] + (y[hi - 1] - half) * (x[hi] - x[hi - 1]) / std::max(y[hi - 1] - y[hi], 1e-30);
    return xr - xl;
}

static double median(std::vector<double> v) {
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

double resolution(const double* x, const double* y, long n) {
    if (n < 5) return 10000.0;
    long imax = 0;
    for (long i = 1; i < n; i++) if (y[i] > y[imax]) imax = i;
    if (y[imax] <= 0) return 10000.0;
    // the 5000 tallest points, tallest first
    std::vector<long> order(n);
    std::iota(order.begin(), order.end(), 0L);
    const size_t m = std::min<size_t>(5000, n);
    auto taller = [&](long p, long q) { return y[p] > y[q] || (y[p] == y[q] && p < q); };
    std::partial_sort(order.begin(), order.begin() + m, order.end(), taller);
    std::vector<double> used, vals;
    const double floor = 0.05 * y[order[0]];
    for (size_t q = 0; q < m; q++) {
        const long i = order[q];
        if (y[i] < floor || vals.size() >= 5) break;
        bool near = false;
        for (double u : used) if (std::fabs(x[i] - u) < 0.2) { near = true; break; }
        if (near) continue;
        used.push_back(x[i]);
        const double w = fwhm_at(x, y, n, i);
        if (std::isfinite(w) && w > 0) vals.push_back(x[i] / w);
    }
    if (!vals.empty()) return median(vals);
    // estimate_peak_width: the tallest peak, else 5 x the median spacing
    double w = fwhm_at(x, y, n, imax);
    if (!(std::isfinite(w) && w > 0)) {
        std::vector<double> d(n - 1);
        for (long i = 0; i + 1 < n; i++) d[i] = x[i + 1] - x[i];
        w = 5 * median(d);
    }
    return x[imax] / std::max(w, 1e-9);
}

// -------------------------------------------------------------- LogBins
void LogBins::add(const double* mz, const double* it, long n) {
    // keys of the points with m/z > 0
    thread_local std::vector<long> keys;
    thread_local std::vector<long> idx;
    keys.clear();
    idx.clear();
    long lo = 0, hi = 0;
    for (long i = 0; i < n; i++) {
        if (!(mz[i] > 0) || !std::isfinite(mz[i])) continue;
        const double kd = std::rint(std::log(mz[i]) / r_);
        if (!(std::fabs(kd) < 2e9)) continue;   // ((long) of a larger value is undefined; long is 32 bit on Windows)
        const long k = (long)kd;
        if (keys.empty()) lo = hi = k;
        else { lo = std::min(lo, k); hi = std::max(hi, k); }
        keys.push_back(k);
        idx.push_back(i);
    }
    if (keys.empty()) return;
    {   // a span no spectrum has (m/z 1e-300 to 1e300 from a damaged file) would ask for terabytes
        const long long top = empty_ ? hi : std::max<long long>(hi, (long long)k0_ + (long long)si_.size() - 1);
        const long long bot = empty_ ? lo : std::min<long long>(lo, k0_);
        if (top - bot >= 100000000LL) throw Error("The m/z values of the spectra span an implausible range (damaged data?)");
    }
    if (empty_) {
        k0_ = lo;
        si_.assign(hi - lo + 1, 0.0);
        sw_.assign(hi - lo + 1, 0.0);
        sm_.assign(hi - lo + 1, 0.0);
        empty_ = false;
    } else if (lo < k0_ || hi > k0_ + (long)si_.size() - 1) {
        const long new0 = std::min(lo, k0_);
        const long nn = std::max(hi, k0_ + (long)si_.size() - 1) - new0 + 1;
        const long off = k0_ - new0;
        for (std::vector<double>* v : {&si_, &sw_, &sm_}) {
            std::vector<double> a(nn, 0.0);
            std::copy(v->begin(), v->end(), a.begin() + off);
            v->swap(a);
        }
        k0_ = new0;
    }
    // this spectrum's bincounts first, then added to the running sums (as
    // numpy does it: the rounding then matches the reference). With the
    // keys in order (m/z sorted) the points of one bin are consecutive and
    // no temporary arrays are needed.
    const long nb = hi - lo + 1;
    const long base = lo - k0_;
    bool sorted = true;
    for (size_t j = 1; j < keys.size(); j++) if (keys[j] < keys[j - 1]) { sorted = false; break; }
    if (sorted) {
        size_t j = 0;
        while (j < keys.size()) {
            const long key = keys[j];
            double ti = 0, tw = 0, tm = 0;
            for (; j < keys.size() && keys[j] == key; j++) {
                const double v = it[idx[j]], w = std::max(v, 1e-30);
                ti += v;
                tw += w;
                tm += mz[idx[j]] * w;
            }
            si_[base + key - lo] += ti;
            sw_[base + key - lo] += tw;
            sm_[base + key - lo] += tm;
        }
        return;
    }
    thread_local std::vector<double> ti, tw, tm;
    ti.assign(nb, 0.0);
    tw.assign(nb, 0.0);
    tm.assign(nb, 0.0);
    for (size_t j = 0; j < keys.size(); j++) {
        const long b = keys[j] - lo;
        const double v = it[idx[j]], w = std::max(v, 1e-30);
        ti[b] += v;
        tw[b] += w;
        tm[b] += mz[idx[j]] * w;
    }
    for (long b = 0; b < nb; b++) {
        si_[base + b] += ti[b];
        sw_[base + b] += tw[b];
        sm_[base + b] += tm[b];
    }
}

void LogBins::result(std::vector<double>& mz, std::vector<double>& it) const {
    mz.clear();
    it.clear();
    if (empty_) return;
    mz.reserve(si_.size());
    it.reserve(si_.size());
    for (size_t b = 0; b < si_.size(); b++) {
        if (sw_[b] > 0) {
            mz.push_back(sm_[b] / sw_[b]);
            it.push_back(si_[b]);
        }
    }
}

}  // namespace ms

// ------------------------------------------------------------ C wrappers
extern "C" {

MS_API int ms_centroid(const double* mz, const double* it, long n, double rel,
                       const double** cmz, const double** cit, long* nc) {
    static thread_local std::vector<double> omz, oit;
    return ms::guarded([&] {
        if (!cmz || !cit || !nc || (n > 0 && (!mz || !it))) throw std::invalid_argument("ms_centroid: null argument");
        ms::centroid(mz, it, n, rel, omz, oit);
        *cmz = omz.data();
        *cit = oit.data();
        *nc = (long)omz.size();
    });
}

MS_API double ms_apex(const double* mz, const double* it, long n, double x) {
    if (!mz || !it || n <= 0) return x;
    return ms::guarded_value(x, [&] { return ms::apex(mz, it, n, x); });
}

MS_API double ms_resolution(const double* mz, const double* it, long n) {
    if (!mz || !it || n <= 0) return 10000.0;
    return ms::guarded_value(10000.0, [&] { return ms::resolution(mz, it, n); });
}

}

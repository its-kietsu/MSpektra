// LC peak integration: port of lcms_integrate.py (smooth, auto_integrate, manual_peak,
// peak_at, split_peak, finish). The scipy pieces it relies on are reproduced step by step:
// savgol_filter(y, w, 2) with the polynomial edge fit (mode 'interp'; the coefficients of
// scipy itself, savgol_tab.h, and the summation order of scipy.ndimage.correlate1d, so the
// smoothed values away from the edges are scipy's), minimum_filter1d / uniform_filter1d
// (window [i - w/2, i + w - w/2 - 1], 'nearest', the running sum of scipy), find_peaks(height,
// prominence, width) with its prominence and half height width definitions, and numpy's
// pairwise summation of the areas. Areas in signal units x seconds (times in minutes).
#include "common.h"
#include "e2_internal.h"   // np_sum
#include "savgol_tab.h"
#include <algorithm>
#include <cfloat>
#include <deque>
#include <limits>

namespace ms {
namespace {

double median(std::vector<double> v) {
    if (v.empty()) return 0.0;
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    return n % 2 ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}

// least squares projection weights of a polynomial of degree p over the points 0 .. w-1:
// H[j * w + k] = weight of point k in the fitted value at point j (orthonormal basis in long double)
std::vector<double> hat_weights(int w, int p) {
    std::vector<std::vector<long double>> q;
    for (int d = 0; d <= p; d++) {
        std::vector<long double> v(w);
        for (int k = 0; k < w; k++) {
            long double x = 1.0L;
            for (int e = 0; e < d; e++) x *= (long double)k;
            v[k] = x;
        }
        for (int pass = 0; pass < 2; pass++) {   // modified Gram-Schmidt, twice
            for (const auto& u : q) {
                long double dot = 0;
                for (int k = 0; k < w; k++) dot += u[k] * v[k];
                for (int k = 0; k < w; k++) v[k] -= dot * u[k];
            }
        }
        long double nn = 0;
        for (int k = 0; k < w; k++) nn += v[k] * v[k];
        nn = std::sqrt(nn);
        for (int k = 0; k < w; k++) v[k] /= nn;
        q.push_back(v);
    }
    std::vector<double> H((size_t)w * w);
    for (int j = 0; j < w; j++)
        for (int k = 0; k < w; k++) {
            long double s = 0;
            for (const auto& u : q) s += u[j] * u[k];
            H[(size_t)j * w + k] = (double)s;
        }
    return H;
}

// scipy.signal.savgol_coeffs(w, p) for windows outside the table: the exact fractions
std::vector<double> savgol_coeffs(int w, int p) {
    std::vector<double> c(w);
    const int m = w / 2;
    for (int j = -m; j <= m; j++) {
        if (p >= 2) c[j + m] = (3.0 * (3.0 * m * m + 3.0 * m - 1) - 15.0 * j * j) / ((2.0 * m + 1) * (4.0 * m * m + 4.0 * m - 3));
        else c[j + m] = 1.0 / w;
    }
    return c;
}

// lcms_integrate.smooth(y, points): Savitzky-Golay of order 2 (order 1 over 3 points);
// points < 3 or fewer than points + 2 values: unchanged; an even count is made odd
std::vector<double> smooth(const std::vector<double>& y, int points) {
    const long n = (long)y.size();
    if (points < 3 || n < points + 2) return y;
    if (points % 2 == 0) points += 1;
    const int w = points, p = w >= 5 ? 2 : 1, m = w / 2;
    std::vector<double> c;
    const double* tab = sg::table(w);
    if (tab) c.assign(tab, tab + w);
    else c = savgol_coeffs(w, p);
    // convolve1d(y, c, mode='constant') = correlate1d with the reversed weights; the
    // summation of ndimage's NI_Correlate1D: symmetric weights add the outermost pair first
    std::vector<double> fw(c.rbegin(), c.rend());
    bool symmetric = true;
    for (int ii = 1; ii <= m; ii++)
        if (std::fabs(fw[m + ii] - fw[m - ii]) > DBL_EPSILON) { symmetric = false; break; }
    std::vector<double> out(n, 0.0);
    for (long i = m; i < n - m; i++) {
        double s;
        if (symmetric) {
            s = y[i] * fw[m];
            for (int jj = -m; jj < 0; jj++) s += (y[i + jj] + y[i - jj]) * fw[m + jj];
        } else {
            s = y[i + m] * fw[m + m];
            for (int jj = -m; jj < m; jj++) s += y[i + jj] * fw[m + jj];
        }
        out[i] = s;
    }
    // edges: the polynomial through the first / last w points (savgol's _fit_edge)
    if (w == 5) {
        // np.polyfit + polyval written as the projection weights (exact fractions / 35)
        out[0] = (31 * y[0] + 9 * y[1] - 3 * y[2] - 5 * y[3] + 3 * y[4]) / 35.0;
        out[1] = (9 * y[0] + 13 * y[1] + 12 * y[2] + 6 * y[3] - 5 * y[4]) / 35.0;
        out[n - 2] = (9 * y[n - 1] + 13 * y[n - 2] + 12 * y[n - 3] + 6 * y[n - 4] - 5 * y[n - 5]) / 35.0;
        out[n - 1] = (31 * y[n - 1] + 9 * y[n - 2] - 3 * y[n - 3] - 5 * y[n - 4] + 3 * y[n - 5]) / 35.0;
    } else {
        const std::vector<double> H = hat_weights(w, p);
        for (int j = 0; j < m; j++) {
            double s = 0, t = 0;
            for (int k = 0; k < w; k++) {
                s += H[(size_t)j * w + k] * y[k];
                t += H[(size_t)(w - m + j) * w + k] * y[n - w + k];
            }
            out[j] = s;
            out[n - m + j] = t;
        }
    }
    return out;
}

// _envelope: rolling minimum then rolling mean, both of width w with the
// ndimage window [i - w/2, i + w - w/2 - 1] and 'nearest' extension
std::vector<double> envelope(const std::vector<double>& y, double window) {
    const long n = (long)y.size();
    const long w = (long)std::max(3.0, std::min(window, (double)n));
    const long s1 = w / 2, s2 = w - s1 - 1;
    auto at = [&](long i) { return y[std::max(0L, std::min(n - 1, i))]; };
    std::vector<double> m(n);
    // rolling minimum with a monotonic deque over the extended line
    std::deque<long> dq;   // indices into the extended line (-s1 .. n-1+s2)
    long next = -s1;
    for (long i = 0; i < n; i++) {
        const long hi = i + s2;
        for (; next <= hi; next++) {
            const double v = at(next);
            while (!dq.empty() && at(dq.back()) >= v) dq.pop_back();
            dq.push_back(next);
        }
        while (dq.front() < i - s1) dq.pop_front();
        m[i] = at(dq.front());
    }
    std::vector<double> out(n);
    auto mat = [&](long i) { return m[std::max(0L, std::min(n - 1, i))]; };
    double tmp = 0.0;
    for (long k = 0; k < w; k++) tmp += mat(k - s1);
    out[0] = tmp / w;
    for (long i = 1; i < n; i++) {
        tmp += mat(i + s2) - mat(i - 1 - s1);
        out[i] = tmp / w;
    }
    return out;
}

double noise_of(const double* y, long n) {
    if (n - 1 < 4) return 0.0;
    std::vector<double> d(n - 1);
    for (long i = 0; i + 1 < n; i++) d[i] = y[i + 1] - y[i];
    const double md = median(d);
    for (double& v : d) v = std::fabs(v - md);
    return 1.4826 * median(d) / std::sqrt(2.0);
}

struct Peak { ms_lc_peak2 p; };

// _make_peak: baseline b0 at i0, b1 at i1; height above it, area in y*s (numpy's pairwise sum)
bool make_peak(const double* t, const double* y, long i0, long i1, double b0, double b1, Peak& out) {
    if (i1 - i0 + 1 < 2) return false;
    const double span = std::max(t[i1] - t[i0], 1e-12);
    long k = i0;
    double dk = 0, prev = 0;
    std::vector<double> prod;
    prod.reserve(i1 - i0);
    for (long i = i0; i <= i1; i++) {
        const double base = b0 + (b1 - b0) * (t[i] - t[i0]) / span;
        double d = y[i] - base;
        if (d < 0) d = 0;
        if (i == i0) dk = d;
        else if (!std::isnan(dk) && (std::isnan(d) || d > dk)) { dk = d; k = i; }   // np.argmax: first maximum
        if (i > i0) prod.push_back((d + prev) * (t[i] * 60.0 - t[i - 1] * 60.0));
        prev = d;
    }
    out.p.rt = t[k];
    out.p.t0 = t[i0];
    out.p.t1 = t[i1];
    out.p.height = dk;
    out.p.area = ms::e2::np_sum(prod.data(), prod.size()) / 2.0;
    out.p.base0 = b0;
    out.p.base1 = b1;
    out.p.apex_y = y[k];
    out.p.i0 = i0;
    out.p.i1 = i1;
    return true;
}

// scipy.signal.find_peaks(x, height=thr, prominence=thr, width=wmin)
std::vector<long> find_peaks(const std::vector<double>& x, double thr, double wmin) {
    const long n = (long)x.size();
    std::vector<long> pk;
    long i = 1;
    const long imax = n - 1;
    while (i < imax) {
        if (x[i - 1] < x[i]) {
            long ahead = i + 1;
            while (ahead < imax && x[ahead] == x[i]) ahead++;
            if (x[ahead] < x[i]) {
                pk.push_back((i + ahead - 1) / 2);
                i = ahead;
            }
        }
        i++;
    }
    std::vector<long> out;
    for (long p : pk) {
        if (!(x[p] >= thr)) continue;           // height
        // prominence: lowest point between the peak and the next higher one on either side
        long lb = p, rb = p;
        double lmin = x[p], rmin = x[p];
        for (long j = p; j >= 0 && x[j] <= x[p]; j--) if (x[j] < lmin) { lmin = x[j]; lb = j; }
        for (long j = p; j <= n - 1 && x[j] <= x[p]; j++) if (x[j] < rmin) { rmin = x[j]; rb = j; }
        const double prom = x[p] - std::max(lmin, rmin);
        if (!(prom >= thr)) continue;
        // width at half prominence (interpolated crossings)
        const double h = x[p] - prom * 0.5;
        long j = p;
        while (lb < j && h < x[j]) j--;
        double lip = (double)j;
        if (x[j] < h) lip += (h - x[j]) / (x[j + 1] - x[j]);
        j = p;
        while (j < rb && h < x[j]) j++;
        double rip = (double)j;
        if (x[j] < h) rip -= (h - x[j]) / (x[j - 1] - x[j]);
        if (!(rip - lip >= wmin)) continue;
        out.push_back(p);
    }
    return out;
}

void finish(std::vector<Peak>& peaks) {
    std::stable_sort(peaks.begin(), peaks.end(), [](const Peak& a, const Peak& b) { return a.p.rt < b.p.rt; });
}

// the borders of the peak at p: back at the baseline or at the valley before the next peak
void borders(const std::vector<double>& d, long p, long lo_i, long hi_i, double lim, double noise, long& i0, long& i1) {
    long i = p;
    while (i > lo_i) {
        if (d[i] <= lim) break;
        if (d[i - 1] > d[i] && d[std::max(i - 3, lo_i)] > d[i] + noise) break;
        i--;
    }
    long j = p;
    while (j < hi_i) {
        if (d[j] <= lim) break;
        if (d[j + 1] > d[j] && d[std::min(j + 3, hi_i)] > d[j] + noise) break;
        j++;
    }
    i0 = i;
    i1 = j;
}

// auto_integrate(t, y, threshold_pct, min_width_s, smooth_points, t_range, baseline_window_min);
// the peaks in the order they are made (finish() sorts them)
std::vector<Peak> auto_integrate(const double* t, const double* y, long n, double threshold_pct, double min_width_s,
                                 int smooth_points = 0, bool use_range = false, double ta = 0, double tb = 0,
                                 double bwin = 1.5) {
    std::vector<Peak> peaks;
    if (n < 8) return peaks;
    const std::vector<double> yv(y, y + n);
    const std::vector<double> ys = smooth(yv, smooth_points ? smooth_points : 5);
    std::vector<double> dt(n - 1);
    for (long i = 0; i + 1 < n; i++) dt[i] = t[i + 1] - t[i];
    const double dt_s = median(dt) * 60.0;
    const std::vector<double> env = envelope(ys, bwin * 60.0 / std::max(dt_s, 1e-6));
    std::vector<double> d(n);
    for (long i = 0; i < n; i++) d[i] = ys[i] - env[i];
    long lo_i = 0, hi_i = n - 1;
    if (use_range) {
        const double a = std::min(ta, tb), b = std::max(ta, tb);
        long first = -1, last = -1, cnt = 0;
        for (long i = 0; i < n; i++)
            if (t[i] >= a && t[i] <= b) { if (first < 0) first = i; last = i; cnt++; }
        if (cnt < 8) return peaks;
        lo_i = first;
        hi_i = last;
    }
    const std::vector<double> seg(d.begin() + lo_i, d.begin() + hi_i + 1);
    double top = seg[0];
    for (double v : seg) {   // seg.max(): a NaN wins
        if (std::isnan(v)) { top = v; break; }
        if (v > top) top = v;
    }
    if (top <= 0) return peaks;
    const double noise = noise_of(y + lo_i, hi_i - lo_i + 1);
    const double thr = std::max(top * threshold_pct / 100.0, 5.0 * noise);
    const double wmin = std::max(1.0, min_width_s / dt_s);
    std::vector<long> pk = find_peaks(seg, thr, wmin);
    if (pk.empty()) return peaks;
    for (long& p : pk) p += lo_i;
    const double stop_level = std::max(noise * 2.0, 0.0);
    std::vector<std::pair<long, long>> bounds;
    for (long p : pk) {
        long i0, i1;
        borders(d, p, lo_i, hi_i, std::max(stop_level, 0.005 * d[p]), noise, i0, i1);
        bounds.push_back({i0, i1});
    }
    // neighbouring peaks meet at the valley of the smoothed signal
    for (size_t q = 1; q < bounds.size(); q++) {
        auto& a = bounds[q - 1];
        auto& b = bounds[q];
        if (b.first <= a.second) {
            long v = pk[q - 1];
            for (long i = pk[q - 1]; i <= pk[q]; i++) if (ys[i] < ys[v]) v = i;
            a.second = v;
            b.first = v;
        }
    }
    // clusters of touching peaks share one baseline (drop lines at the valleys)
    size_t q = 0;
    while (q < bounds.size()) {
        size_t m = q;
        while (m + 1 < bounds.size() && bounds[m + 1].first <= bounds[m].second) m++;
        const long c0 = bounds[q].first, c1 = bounds[m].second;
        const double y0 = ys[c0], y1 = ys[c1];
        const double span = std::max(t[c1] - t[c0], 1e-12);
        for (size_t r = q; r <= m; r++) {
            const long i0 = bounds[r].first, i1 = bounds[r].second;
            if (i1 - i0 < 2) continue;
            const double f0 = (t[i0] - t[c0]) / span;
            const double f1 = (t[i1] - t[c0]) / span;
            Peak pkd;
            if (make_peak(t, y, i0, i1, y0 + (y1 - y0) * f0, y0 + (y1 - y0) * f1, pkd) && pkd.p.height > 0)
                peaks.push_back(pkd);
        }
        q = m + 1;
    }
    return peaks;
}

bool manual_peak(const double* t, const double* y, long n, double t0, double t1, Peak& out) {
    const double a = std::min(t0, t1), b = std::max(t0, t1);
    long i0 = -1, i1 = -1, cnt = 0;
    for (long i = 0; i < n; i++) {
        if (t[i] >= a && t[i] <= b) {
            if (i0 < 0) i0 = i;
            i1 = i;
            cnt++;
        }
    }
    if (cnt < 3) return false;
    const std::vector<double> ys = smooth(std::vector<double>(y, y + n), 5);
    return make_peak(t, y, i0, i1, ys[i0], ys[i1], out);
}

// first index of the smallest |t - x| (np.argmin)
long nearest(const double* t, long n, double x) {
    long k = 0;
    double best = std::fabs(t[0] - x);
    for (long i = 1; i < n && !std::isnan(best); i++) {
        const double d = std::fabs(t[i] - x);
        if (std::isnan(d) || d < best) { best = d; k = i; }
    }
    return k;
}

// peak_at: the peak whose top is nearest to x (within +- search_s seconds)
bool peak_at(const double* t, const double* y, long n, double x, double search_s, double bwin, Peak& out) {
    if (n < 8) return false;
    const std::vector<double> ys = smooth(std::vector<double>(y, y + n), 5);
    std::vector<double> dt(n - 1);
    for (long i = 0; i + 1 < n; i++) dt[i] = t[i + 1] - t[i];
    const double dt_s = median(dt) * 60.0;
    const std::vector<double> env = envelope(ys, bwin * 60.0 / std::max(dt_s, 1e-6));
    std::vector<double> d(n);
    for (long i = 0; i < n; i++) d[i] = ys[i] - env[i];
    const double noise = noise_of(y, n);
    const long i = nearest(t, n, x);
    const long w = std::max(2L, (long)std::nearbyint(search_s / std::max(dt_s, 1e-6)));
    const long lo = std::max(0L, i - w), hi = std::min(n - 1, i + w);
    long p = lo;
    for (long j = lo + 1; j <= hi; j++) {
        if (std::isnan(d[p])) break;
        if (std::isnan(d[j]) || d[j] > d[p]) p = j;
    }
    const double md = median(d);
    std::vector<double> ad(n);
    for (long j = 0; j < n; j++) ad[j] = std::fabs(d[j] - md);
    const double spread = 1.4826 * median(ad);
    if (d[p] <= std::max(4.0 * noise, 3.0 * spread) || d[p] <= 0) return false;
    long i0, i1;
    borders(d, p, 0, n - 1, std::max(2.0 * noise, 0.005 * d[p]), noise, i0, i1);
    if (i1 - i0 < 2) return false;
    return make_peak(t, y, i0, i1, ys[i0], ys[i1], out);
}

bool split_peak(const double* t, const double* y, long n, long i0, long i1, double b0, double b1, double x,
                Peak& a, Peak& b) {
    const long k = nearest(t, n, x);
    if (!(i0 + 1 < k && k < i1 - 1)) return false;
    const double ta = t[i0], tb = t[i1];
    const double bm = b0 + (b1 - b0) * (t[k] - ta) / std::max(tb - ta, 1e-12);
    return make_peak(t, y, i0, k, b0, bm, a) && make_peak(t, y, k, i1, bm, b1, b);
}

ms_lc_peak to_old(const ms_lc_peak2& p) {
    ms_lc_peak o;
    o.t0 = p.t0; o.t1 = p.t1; o.rt = p.rt; o.height = p.height; o.area = p.area; o.base0 = p.base0; o.base1 = p.base1;
    return o;
}

}  // namespace
}  // namespace ms

extern "C" {

MS_API int ms_integrate_auto(const double* t, const double* y, long n, double threshold_pct, double min_width_s,
                             const ms_lc_peak** peaks, long* np) {
    static thread_local std::vector<ms_lc_peak> out;
    return ms::guarded([&] {
        if (!peaks || !np || (n > 0 && (!t || !y))) throw std::invalid_argument("ms_integrate_auto: null argument");
        std::vector<ms::Peak> pk = ms::auto_integrate(t, y, n, threshold_pct, min_width_s);
        ms::finish(pk);
        out.clear();
        for (const ms::Peak& p : pk) out.push_back(ms::to_old(p.p));
        *peaks = out.data();
        *np = (long)out.size();
    });
}

MS_API int ms_integrate_range(const double* t, const double* y, long n, double t0, double t1, ms_lc_peak* out) {
    bool few = false;
    const int rc = ms::guarded([&] {   // (no exception, e.g. bad_alloc, may leave the C interface)
        if (!t || !y || !out || n <= 0) throw std::invalid_argument("ms_integrate_range: null argument");
        ms::Peak p;
        if (!ms::manual_peak(t, y, n, t0, t1, p)) { few = true; return; }
        *out = ms::to_old(p.p);
    });
    if (rc == MS_OK && few) { ms::set_error("the range holds fewer than 3 points"); return MS_NOT_FOUND; }
    return rc;
}

// ---- 3.1 additions (agent P) --------------------------------------------------------

MS_API int ms_smooth(const double* y, long n, int points, double* out) {
    return ms::guarded([&] {
        if (n < 0 || (n > 0 && (!y || !out))) throw std::invalid_argument("ms_smooth: bad arguments");
        if (points > 100001) throw std::invalid_argument("ms_smooth: at most 100001 points");
        const std::vector<double> r = ms::smooth(std::vector<double>(y, y + n), points);
        std::copy(r.begin(), r.end(), out);
    });
}

MS_API int ms_lc_integrate(const double* t, const double* y, long n, double threshold_pct, double min_width_s,
                           int smooth_points, int use_range, double t_lo, double t_hi, double baseline_window_min,
                           const ms_lc_peak2** peaks, long* np) {
    static thread_local std::vector<ms_lc_peak2> out;
    return ms::guarded([&] {
        if (!peaks || !np || n < 0 || (n > 0 && (!t || !y))) throw std::invalid_argument("ms_lc_integrate: bad arguments");
        if (smooth_points > 100001) throw std::invalid_argument("ms_lc_integrate: at most 100001 smoothing points");
        if (use_range && (std::isnan(t_lo) || std::isnan(t_hi))) throw std::invalid_argument("ms_lc_integrate: NaN time range");
        const std::vector<ms::Peak> pk = ms::auto_integrate(t, y, n, threshold_pct, min_width_s, smooth_points,
                                                            use_range != 0, t_lo, t_hi, baseline_window_min);
        out.clear();
        for (const ms::Peak& p : pk) out.push_back(p.p);
        *peaks = out.data();
        *np = (long)out.size();
    });
}

MS_API int ms_lc_manual(const double* t, const double* y, long n, double t0, double t1, ms_lc_peak2* out) {
    bool none = false;
    const int rc = ms::guarded([&] {
        if (!out || n < 0 || (n > 0 && (!t || !y))) throw std::invalid_argument("ms_lc_manual: bad arguments");
        ms::Peak p;
        if (n == 0 || !ms::manual_peak(t, y, n, t0, t1, p)) { none = true; return; }
        *out = p.p;
    });
    return (rc == MS_OK && none) ? MS_NOT_FOUND : rc;
}

MS_API int ms_lc_peak_at(const double* t, const double* y, long n, double x, double search_s, double baseline_window_min,
                         ms_lc_peak2* out) {
    bool none = false;
    const int rc = ms::guarded([&] {
        if (!out || n < 0 || (n > 0 && (!t || !y))) throw std::invalid_argument("ms_lc_peak_at: bad arguments");
        ms::Peak p;
        if (!ms::peak_at(t, y, n, x, search_s, baseline_window_min, p)) { none = true; return; }
        *out = p.p;
    });
    return (rc == MS_OK && none) ? MS_NOT_FOUND : rc;
}

MS_API int ms_lc_split(const double* t, const double* y, long n, long i0, long i1, double b0, double b1, double x,
                       ms_lc_peak2* out2) {
    bool none = false;
    const int rc = ms::guarded([&] {
        if (!out2 || n < 0 || (n > 0 && (!t || !y))) throw std::invalid_argument("ms_lc_split: bad arguments");
        if (i0 < 0 || i1 >= n || i0 > i1) throw std::invalid_argument("ms_lc_split: peak borders out of range");
        ms::Peak a, b;
        if (!ms::split_peak(t, y, n, i0, i1, b0, b1, x, a, b)) { none = true; return; }
        out2[0] = a.p;
        out2[1] = b.p;
    });
    return (rc == MS_OK && none) ? MS_NOT_FOUND : rc;
}

}

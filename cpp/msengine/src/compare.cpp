// Compare view of LCMS Analysis: port of the calculations of lcms_compare_core.py (subtract_blank,
// find_apex, baseline, process, stack, the numbers of peak_rows) and of the lambda max of
// lcms_compare.CompareTab.on_lambda_max. The smoothing, the rolling baseline, the integration and
// the PDA traces are the library's own (ms_smooth, ms_subtract_baseline, ms_lc_integrate,
// ms_pda_*), called as the Python code calls them. numpy's rules are kept where they decide a
// result: np.median, np.nanmin / np.nanmax (NaN ignored, all NaN gives NaN), np.argmax /
// np.argmin (the first NaN wins), np.interp, Python's max() / min() (the first of equal values,
// NaN never replaces), round() half to even and the compensated float sum() of Python 3.12.
#include "common.h"
#include <algorithm>
#include <cfenv>
#include <limits>
#include <stdexcept>
#include <vector>

namespace ms {
namespace {

const double NaN = std::numeric_limits<double>::quiet_NaN();

void raise_if(int rc) {
    // an inner library call failed: its message is already set (ms_last_error); keep it
    if (rc != MS_OK) throw std::runtime_error(ms_last_error());
}

double np_median(std::vector<double> v) {
    if (v.empty()) return NaN;
    for (double x : v)
        if (std::isnan(x)) return NaN;
    std::sort(v.begin(), v.end());
    const size_t n = v.size();
    return n % 2 ? v[n / 2] : (v[n / 2 - 1] + v[n / 2]) / 2.0;
}

double nanmin(const double* y, long n) {
    double m = NaN;
    for (long i = 0; i < n; i++)
        if (!std::isnan(y[i]) && (std::isnan(m) || y[i] < m)) m = y[i];
    return m;
}

double nanmax(const double* y, long n) {
    double m = NaN;
    for (long i = 0; i < n; i++)
        if (!std::isnan(y[i]) && (std::isnan(m) || y[i] > m)) m = y[i];
    return m;
}

// np.max / np.min of y[lo:hi] (a NaN gives NaN)
double np_max(const double* y, long lo, long hi) {
    double m = y[lo];
    for (long i = lo; i < hi; i++) {
        if (std::isnan(y[i])) return NaN;
        if (y[i] > m) m = y[i];
    }
    return m;
}

double np_min(const double* y, long lo, long hi) {
    double m = y[lo];
    for (long i = lo; i < hi; i++) {
        if (std::isnan(y[i])) return NaN;
        if (y[i] < m) m = y[i];
    }
    return m;
}

// Python's max(a, b) / min(a, b): b only when it compares greater / smaller
inline double py_max(double a, double b) { return b > a ? b : a; }

// np.argmax / np.argmin over idx (the first NaN wins, else the first extreme)
long argext(const double* y, const std::vector<long>& idx, bool want_max) {
    long k = idx[0];
    if (std::isnan(y[k])) return k;
    for (size_t q = 1; q < idx.size(); q++) {
        const double v = y[idx[q]];
        if (std::isnan(v)) return idx[q];
        if (want_max ? v > y[k] : v < y[k]) k = idx[q];
    }
    return k;
}

double py_round(double x) {
    // round() of Python: half to even (the default rounding mode of nearbyint)
    const int old = std::fegetround();
    std::fesetround(FE_TONEAREST);
    const double r = std::nearbyint(x);
    std::fesetround(old);
    return r;
}

// sum() of floats of Python 3.12 (Neumaier compensation)
double py_sum(const std::vector<double>& v) {
    double s = 0.0, c = 0.0;
    for (double x : v) {
        const double t = s + x;
        if (std::fabs(s) >= std::fabs(x)) c += (s - t) + x;
        else c += (x - t) + s;
        s = t;
    }
    return s + c;
}

std::vector<double> smoothed(const double* y, long n, int points) {
    std::vector<double> out(y, y + n);
    if (points >= 3 && n > 0) raise_if(ms_smooth(y, n, points, out.data()));
    return out;
}

// find_apex(t, y, tc, win)
double find_apex(const double* t, const double* y, long n, double tc, double win) {
    std::vector<long> m;
    for (long i = 0; i < n; i++)
        if (t[i] >= tc - win && t[i] <= tc + win) m.push_back(i);
    if (m.size() < 3) return NaN;
    double dt = 0.0;
    if (n > 1) {
        std::vector<double> d(n - 1);
        for (long i = 0; i + 1 < n; i++) d[i] = t[i + 1] - t[i];
        dt = np_median(d);
    }
    long k = 2;
    if (dt > 0) k = std::max(2L, (long)py_round(0.03 / dt));
    double noise = 0.0;
    if (n > 2) {
        std::vector<double> d(n - 1);
        for (long i = 0; i + 1 < n; i++) d[i] = std::fabs(y[i + 1] - y[i]);
        noise = 1.4826 * np_median(d) / std::sqrt(2.0);
    }
    // the lowest point of the window left of i (y[m0 .. i]) and right of it (y[i .. m_last])
    const long m0 = m.front(), m1 = m.back();
    std::vector<double> pre(m1 - m0 + 1), suf(m1 - m0 + 1);
    for (long i = m0; i <= m1; i++) {
        const double v = y[i];
        const double p = i == m0 ? v : pre[i - 1 - m0];
        pre[i - m0] = (std::isnan(p) || std::isnan(v)) ? NaN : std::min(p, v);
    }
    for (long i = m1; i >= m0; i--) {
        const double v = y[i];
        const double s = i == m1 ? v : suf[i + 1 - m0];
        suf[i - m0] = (std::isnan(s) || std::isnan(v)) ? NaN : std::min(s, v);
    }
    std::vector<long> cand;
    for (size_t q = 1; q + 1 < m.size(); q++) {
        const long i = m[q];
        const long lo = std::max(0L, i - k), hi = std::min(n, i + k + 1);
        if (y[i] >= np_max(y, lo, hi) && y[i] > np_min(y, lo, hi)) {
            const double rise = y[i] - py_max(pre[i - m0], suf[i - m0]);
            if (rise > 10.0 * noise) cand.push_back(i);
        }
    }
    if (cand.empty()) return NaN;
    long i = cand[0];
    for (size_t q = 1; q < cand.size(); q++)
        if (y[cand[q]] > y[i]) i = cand[q];
    if (0 < i && i < n - 1) {
        const double y0 = y[i - 1], y1 = y[i], y2 = y[i + 1];
        const double d = y0 - 2 * y1 + y2;
        if (d < 0) {
            const double off = 0.5 * (y0 - y2) / d;
            if (std::fabs(off) <= 1) return t[i] + off * (t[i + 1] - t[i - 1]) / 2.0;
        }
    }
    return t[i];
}

// baseline(t, y, kind, rolling_min) in place: 0 none, 1 offset, 2 drift, 3 rolling
void baseline(const double* t, std::vector<double>& y, int kind, double rolling_min) {
    const long n = (long)y.size();
    if (kind == 1 && n) {
        const double lo = nanmin(y.data(), n);
        for (double& v : y) v = v - lo;
    } else if (kind == 2 && n > 10) {
        const long k = std::max(2L, n / 33);
        const double a = np_median(std::vector<double>(y.begin(), y.begin() + k));
        const double b = np_median(std::vector<double>(y.end() - k, y.end()));
        const double ta = np_median(std::vector<double>(t, t + k));
        const double tb = np_median(std::vector<double>(t + n - k, t + n));
        double den = tb - ta;
        if (den == 0.0) den = 1.0;   // (tb - ta) or 1.0
        for (long i = 0; i < n; i++) y[i] = y[i] - (a + (b - a) * (t[i] - ta) / den);
    } else if (kind == 3 && n > 20) {
        std::vector<double> out(n);
        raise_if(ms_subtract_baseline(t, y.data(), n, py_max(rolling_min, 0.01), 15.0, out.data()));
        y.swap(out);
    }
}

void check_xy(const double* t, const double* y, long n, const char* fn) {
    if (n < 0 || (n > 0 && (!t || !y))) throw std::invalid_argument(std::string(fn) + ": bad arguments");
}

// np.interp(x, xp, fp) at one point (xp rising, n >= 2; left / right: the end values)
double np_interp(double x, const double* xp, const double* fp, long n) {
    if (std::isnan(x)) return x;
    if (x < xp[0]) return fp[0];
    if (x > xp[n - 1]) return fp[n - 1];
    const long j = (long)(std::upper_bound(xp, xp + n, x) - xp) - 1;
    if (j >= n - 1 || xp[j] == x) return fp[j];
    const double slope = (fp[j + 1] - fp[j]) / (xp[j + 1] - xp[j]);
    double r = slope * (x - xp[j]) + fp[j];
    if (std::isnan(r)) {
        r = slope * (x - xp[j + 1]) + fp[j + 1];
        if (std::isnan(r) && fp[j] == fp[j + 1]) r = fp[j];
    }
    return r;
}

// np.sum of a contiguous float64 array (numpy's pairwise summation)
double np_pairwise_sum(const double* a, long n) {
    if (n < 8) {
        double res = 0.0;
        for (long i = 0; i < n; i++) res += a[i];
        return res;
    }
    if (n <= 128) {
        double r[8];
        for (int j = 0; j < 8; j++) r[j] = a[j];
        long i = 8;
        for (; i < n - (n % 8); i += 8)
            for (int j = 0; j < 8; j++) r[j] += a[i + j];
        double res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        for (; i < n; i++) res += a[i];
        return res;
    }
    long n2 = n / 2;
    n2 -= n2 % 8;
    return np_pairwise_sum(a, n2) + np_pairwise_sum(a + n2, n - n2);
}

// _region_peak(t, y, a, b): the samples from a to b (both ends interpolated), the straight baseline between
// the signal at the ends, the points where the signal crosses it inserted. false: no data there (fewer than
// 2 points, a value that is not finite, time not rising, or the region outside the run).
bool region_peak(const double* t, const double* y, long n, double a, double b, std::vector<double>& ts,
                 std::vector<double>& ys, std::vector<double>& base, double* res) {
    if (n < 2) return false;
    for (long i = 0; i < n; i++)
        if (!std::isfinite(t[i]) || !std::isfinite(y[i])) return false;
    for (long i = 0; i + 1 < n; i++)
        if (t[i + 1] - t[i] <= 0) return false;
    const double lo = t[0] > a ? t[0] : a;          // max(a, t[0]) and min(b, t[-1]) of Python
    const double hi = t[n - 1] < b ? t[n - 1] : b;
    if (!std::isfinite(lo) || !std::isfinite(hi) || hi <= lo) return false;
    ts.clear();
    ts.push_back(lo);
    for (long i = 0; i < n; i++)
        if (t[i] > lo && t[i] < hi) ts.push_back(t[i]);
    ts.push_back(hi);
    auto fill = [&]() {
        const size_t m = ts.size();
        ys.resize(m);
        base.resize(m);
        for (size_t i = 0; i < m; i++) ys[i] = np_interp(ts[i], t, y, n);
        for (size_t i = 0; i < m; i++) base[i] = ys[0] + (ys[m - 1] - ys[0]) * (ts[i] - lo) / (hi - lo);
    };
    fill();
    // the points where the signal crosses the baseline (the partial trapezoids exact)
    std::vector<double> zeros;
    for (size_t i = 0; i + 1 < ts.size(); i++) {
        const double d0 = ys[i] - base[i], d1 = ys[i + 1] - base[i + 1];
        if (d0 * d1 < 0) zeros.push_back(ts[i] - d0 * (ts[i + 1] - ts[i]) / (d1 - d0));
    }
    if (!zeros.empty()) {
        ts.insert(ts.end(), zeros.begin(), zeros.end());
        std::sort(ts.begin(), ts.end());
        fill();
    }
    const long m = (long)ts.size();
    std::vector<double> d(m);
    for (long i = 0; i < m; i++) {
        const double v = ys[i] - base[i];
        d[i] = std::isnan(v) ? v : (v > 0.0 ? v : 0.0);   // np.maximum(v, 0.0)
    }
    std::vector<long> all(m);
    for (long i = 0; i < m; i++) all[i] = i;
    const long k = argext(d.data(), all, true);
    std::vector<double> tr(m - 1);
    for (long i = 0; i + 1 < m; i++) tr[i] = (d[i + 1] + d[i]) * (ts[i + 1] - ts[i]);
    res[0] = lo;
    res[1] = hi;
    res[2] = ts[k];
    res[3] = d[k];
    res[4] = np_pairwise_sum(tr.data(), m - 1) * 30.0;
    res[5] = base[0];
    res[6] = base[m - 1];
    return true;
}

}  // namespace
}  // namespace ms

extern "C" {

MS_API int ms_cmp_subtract_blank(const double* t, const double* y, long n, const double* bt, const double* by,
                                 long nb, double* out) {
    return ms::guarded([&] {
        ms::check_xy(t, y, n, "ms_cmp_subtract_blank");
        if (!out || nb < 0 || (nb > 0 && (!bt || !by))) throw std::invalid_argument("ms_cmp_subtract_blank: bad arguments");
        if (nb < 2) { std::copy(y, y + n, out); return; }
        for (long i = 0; i < n; i++) {
            const double x = t[i];
            double b = ms::NaN;   // np.interp(t, bt, by, left=NaN, right=NaN)
            if (x >= bt[0] && x <= bt[nb - 1]) {
                if (x == bt[nb - 1]) b = by[nb - 1];
                else {
                    const long j = (long)(std::upper_bound(bt, bt + nb, x) - bt) - 1;
                    if (x == bt[j]) b = by[j];
                    else {
                        const double slope = (by[j + 1] - by[j]) / (bt[j + 1] - bt[j]);
                        b = slope * (x - bt[j]) + by[j];
                        if (std::isnan(b)) {   // numpy: an infinite slope at a defined end
                            b = slope * (x - bt[j + 1]) + by[j + 1];
                            if (std::isnan(b) && by[j] == by[j + 1]) b = by[j];
                        }
                    }
                }
            }
            out[i] = y[i] - (std::isfinite(b) ? b : 0.0);
        }
    });
}

MS_API int ms_cmp_find_apex(const double* t, const double* y, long n, double tc, double win, int smooth_points,
                            double* apex) {
    return ms::guarded([&] {
        ms::check_xy(t, y, n, "ms_cmp_find_apex");
        if (!apex) throw std::invalid_argument("ms_cmp_find_apex: null argument");
        if (smooth_points > 100001) throw std::invalid_argument("ms_cmp_find_apex: at most 100001 smoothing points");
        *apex = ms::NaN;
        if (n == 0) return;
        const std::vector<double> ys = ms::smoothed(y, n, smooth_points);
        *apex = ms::find_apex(t, ys.data(), n, tc, win);
    });
}

MS_API int ms_cmp_baseline(const double* t, const double* y, long n, int kind, double rolling_min, double* out) {
    return ms::guarded([&] {
        ms::check_xy(t, y, n, "ms_cmp_baseline");
        if (!out && n > 0) throw std::invalid_argument("ms_cmp_baseline: null argument");
        std::vector<double> v(y, y + n);
        ms::baseline(t, v, kind, rolling_min);
        std::copy(v.begin(), v.end(), out);
    });
}

MS_API int ms_cmp_process(const double* t, const double* y, long n, const ms_cmp_opts* o, double* t_out,
                          double* y_out, long* n_out, double* factor, double* top, int* ref_missing) {
    return ms::guarded([&] {
        ms::check_xy(t, y, n, "ms_cmp_process");
        if (!o || !n_out || !factor || !top || !ref_missing || (n > 0 && (!t_out || !y_out)))
            throw std::invalid_argument("ms_cmp_process: null argument");
        if (o->smooth > 100001) throw std::invalid_argument("ms_cmp_process: at most 100001 smoothing points");
        const std::vector<double> ys = ms::smoothed(y, n, o->smooth);
        std::vector<double> tw, yw;
        for (long i = 0; i < n; i++) {
            const double ti = t[i] + o->shift;
            if (!std::isnan(o->t0) && !(ti >= o->t0)) continue;
            if (!std::isnan(o->t1) && !(ti <= o->t1)) continue;
            tw.push_back(ti);
            yw.push_back(ys[i]);
        }
        const long m = (long)tw.size();
        ms::baseline(tw.data(), yw, o->baseline, o->rolling_min);
        double f = 1.0;
        if (o->scale == 1 && m) {
            const double tp = ms::nanmax(yw.data(), m);
            f = tp > 0 ? 100.0 / tp : 1.0;
        }
        int missing = 0;
        if (o->scale == 2 && m && !std::isnan(o->ref_t)) {
            const double win = ms::py_max(o->ref_win, 0.02);
            std::vector<double> sel;
            for (long i = 0; i < m; i++)
                if (tw[i] >= o->ref_t - win && tw[i] <= o->ref_t + win) sel.push_back(yw[i]);
            const double tp = sel.empty() ? 0.0 : ms::nanmax(sel.data(), (long)sel.size());
            f = tp > 0 ? 100.0 / tp : 1.0;
            missing = !(tp > 0);
        }
        for (double& v : yw) v = v * f;
        std::copy(tw.begin(), tw.end(), t_out);
        std::copy(yw.begin(), yw.end(), y_out);
        *n_out = m;
        *factor = f;
        *top = m ? ms::nanmax(yw.data(), m) : 0.0;
        *ref_missing = missing;
    });
}

MS_API int ms_cmp_stack(const double* tops, const double* spans, const int* ok, long n, double spacing,
                        int reverse, double skew, int stacked, double* dx, double* dy) {
    return ms::guarded([&] {
        if (n < 0 || (n > 0 && (!tops || !spans || !ok || !dx || !dy)))
            throw std::invalid_argument("ms_cmp_stack: bad arguments");
        if (!stacked || n == 0) {
            std::fill(dx, dx + n, 0.0);
            std::fill(dy, dy + n, 0.0);
            return;
        }
        // max([top of the traces drawn] + [0.0]) and max(spans) of Python: the first of equal values
        bool first = true;
        double h = 0.0, span = 0.0;
        bool any = false;
        for (long i = 0; i < n; i++) {
            if (!ok[i]) continue;
            if (first) { h = tops[i]; first = false; }
            else h = ms::py_max(h, tops[i]);
            span = any ? ms::py_max(span, spans[i]) : spans[i];
            any = true;
        }
        h = first ? 0.0 : ms::py_max(h, 0.0);
        const double step = h * spacing / 100.0;
        const double dxs = (any ? span : 0.0) * skew / 100.0;
        for (long i = 0; i < n; i++) {
            const long level = reverse ? i : n - 1 - i;
            dx[i] = level * dxs;
            dy[i] = level * step;
        }
    });
}

MS_API int ms_cmp_peaks(const double* t, const double* y, long n, double factor, double thr_pct,
                        const double* guides, long ng, const ms_lc_peak2** peaks, long* np, long* main,
                        double* main_pct, double* guide_pct) {
    static thread_local std::vector<ms_lc_peak2> out;
    return ms::guarded([&] {
        ms::check_xy(t, y, n, "ms_cmp_peaks");
        if (!peaks || !np || !main || !main_pct || ng < 0 || (ng > 0 && (!guides || !guide_pct)))
            throw std::invalid_argument("ms_cmp_peaks: bad arguments");
        out.clear();
        *main = -1;
        *main_pct = ms::NaN;
        std::fill(guide_pct, guide_pct + ng, ms::NaN);
        if (n >= 8) {
            // the signal itself: y drawn / its scale (factor 0 counts as 1, as `or 1.0`)
            const double f = factor != 0.0 ? factor : 1.0;
            std::vector<double> sig(n);
            for (long i = 0; i < n; i++) sig[i] = y[i] / f;
            const ms_lc_peak2* pk = nullptr;
            long k = 0;
            ms::raise_if(ms_lc_integrate(t, sig.data(), n, thr_pct, 2.0, 0, 0, 0.0, 0.0, 1.5, &pk, &k));
            out.assign(pk, pk + k);
            std::stable_sort(out.begin(), out.end(), [](const ms_lc_peak2& a, const ms_lc_peak2& b) { return a.rt < b.rt; });
        }
        *peaks = out.data();
        *np = (long)out.size();
        if (out.empty()) return;
        std::vector<double> areas;
        for (const auto& p : out) areas.push_back(p.area);
        const double tot = ms::py_sum(areas);
        long mi = 0;
        for (long q = 1; q < (long)out.size(); q++)
            if (out[q].height > out[mi].height) mi = q;
        *main = mi;
        if (tot > 0) *main_pct = 100.0 * out[mi].area / tot;
        for (long g = 0; g < ng; g++) {
            const double x = guides[g];
            long hit = -1;
            for (long q = 0; q < (long)out.size() && hit < 0; q++)
                if (out[q].t0 <= x && x <= out[q].t1) hit = q;
            if (hit < 0) {
                long near = 0;
                for (long q = 1; q < (long)out.size(); q++)
                    if (std::fabs(out[q].rt - x) < std::fabs(out[near].rt - x)) near = q;
                if (std::fabs(out[near].rt - x) <= 0.1) hit = near;
            }
            if (hit >= 0 && tot > 0) guide_pct[g] = 100.0 * out[hit].area / tot;
        }
    });
}

MS_API int ms_cmp_lambda_max(const float* a, long nt, long nw, const double* wl, const double* times,
                             double interval_s, double wl0, double bw, int smooth, double shift, double t0,
                             double t1, double* lmax, double* t_apex, int* bg_used) {
    bool none = false;
    const int rc = ms::guarded([&] {
        if (!a || nt <= 0 || nw <= 0 || !wl || !times || !lmax || !t_apex || !bg_used)
            throw std::invalid_argument("ms_cmp_lambda_max: bad arguments");
        if (smooth > 100001) throw std::invalid_argument("ms_cmp_lambda_max: at most 100001 smoothing points");
        // the trace shown (the PDA chromatogram, else the max plot), smoothed over at least 5 points
        std::vector<float> f(nt);
        if (!std::isnan(wl0)) ms::raise_if(ms_pda_chromatogram(a, nt, nw, wl, wl0, std::max(bw, 0.0), f.data()));
        else ms::raise_if(ms_pda_max_plot(a, nt, nw, f.data()));
        std::vector<double> y(f.begin(), f.end());
        y = ms::smoothed(y.data(), nt, std::max(5, smooth));
        // its tallest point in the time window, after the first 8 % of the run (the injection peak)
        std::vector<long> m;
        const double skip = times[0] + 0.08 * (times[nt - 1] - times[0]);
        for (long i = 0; i < nt; i++) {
            const double ts = times[i] + shift;
            if (!std::isnan(t0) && !(ts >= t0)) continue;
            if (!std::isnan(t1) && !(ts <= t1)) continue;
            if (std::isnan(t0) && !(times[i] >= skip)) continue;
            m.push_back(i);
        }
        if (m.empty()) { none = true; return; }
        const long k = ms::argext(y.data(), m, true);
        const double ta = times[k];
        // background: the lowest point of the max plot 0.15 to 1 min before the peak
        std::vector<float> mp(nt);
        ms::raise_if(ms_pda_max_plot(a, nt, nw, mp.data()));
        const std::vector<double> ym(mp.begin(), mp.end());
        std::vector<long> before;
        for (long i = 0; i < nt; i++)
            if (times[i] >= ta - 1.0 && times[i] <= ta - 0.15) before.push_back(i);
        std::vector<double> sp(nw);
        int mode = 0;
        double b0 = 0.0, b1 = 0.0;
        if (!before.empty()) {
            const double tb = times[ms::argext(ym.data(), before, false)];
            mode = 2;
            b0 = tb - 0.02;
            b1 = tb + 0.02;
        }
        ms::raise_if(ms_pda_spectrum(a, nt, nw, times, interval_s, ta, ms::NaN, mode, b0, b1, sp.data()));
        // the highest local maximum of the spectrum above 210 nm (else its highest point)
        std::vector<double> ww, aa;
        for (long j = 0; j < nw; j++)
            if (wl[j] >= 210.0) { ww.push_back(wl[j]); aa.push_back(sp[j]); }
        if (ww.size() < 3) { none = true; return; }
        long best = -1;
        for (long i = 1; i + 1 < (long)aa.size(); i++)
            if (aa[i] >= aa[i - 1] && aa[i] > aa[i + 1] && (best < 0 || aa[i] > aa[best])) best = i;
        if (best < 0) {
            std::vector<long> all(aa.size());
            for (size_t i = 0; i < aa.size(); i++) all[i] = (long)i;
            best = ms::argext(aa.data(), all, true);
        }
        *lmax = ww[best];
        *t_apex = ta;
        *bg_used = mode == 2;
    });
    if (rc == MS_OK && none) { ms::set_error("no peak with a UV spectrum above 210 nm"); return MS_NOT_FOUND; }
    return rc;
}

MS_API int ms_cmp_region(const double* t, const double* y, long n, double a, double b, double* t_out,
                         double* y_out, double* base_out, long cap, long* n_out, double* res) {
    bool none = false;
    const int rc = ms::guarded([&] {
        ms::check_xy(t, y, n, "ms_cmp_region");
        if (!t_out || !y_out || !base_out || !n_out || !res || cap < 0)
            throw std::invalid_argument("ms_cmp_region: bad arguments");
        *n_out = 0;
        std::vector<double> ts, ys, base;
        if (!ms::region_peak(t, y, n, a, b, ts, ys, base, res)) { none = true; return; }
        if ((long)ts.size() > cap) throw std::invalid_argument("ms_cmp_region: the output holds fewer than 2 n points");
        std::copy(ts.begin(), ts.end(), t_out);
        std::copy(ys.begin(), ys.end(), y_out);
        std::copy(base.begin(), base.end(), base_out);
        *n_out = (long)ts.size();
    });
    if (rc == MS_OK && none) { ms::set_error("no data in the region"); return MS_NOT_FOUND; }
    return rc;
}

}  // extern "C"

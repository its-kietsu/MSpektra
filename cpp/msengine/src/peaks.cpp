// Peak picking of mass spectra (port of pick_mass_peaks, group_isotopes and
// refine_peak_mass from ms_deconv.py) and the small numeric helpers shared
// with maxent.cpp.
#include "e2_internal.h"
#include "f2_internal.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <cfenv>
#include <functional>
#include <deque>
#include <limits>
#include <map>

namespace ms {
namespace e2 {

// numpy's pairwise summation (blocks of 8, split above 128 elements)
double np_sum(const double* a, size_t n) {
    if (n < 8) {
        double r = 0.0;
        for (size_t i = 0; i < n; i++) r += a[i];
        return r;
    }
    if (n <= 128) {
        double r[8];
        for (int j = 0; j < 8; j++) r[j] = a[j];
        size_t i = 8;
        for (; i + 8 <= n; i += 8)
            for (int j = 0; j < 8; j++) r[j] += a[i + j];
        double res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        for (; i < n; i++) res += a[i];
        return res;
    }
    size_t n2 = n / 2;
    n2 -= n2 % 8;
    return np_sum(a, n2) + np_sum(a + n2, n - n2);
}

double median(std::vector<double> v) {
    size_t n = v.size();
    if (n == 0) return std::nan("");
    size_t h = n / 2;
    std::nth_element(v.begin(), v.begin() + h, v.end());
    double m = v[h];
    if (n % 2 == 0) {
        double lo = *std::max_element(v.begin(), v.begin() + h);
        m = 0.5 * (lo + m);
    }
    return m;
}

double median_diff(const double* x, size_t n) {
    if (n < 2) return 1.0;
    std::vector<double> d(n - 1);
    for (size_t i = 0; i + 1 < n; i++) d[i] = x[i + 1] - x[i];
    return median(std::move(d));
}

double py_round(double v) {
    double r = std::nearbyint(v);   // default rounding mode: to nearest, ties to even
    return r;
}

// indices sorted by descending y (ties: lower index first)
static std::vector<size_t> order_desc(const double* y, size_t n) {
    std::vector<size_t> o(n);
    std::iota(o.begin(), o.end(), 0);
    std::stable_sort(o.begin(), o.end(), [y](size_t a, size_t b) { return y[a] > y[b]; });
    return o;
}

std::vector<Peak> pick_mass_peaks(const double* x, const double* y, size_t n, double window, double threshold) {
    std::vector<Peak> peaks;
    if (n < 3) return peaks;
    double top = *std::max_element(y, y + n);
    if (top <= 0) return peaks;
    double dx = n > 1 ? median_diff(x, n) : 1.0;
    // ((long) of NaN or of 1e300 is undefined: the window in points is clamped first; a window beyond
    // 2 n + 2 points covers the whole spectrum for the peak and for its area alike)
    const double wd = py_round(window / std::max(dx, 1e-9));
    long w = 1;
    if (wd > 1) w = wd < 2.0 * (double)n + 2 ? (long)wd : 2 * (long)n + 2;
    std::vector<size_t> order = order_desc(y, n);
    std::vector<char> taken(n, 0);
    for (size_t i : order) {
        if (y[i] <= 0 || y[i] < threshold * top) break;
        if (taken[i]) continue;
        size_t lo = (long)i - w > 0 ? i - w : 0;
        size_t hi = std::min(n, i + w + 1);
        double mx = y[lo];
        for (size_t j = lo + 1; j < hi; j++) mx = std::max(mx, y[j]);
        if (y[i] < mx) continue;
        for (size_t j = lo; j < hi; j++) taken[j] = 1;
        long hw = w / 2;
        size_t a0 = (long)i - hw > 0 ? i - hw : 0;
        size_t a1 = std::min(n, i + hw + 1);
        double sxy = 0, sy = 0, sall = 0;
        bool any = false;
        for (size_t j = a0; j < a1; j++) {
            sall += y[j];
            if (y[j] >= 0.5 * y[i]) { sxy += x[j] * y[j]; sy += y[j]; any = true; }
        }
        Peak p;
        p.mass = any ? sxy / sy : x[i];
        p.apex = x[i];
        p.height = y[i];
        p.area = sall * dx;
        p.n_iso = 0;
        peaks.push_back(p);
    }
    std::sort(peaks.begin(), peaks.end(), [](const Peak& a, const Peak& b) { return a.mass < b.mass; });
    return peaks;
}

// standard deviation (Da) of the averagine isotope distribution (ms_deconv.averagine_sigma)
static double iso_sigma(double mass) {
    const double nr = mass / 111.1254;
    struct El { double cnt; double iso[3][2]; int k; };
    const El els[5] = {{4.9384 * nr, {{0, 0.9893}, {1.00336, 0.0107}, {0, 0}}, 2},
                       {7.7583 * nr, {{0, 0.999885}, {1.00628, 0.000115}, {0, 0}}, 2},
                       {1.3577 * nr, {{0, 0.99636}, {0.99704, 0.00364}, {0, 0}}, 2},
                       {1.4773 * nr, {{0, 0.99757}, {1.00422, 0.00038}, {2.00425, 0.00205}}, 3},
                       {0.0417 * nr, {{0, 0.9499}, {0.99939, 0.0075}, {1.99580, 0.0425}}, 3}};
    double var = 0;
    for (const El& e : els) {
        double m1 = 0, m2 = 0;
        for (int i = 0; i < e.k; i++) { m1 += e.iso[i][0] * e.iso[i][1]; m2 += e.iso[i][0] * e.iso[i][0] * e.iso[i][1]; }
        var += e.cnt * (m2 - m1 * m1);
    }
    return std::sqrt(var);
}

std::vector<Peak> group_isotopes(const double* x, const double* y, size_t n, double threshold, double spacing) {
    return group_isotopes_ranges(x, y, n, threshold, nullptr, spacing);
}

// group_isotopes, and the mass of the first and last isotope peak of each species returned
std::vector<Peak> group_isotopes_ranges(const double* x, const double* y, size_t n, double threshold,
                                        std::vector<std::pair<double, double>>* ranges, double spacing) {
    std::vector<Peak> out;
    std::vector<Peak> iso = pick_mass_peaks(x, y, n, 0.35, threshold * 0.2);
    if (iso.empty()) return out;
    // chains of peaks about one isotope spacing apart
    auto chains = [&](const std::vector<Peak>& pk) {
        std::vector<std::vector<Peak>> groups;
        std::vector<Peak> cur{pk[0]};
        for (size_t i = 1; i < pk.size(); i++) {
            const Peak& q = pk[i];
            double dm = q.mass - cur.back().mass;
            if (std::fabs(dm - spacing) <= 0.12 || dm < 0.6) cur.push_back(q);
            else { groups.push_back(cur); cur.assign(1, q); }
        }
        groups.push_back(cur);
        return groups;
    };
    // a deep dip inside a chain: two species next to each other
    auto dips = [&](const std::vector<std::vector<Peak>>& groups) {
        std::vector<std::vector<Peak>> split;
        for (const auto& g0 : groups) {
            std::vector<Peak> g = g0;
            while (g.size() >= 5) {
                std::vector<double> h(g.size());
                for (size_t j = 0; j < g.size(); j++) h[j] = g[j].height;
                long cut = -1;
                for (size_t j = 1; j + 1 < g.size(); j++) {
                    double lmax = *std::max_element(h.begin(), h.begin() + j);
                    double rmax = *std::max_element(h.begin() + j + 1, h.end());
                    if (h[j] <= h[j - 1] && h[j] <= h[j + 1] && h[j] <= 0.5 * std::min(lmax, rmax)) { cut = (long)j; break; }
                }
                if (cut < 0) break;
                split.emplace_back(g.begin(), g.begin() + cut);
                if (h[cut + 1] > h[cut]) g.erase(g.begin(), g.begin() + cut + 1);
                else g.erase(g.begin(), g.begin() + cut);
            }
            if (!g.empty()) split.push_back(g);
        }
        return split;
    };
    // each species is limited to its isotope envelope (most abundant isotope +- 2.5 averagine sigma); what
    // lies outside is grouped again: another species when it has its own maximum, dropped as the tail of
    // this one when it only falls away from it (ms_deconv.group_isotopes)
    std::vector<std::vector<Peak>> split;
    std::function<void(const std::vector<Peak>&, int)> limit = [&](const std::vector<Peak>& g, int side) {
        size_t t = 0;
        for (size_t j = 1; j < g.size(); j++) if (g[j].height > g[t].height) t = j;
        if ((side < 0 && t == g.size() - 1) || (side > 0 && t == 0)) return;
        const double sig = iso_sigma(g[t].mass);
        const double lo = g[t].mass - 2.5 * sig - 0.5, hi = g[t].mass + 2.5 * sig + 0.5;
        std::vector<Peak> inner, left, right;
        for (const Peak& q : g) (q.mass < lo ? left : q.mass > hi ? right : inner).push_back(q);
        split.push_back(inner);
        if (!left.empty()) for (const auto& part : dips(chains(left))) limit(part, -1);
        if (!right.empty()) for (const auto& part : dips(chains(right))) limit(part, 1);
    };
    for (const auto& g : dips(chains(iso))) limit(g, 0);
    std::stable_sort(split.begin(), split.end(), [](const std::vector<Peak>& a, const std::vector<Peak>& b) { return a.front().mass < b.front().mass; });
    double dx = n > 1 ? median_diff(x, n) : 1.0;
    for (auto& g : split) {
        double lo = g.front().mass - 0.5, hi = g.back().mass + 0.5;
        size_t a = std::lower_bound(x, x + n, lo) - x;     // x >= lo
        size_t b = std::upper_bound(x, x + n, hi) - x;     // x <= hi
        double sxy = 0, sy = 0;
        for (size_t j = a; j < b; j++) { sxy += x[j] * y[j]; sy += y[j]; }
        Peak p;
        p.mass = (b > a && sy > 0) ? sxy / sy : g.front().mass;
        size_t t = 0;
        for (size_t j = 1; j < g.size(); j++) if (g[j].height > g[t].height) t = j;
        p.apex = g[t].mass;
        p.height = g[t].height;
        p.area = b > a ? sy * dx : 0.0;
        p.n_iso = (int)g.size();
        out.push_back(p);
    }
    double tallest = 0;
    for (auto& q : out) tallest = std::max(tallest, q.height);
    std::vector<Peak> kept;
    for (size_t k = 0; k < out.size(); k++) {
        if (!(out[k].height >= threshold * tallest)) continue;
        kept.push_back(out[k]);
        if (ranges) ranges->push_back({split[k].front().mass, split[k].back().mass});
    }
    return kept;
}

// Householder QR least squares for small dense systems (m >= n)
void lstsq_small(const std::vector<double>& A0, const std::vector<double>& b0, int m, int n, std::vector<double>& x) {
    std::vector<double> A = A0, b = b0;
    x.assign(n, 0.0);
    for (int k = 0; k < n; k++) {
        double norm = 0;
        for (int i = k; i < m; i++) norm += A[i * n + k] * A[i * n + k];
        norm = std::sqrt(norm);
        if (norm == 0) continue;
        double alpha = A[k * n + k] > 0 ? -norm : norm;
        std::vector<double> v(m - k);
        for (int i = k; i < m; i++) v[i - k] = A[i * n + k];
        v[0] -= alpha;
        double vn = 0;
        for (double t : v) vn += t * t;
        if (vn == 0) continue;
        for (int j = k; j < n; j++) {
            double s = 0;
            for (int i = k; i < m; i++) s += v[i - k] * A[i * n + j];
            s = 2 * s / vn;
            for (int i = k; i < m; i++) A[i * n + j] -= s * v[i - k];
        }
        double s = 0;
        for (int i = k; i < m; i++) s += v[i - k] * b[i];
        s = 2 * s / vn;
        for (int i = k; i < m; i++) b[i] -= s * v[i - k];
    }
    for (int k = n - 1; k >= 0; k--) {
        double s = b[k];
        for (int j = k + 1; j < n; j++) s -= A[k * n + j] * x[j];
        double d = A[k * n + k];
        x[k] = d != 0 ? s / d : 0.0;
    }
}

double refine_peak_mass(const double* x, const double* y, size_t n, double m0) {
    if (n < 3) return m0;
    size_t k = 0;
    double best = std::fabs(x[0] - m0);
    for (size_t i = 1; i < n; i++) { double d = std::fabs(x[i] - m0); if (d < best) { best = d; k = i; } }
    double top = y[k];
    if (top <= 0 || k == 0 || k == n - 1 || std::max(y[k - 1], y[k + 1]) > top) return m0;
    size_t a = k, b = k;
    while (a > 0 && y[a - 1] >= 0.5 * top && y[a - 1] <= y[a]) a--;
    while (b < n - 1 && y[b + 1] >= 0.5 * top && y[b + 1] <= y[b]) b++;
    if (b - a >= 2) {
        // weighted parabola through log y (np.polyfit with w = y / top)
        int m = (int)(b - a + 1);
        std::vector<double> A(m * 3), rhs(m);
        for (int i = 0; i < m; i++) {
            double xx = x[a + i] - x[k];
            double w = y[a + i] / top;
            A[i * 3 + 0] = xx * xx * w;
            A[i * 3 + 1] = xx * w;
            A[i * 3 + 2] = w;
            rhs[i] = std::log(std::max(y[a + i], 1e-300)) * w;
        }
        // column scaling as polyfit does (conditioning only)
        double sc[3];
        for (int j = 0; j < 3; j++) {
            double s = 0;
            for (int i = 0; i < m; i++) s += A[i * 3 + j] * A[i * 3 + j];
            sc[j] = std::sqrt(s);
            if (sc[j] > 0) for (int i = 0; i < m; i++) A[i * 3 + j] /= sc[j];
        }
        std::vector<double> c;
        lstsq_small(A, rhs, m, 3, c);
        for (int j = 0; j < 3; j++) if (sc[j] > 0) c[j] /= sc[j];
        double c2 = c[0], c1 = c[1];
        if (c2 < 0 && std::fabs(c1 / (2 * c2)) <= (x[b] - x[a])) return x[k] - c1 / (2 * c2);
    }
    if (y[k - 1] > 0 && y[k + 1] > 0) {
        double la = std::log(y[k - 1]), lb = std::log(top), lc = std::log(y[k + 1]);
        double den = la - 2 * lb + lc;
        if (den < 0) {
            double off = std::max(-0.5, std::min(0.5, 0.5 * (la - lc) / den));
            double step = off >= 0 ? x[k + 1] - x[k] : x[k] - x[k - 1];
            return x[k] + off * step;
        }
    }
    return x[k];
}

double fwhm_at(const double* x, const double* y, size_t n, size_t i) {
    double half = y[i] / 2.0;
    size_t lo = i, hi = i;
    while (lo > 0 && y[lo] > half) lo--;
    while (hi < n - 1 && y[hi] > half) hi++;
    if (y[lo] > half || y[hi] > half || lo == i || hi == i) return -1.0;
    double xl = x[lo] + (half - y[lo]) * (x[lo + 1] - x[lo]) / std::max(y[lo + 1] - y[lo], 1e-30);
    double xr = x[hi - 1] + (y[hi - 1] - half) * (x[hi] - x[hi - 1]) / std::max(y[hi - 1] - y[hi], 1e-30);
    return xr - xl;
}

double estimate_peak_width(const double* x, const double* y, size_t n) {
    if (n == 0) return 1.0;
    size_t i = std::max_element(y, y + n) - y;
    double w = n > 1 ? fwhm_at(x, y, n, i) : -1.0;
    if (!(w > 0) || !std::isfinite(w)) w = 5.0 * median_diff(x, n);
    return w;
}

double estimate_resolution(const double* x, const double* y, size_t n, int npk) {
    if (n < 5) return 10000.0;
    double top = *std::max_element(y, y + n);
    if (top <= 0) return 10000.0;
    std::vector<size_t> order = order_desc(y, n);
    std::vector<double> used, vals;
    size_t lim = std::min<size_t>(5000, order.size());
    for (size_t t = 0; t < lim; t++) {
        size_t i = order[t];
        if (y[i] < 0.05 * top || (int)vals.size() >= npk) break;
        bool near = false;
        for (double u : used) if (std::fabs(x[i] - u) < 0.2) { near = true; break; }
        if (near) continue;
        used.push_back(x[i]);
        double w = fwhm_at(x, y, n, i);
        if (w > 0) vals.push_back(x[i] / w);
    }
    if (!vals.empty()) return median(vals);
    return x[order[0]] / std::max(estimate_peak_width(x, y, n), 1e-9);
}

}  // namespace e2
}  // namespace ms

// ------------------------------------------------------------ C interface
namespace {
struct PeakOut {
    std::vector<double> mass, apex, height, area;
    std::vector<int> n_iso;
    void set(const std::vector<ms::e2::Peak>& p) {
        mass.clear(); apex.clear(); height.clear(); area.clear(); n_iso.clear();
        for (auto& q : p) { mass.push_back(q.mass); apex.push_back(q.apex); height.push_back(q.height); area.push_back(q.area); n_iso.push_back(q.n_iso); }
    }
};
thread_local PeakOut g_pick, g_group;
}

extern "C" {

MS_API int ms_pick_peaks(const double* mass, const double* it, long n, double window, double threshold,
                         const double** pmass, const double** height, const double** area, long* np) {
    return ms::guarded([&] {
        if (!mass || !it || n < 0 || std::isnan(window) || std::isnan(threshold)) throw std::invalid_argument("ms_pick_peaks: bad arguments");
        g_pick.set(ms::e2::pick_mass_peaks(mass, it, (size_t)n, window, threshold));
        if (pmass) *pmass = g_pick.mass.data();
        if (height) *height = g_pick.height.data();
        if (area) *area = g_pick.area.data();
        if (np) *np = (long)g_pick.mass.size();
    });
}

MS_API int ms_group_isotopes(const double* mass, const double* it, long n, double threshold,
                             const double** avg, const double** apex, const double** height,
                             const double** area, const int** n_iso, long* np) {
    return ms::guarded([&] {
        if (!mass || !it || n < 0 || std::isnan(threshold)) throw std::invalid_argument("ms_group_isotopes: bad arguments");
        g_group.set(ms::e2::group_isotopes(mass, it, (size_t)n, threshold));
        if (avg) *avg = g_group.mass.data();
        if (apex) *apex = g_group.apex.data();
        if (height) *height = g_group.height.data();
        if (area) *area = g_group.area.data();
        if (n_iso) *n_iso = g_group.n_iso.data();
        if (np) *np = (long)g_group.mass.size();
    });
}

}

// ---- 3.1 additions (agent P): the label, species and data helpers of ms_deconv ---------
namespace ms {
namespace p31 {
using ms::e2::np_sum;

// np.argmax of y[lo:hi] (first maximum; a NaN wins)
size_t argmax(const double* y, size_t lo, size_t hi) {
    size_t k = lo;
    for (size_t i = lo + 1; i < hi && !std::isnan(y[k]); i++)
        if (std::isnan(y[i]) || y[i] > y[k]) k = i;
    return k;
}

// ms_deconv._find_peaks: local maxima at least `height` high and `distance` points apart
std::vector<size_t> find_peaks_np(const double* y, size_t n, double height, long distance) {
    std::vector<size_t> cand;
    if (n < 3) return cand;
    for (size_t i = 1; i + 1 < n; i++)
        if (y[i] > y[i - 1] && y[i] >= y[i + 1] && y[i] >= height) cand.push_back(i);
    if (cand.size() > 1) {   // a flat top counts once (compared with the previous candidate)
        std::vector<size_t> kept{cand[0]};
        for (size_t k = 1; k < cand.size(); k++)
            if (!(cand[k] - cand[k - 1] == 1 && y[cand[k]] == y[cand[k - 1]])) kept.push_back(cand[k]);
        cand.swap(kept);
    }
    if (distance <= 1 || cand.size() < 2) return cand;
    std::vector<size_t> order(cand);
    std::stable_sort(order.begin(), order.end(), [y](size_t a, size_t b) { return -y[a] < -y[b]; });
    std::vector<char> taken(n, 0);
    std::vector<size_t> out;
    for (size_t i : order) {
        if (taken[i]) continue;
        out.push_back(i);
        const long a = std::max(0L, (long)i - distance + 1);
        const long b = std::min((long)n, (long)i + distance);
        for (long j = a; j < b; j++) taken[j] = 1;
    }
    std::sort(out.begin(), out.end());
    return out;
}

// ms_deconv.isotope_peaks
std::vector<std::pair<double, double>> isotope_peaks(const double* x, const double* y, size_t n, double threshold, long max_n) {
    std::vector<std::pair<double, double>> out;
    if (n < 5) return out;
    double top = y[0];
    for (size_t i = 1; i < n; i++) if (y[i] > top) top = y[i];
    if (top <= 0) return out;
    const double dx = ms::e2::median_diff(x, n);
    const double dd = 0.4 / std::max(dx, 1e-9);
    const long dist = std::max(1L, dd < 4e18 ? (long)dd : (long)4e18);
    const std::vector<size_t> idx = find_peaks_np(y, n, threshold * top, dist);
    std::vector<double> xy;
    for (size_t i : idx) {
        size_t lo = i, hi = i;
        while (lo > 0 && y[lo - 1] >= 0.5 * y[i] && y[lo - 1] <= y[lo]) lo--;
        while (hi < n - 1 && y[hi + 1] >= 0.5 * y[i] && y[hi + 1] <= y[hi]) hi++;
        double cen = x[i];
        if (hi > lo) {
            xy.resize(hi - lo + 1);
            for (size_t j = lo; j <= hi; j++) xy[j - lo] = x[j] * y[j];
            cen = np_sum(xy.data(), xy.size()) / np_sum(y + lo, hi - lo + 1);
        }
        out.push_back({cen, y[i]});
    }
    std::stable_sort(out.begin(), out.end(), [](const std::pair<double, double>& a, const std::pair<double, double>& b) {
        return -a.second < -b.second; });
    if (max_n >= 0 && (long)out.size() > max_n) out.resize(max_n);
    std::sort(out.begin(), out.end());
    return out;
}

// _raw_isotope_areas: area of each isotope peak in the m/z data, summed over the charges
std::vector<double> raw_isotope_areas(const double* x, const double* y, size_t n, const std::vector<double>& masses,
                                      double a, int z0, int z1) {
    std::vector<double> out(masses.size(), 0.0);
    std::vector<double> prod;
    for (int z = std::max(1, z0); z <= std::max(z0, z1); z++) {
        const double hw = 0.45 / z;
        for (size_t i = 0; i < masses.size(); i++) {
            const double c = (masses[i] + z * a) / z;
            const size_t lo = std::lower_bound(x, x + n, c - hw) - x;
            const size_t hi = std::lower_bound(x, x + n, c + hw) - x;
            if (hi < lo + 2) continue;
            double floor = y[lo];
            for (size_t j = lo + 1; j < hi; j++) if (y[j] < floor || std::isnan(y[j])) { floor = y[j]; if (std::isnan(floor)) break; }
            prod.resize(hi - lo - 1);
            for (size_t j = lo + 1; j < hi; j++) prod[j - lo - 1] = 0.5 * (y[j] + y[j - 1] - 2 * floor) * (x[j] - x[j - 1]);
            out[i] += np_sum(prod.data(), prod.size());
        }
    }
    return out;
}

// refine_species for one species: out = status (1 checked, 0 left alone), apex, apex2 (NaN: none),
// apex2_rel, source (0 deconvolution, 1 m/z data), best (-1, 0 or 1: the isotope taken as the apex)
void refine_species(const double* x, const double* y, size_t n, const double* dx, const double* dy, size_t nd,
                    double carrier, int z0, int z1, double spacing, double ap, double* out) {
    out[0] = 0; out[1] = ap; out[2] = std::numeric_limits<double>::quiet_NaN(); out[3] = 0; out[4] = 0; out[5] = 0;
    std::map<int, double> pos, hdec;
    for (int k = -3; k <= 3; k++) {
        const double c = ap + k * spacing;
        const size_t lo = std::lower_bound(x, x + n, c - 0.3) - x;
        const size_t hi = std::lower_bound(x, x + n, c + 0.3) - x;
        if (hi > lo) {
            const size_t t = argmax(y, lo, hi);
            hdec[k] = y[t];
            pos[k] = k == 0 ? ap : ms::e2::refine_peak_mass(x, y, n, x[t]);
        }
    }
    if (!(hdec.count(0) && hdec[0] > 0)) return;
    std::vector<int> ks;
    for (int k = -3; k <= 3; k++) if (pos.count(k)) ks.push_back(k);
    std::map<int, double> heights = hdec;
    int source = 0;
    if (nd > 10) {
        std::vector<double> ms_;
        for (int k : ks) ms_.push_back(pos[k]);
        const std::vector<double> areas = raw_isotope_areas(dx, dy, nd, ms_, carrier, z0, z1);
        double amax = areas[0];
        for (double v : areas) if (std::isnan(amax) || std::isnan(v) || v > amax) { amax = v; if (std::isnan(v)) break; }
        const size_t i0 = std::find(ks.begin(), ks.end(), 0) - ks.begin();
        if (amax > 0 && areas[i0] > 0) {
            heights.clear();
            for (size_t i = 0; i < ks.size(); i++) heights[ks[i]] = areas[i];
            source = 1;
        }
    }
    int best = 0;
    bool have = false;
    for (int k : {-1, 0, 1}) {
        if (!heights.count(k)) continue;
        if (!have || heights[k] > heights[best]) { best = k; have = true; }
    }
    out[0] = 1;
    out[4] = source;
    out[5] = best;
    out[1] = best != 0 ? pos[best] : ap;
    const double hb = heights[best];
    int k2 = 0;
    bool c2 = false;
    for (int k : {best - 1, best + 1}) {
        if (!heights.count(k) || !(heights[k] >= 0.95 * hb)) continue;
        if (!c2 || heights[k] > heights[k2]) { k2 = k; c2 = true; }
    }
    if (c2) {
        out[2] = pos[k2];
        out[3] = 100.0 * heights[k2] / hb;
    }
}

// np.minimum / np.maximum (a NaN wins)
inline double np_min(double a, double b) { return std::isnan(a) || std::isnan(b) ? std::numeric_limits<double>::quiet_NaN() : (b < a ? b : a); }

// scipy.ndimage minimum_filter1d / maximum_filter1d (size k, mode 'nearest'): scipy's own algorithm
// (NI_MinOrMaxFilter1D, R. Harter's ascending minima in a ring of k pairs), so that NaN behaves as
// there too: a NaN never wins a comparison and only comes out while it is the oldest entry left, so a
// short run of NaN stays local (the van Herk blocks used before carried it over whole blocks, and the
// running sum of the smoothing then over the rest of the trace). Exact (no rounding): finite data give
// the same values as before.
template <bool IS_MAX>
void minmax_filter_t(const std::vector<double>& in, std::vector<double>& out, long k) {
    const long n = (long)in.size();
    out.assign(n, 0.0);
    if (n == 0) return;
    if (k <= 1) { out = in; return; }
    const long s1 = k / 2;
    const long m = n + k - 1;   // the line extended by its end values ('nearest'): ext[j] = in[clip(j - s1)]
    auto ext = [&](long j) { return in[std::max(0L, std::min(n - 1, j - s1))]; };
    struct Pair { double value; long death; };
    std::vector<Pair> ring(k);
    Pair* const begin = ring.data();
    Pair* const end = begin + k;
    auto inc = [&](Pair*& p) { if (++p >= end) p = begin; };
    auto dec = [&](Pair*& p) { if (p == begin) p = end; --p; };
    auto better = [](double a, double b) { return IS_MAX ? a >= b : a <= b; };   // a replaces b (false with NaN)
    Pair* best = begin;
    best->value = ext(0);
    best->death = k;
    Pair* last = begin;
    long o = 0;
    for (long l = 1; l < m; l++) {
        const double v = ext(l);
        if (best->death == l) inc(best);
        if (better(v, best->value)) {
            best->value = v;
            best->death = l + k;
            last = best;
        } else {
            while (IS_MAX ? last->value <= v : last->value >= v) dec(last);   // (stops at best at the latest)
            inc(last);
            last->value = v;
            last->death = l + k;
        }
        if (l >= k - 1) out[o++] = best->value;
    }
}

void minmax_filter(const std::vector<double>& in, std::vector<double>& out, long k, bool is_max) {
    if (is_max) minmax_filter_t<true>(in, out, k);
    else minmax_filter_t<false>(in, out, k);
}

// scipy.ndimage.uniform_filter1d (mode 'nearest'): running sum, divided at each point
void uniform_filter(const std::vector<double>& in, std::vector<double>& out, long k) {
    const long n = (long)in.size();
    out.assign(n, 0.0);
    const long s1 = k / 2, s2 = k - s1 - 1;
    auto at = [&](long i) { return in[std::max(0L, std::min(n - 1, i))]; };
    double tmp = 0.0;
    for (long j = 0; j < k; j++) tmp += at(j - s1);
    out[0] = tmp / k;
    for (long i = 1; i < n; i++) {
        tmp += at(i + s2) - at(i - 1 - s1);
        out[i] = tmp / k;
    }
}

// ms_deconv.subtract_baseline (the intensities)
void subtract_baseline(const double* x, const double* y, size_t n, double width, double factor, double* out) {
    if (n < 20) { std::copy(y, y + n, out); return; }
    if (!(width != 0) || std::isnan(width)) width = factor * ms::e2::estimate_peak_width(x, y, n);
    const double dx = ms::e2::median_diff(x, n);
    const double r = ms::e2::py_round(width / std::max(dx, 1e-12));
    const double kk = std::max(5.0, std::min((double)(n / 4), r));
    const long k = (long)kk;
    const std::vector<double> yv(y, y + n);
    std::vector<double> b1, b2, b3;
    minmax_filter(yv, b1, k, false);
    minmax_filter(b1, b2, k, true);
    uniform_filter(b2, b3, std::max(3L, k / 2));
    for (size_t i = 0; i < n; i++) {
        const double v = y[i] - np_min(b3[i], y[i]);
        out[i] = std::isnan(v) ? v : (v < 0 ? 0.0 : v);   // np.clip(.., 0, None)
    }
}

// ms_deconv.peak_mask
long peak_mask(const double* x, const double* y, long n, double thr, unsigned char* keep) {
    std::fill(keep, keep + n, 0);
    if (n < 3) {
        long c = 0;
        for (long i = 0; i < n; i++) if (y[i] >= thr) { keep[i] = 1; c++; }
        return c;
    }
    // find_peaks(y, height=thr, plateau_size=(1, None)): every local maximum (the middle of a flat top)
    std::vector<long> idx;
    long i = 1;
    const long imax = n - 1;
    while (i < imax) {
        if (y[i - 1] < y[i]) {
            long ahead = i + 1;
            while (ahead < imax && y[ahead] == y[i]) ahead++;
            if (y[ahead] < y[i]) {
                const long mid = (i + ahead - 1) / 2;
                if (y[mid] >= thr) idx.push_back(mid);
                i = ahead;
            }
        }
        i++;
    }
    if (idx.empty()) return 0;
    const double fw = ms::e2::estimate_peak_width(x, y, (size_t)n);
    const long i0 = (long)argmax(y, 0, (size_t)n);
    const long lo = std::max(0L, i0 - 50), hi = std::min(n - 1, i0 + 50);
    const double dx = hi > lo ? ms::e2::median_diff(x + lo, (size_t)(hi - lo + 1)) : ms::e2::median_diff(x, (size_t)n);
    const double n_fw = dx > 0 ? fw / dx : 5.0;
    const double wv = std::max(7.0, std::min(10 * n_fw, (double)n));
    const long wlen = ((long)wv) | 1;
    // peak_prominences(y, idx, wlen): the bases of each peak
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
    return (long)idx.size();
}

// unilcms.label_maxima: the nlab highest local maxima at least rel x the tallest one, further
// apart than sep (NaN: a fortieth of the m/z span); ties: the later point first (argsort()[::-1])
std::vector<long> label_maxima(const double* x, const double* y, long n, int nlab, double sep, double rel) {
    std::vector<long> chosen;
    if (n < 3) return chosen;
    std::vector<long> lm;
    for (long i = 1; i + 1 < n; i++) if (y[i] >= y[i - 1] && y[i] > y[i + 1]) lm.push_back(i);
    if (lm.empty()) return chosen;
    // the selection stops at the first maximum below rel x the tallest: only those above are sorted
    double top = y[lm[0]];
    for (long i : lm) if (y[i] > top) top = y[i];
    std::vector<long> order;
    for (auto it = lm.rbegin(); it != lm.rend(); ++it) if (y[*it] >= rel * top) order.push_back(*it);
    std::stable_sort(order.begin(), order.end(), [y](long a, long b) { return y[a] > y[b]; });
    if (order.empty()) return chosen;
    if (std::isnan(sep)) sep = (x[n - 1] - x[0]) / 40.0;
    for (long i : order) {
        if (y[i] < rel * top || (long)chosen.size() >= nlab) break;
        bool ok = true;
        for (long j : chosen) if (!(std::fabs(x[i] - x[j]) > sep)) { ok = false; break; }
        if (ok) chosen.push_back(i);
    }
    return chosen;
}

struct Out31 {
    std::vector<double> a, b, c, d, e, f;
    std::vector<int> ni;
};
thread_local Out31 g_isop, g_group2, g_restrict;

}  // namespace p31
}  // namespace ms

extern "C" {

MS_API int ms_isotope_peaks(const double* mass, const double* it, long n, double threshold, long max_n,
                            const double** apex, const double** height, long* np) {
    return ms::guarded([&] {
        if (!apex || !height || !np || n < 0 || (n > 0 && (!mass || !it))) throw std::invalid_argument("ms_isotope_peaks: bad arguments");
        const auto r = ms::p31::isotope_peaks(mass, it, (size_t)n, threshold, max_n);
        ms::p31::Out31& o = ms::p31::g_isop;
        o.a.clear(); o.b.clear();
        for (const auto& q : r) { o.a.push_back(q.first); o.b.push_back(q.second); }
        *apex = o.a.data();
        *height = o.b.data();
        *np = (long)r.size();
    });
}

MS_API double ms_refine_peak_mass(const double* mass, const double* it, long n, double m0) {
    return ms::guarded_value(m0, [&] {
        if (!mass || !it || n < 0) throw std::invalid_argument("ms_refine_peak_mass: bad arguments");
        return ms::e2::refine_peak_mass(mass, it, (size_t)n, m0);
    });
}

MS_API int ms_refine_species(const double* mass, const double* it, long n, const double* mz, const double* mz_it, long nd,
                             double carrier, int z_lo, int z_hi, double spacing, double apex, double* out) {
    return ms::guarded([&] {
        if (!mass || !it || !out || n < 0 || nd < 0 || (nd > 0 && (!mz || !mz_it)))
            throw std::invalid_argument("ms_refine_species: bad arguments");
        if (std::isnan(apex) || std::isnan(spacing)) throw std::invalid_argument("ms_refine_species: NaN apex or spacing");
        ms::p31::refine_species(mass, it, (size_t)n, mz, mz_it, (size_t)nd, carrier, z_lo, z_hi, spacing, apex, out);
    });
}

MS_API int ms_group_isotopes2(const double* mass, const double* it, long n, double threshold, double spacing,
                              const double** avg, const double** apex, const double** height, const double** area,
                              const int** n_iso, const double** first, const double** last, long* np) {
    return ms::guarded([&] {
        if (!avg || !apex || !height || !area || !n_iso || !first || !last || !np || n < 0 || (n > 0 && (!mass || !it)))
            throw std::invalid_argument("ms_group_isotopes2: bad arguments");
        if (std::isnan(threshold) || std::isnan(spacing)) throw std::invalid_argument("ms_group_isotopes2: NaN argument");
        std::vector<std::pair<double, double>> ranges;
        const std::vector<ms::e2::Peak> r = ms::e2::group_isotopes_ranges(mass, it, (size_t)n, threshold, &ranges, spacing);
        ms::p31::Out31& o = ms::p31::g_group2;
        o.a.clear(); o.b.clear(); o.c.clear(); o.d.clear(); o.e.clear(); o.f.clear(); o.ni.clear();
        for (size_t k = 0; k < r.size(); k++) {
            o.a.push_back(r[k].mass); o.b.push_back(r[k].apex); o.c.push_back(r[k].height); o.d.push_back(r[k].area);
            o.ni.push_back(r[k].n_iso); o.e.push_back(ranges[k].first); o.f.push_back(ranges[k].second);
        }
        *avg = o.a.data(); *apex = o.b.data(); *height = o.c.data(); *area = o.d.data(); *n_iso = o.ni.data();
        *first = o.e.data(); *last = o.f.data();
        *np = (long)r.size();
    });
}

// ms_deconv._restrict: finite rows with m/z > 0 (negative intensities: MS_BAD_ARG), sorted by m/z
// (stable), within lo .. hi when hi > lo
MS_API int ms_restrict(const double* mz, const double* it, long n, double lo, double hi, const double** omz,
                       const double** oit, long* nout) {
    return ms::guarded([&] {
        if (!omz || !oit || !nout || n < 0 || (n > 0 && (!mz || !it))) throw std::invalid_argument("ms_restrict: bad arguments");
        std::vector<long> idx;
        idx.reserve(n);
        bool sorted = true;
        for (long i = 0; i < n; i++) {
            if (!std::isfinite(mz[i]) || !std::isfinite(it[i]) || !(mz[i] > 0)) continue;
            if (it[i] < 0) throw std::invalid_argument("Spectrum intensities must not be negative");
            if (!idx.empty() && mz[i] < mz[idx.back()]) sorted = false;
            idx.push_back(i);
        }
        if (!sorted) std::stable_sort(idx.begin(), idx.end(), [&](long a, long b) { return mz[a] < mz[b]; });
        const bool range = hi > lo;
        ms::p31::Out31& o = ms::p31::g_restrict;
        o.a.clear(); o.b.clear();
        for (long i : idx) {
            if (range && !(mz[i] >= lo && mz[i] <= hi)) continue;
            o.a.push_back(mz[i]);
            o.b.push_back(it[i]);
        }
        *omz = o.a.data();
        *oit = o.b.data();
        *nout = (long)o.a.size();
    });
}

MS_API int ms_peak_mask(const double* mz, const double* it, long n, double thr, unsigned char* mask, long* npeaks) {
    return ms::guarded([&] {
        if (!mask || !npeaks || n < 0 || (n > 0 && (!mz || !it))) throw std::invalid_argument("ms_peak_mask: bad arguments");
        if (std::isnan(thr)) throw std::invalid_argument("ms_peak_mask: NaN threshold");
        *npeaks = ms::p31::peak_mask(mz, it, n, thr, mask);
    });
}

MS_API int ms_subtract_baseline(const double* mz, const double* it, long n, double width, double factor, double* out) {
    return ms::guarded([&] {
        if (!out || n < 0 || (n > 0 && (!mz || !it))) throw std::invalid_argument("ms_subtract_baseline: bad arguments");
        if (!std::isfinite(factor)) throw std::invalid_argument("ms_subtract_baseline: the factor must be finite");
        ms::p31::subtract_baseline(mz, it, (size_t)n, width, factor, out);
    });
}

MS_API int ms_label_maxima(const double* x, const double* y, long n, int nlab, double sep, double rel, long* idx, long* nout) {
    return ms::guarded([&] {
        if (!idx || !nout || n < 0 || nlab < 0 || (n > 0 && (!x || !y))) throw std::invalid_argument("ms_label_maxima: bad arguments");
        const std::vector<long> r = ms::p31::label_maxima(x, y, n, nlab, sep, rel);
        std::copy(r.begin(), r.end(), idx);
        *nout = (long)r.size();
    });
}

}

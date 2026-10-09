// Internal m/z calibration: port of hrms_calib.py (Calibration, _hpc_fit,
// _peak_position, find_calibrants).
//
// Models (C numbering of msengine.h; the Python class numbers them 0..5 with
// 0 = single point, so C model = Python model - 1):
//   0 linear      m' = m + c0 + c1 x                 x = (m - x0) / s
//   1 quadratic   m' = m + c0 + c1 x + c2 x^2
//   2 cubic       m' = m + c0 + c1 x + c2 x^2 + c3 x^3
//   3 TOF         sqrt(m') = a + b sqrt(m) + c m
//   4 HPC         m1 = quadratic as above; m' = m1 (1 + sum_j h_j P_j(t)),
//                 P_j Legendre polynomials, t = 2 (m - lo) / (hi - lo) - 1 clipped to -1..1
//   5 single point (ppm shift, Python model 0; an addition to the header list)  m' = k m
// with x0 = mean of the measured calibrant m/z, s = max(max - min, 1), lo / hi =
// min / max of the measured calibrant m/z (the m_lo / m_hi of ms_calibration_apply).
//
// Coefficient array written by ms_calibration_fit and read by ms_calibration_apply:
//   coef[0] = x0, coef[1] = s, then the model coefficients:
//   model 0: c0 c1            (ncoef 4)      model 3: a b c            (ncoef 5)
//   model 1: c0 c1 c2         (ncoef 5)      model 4: c0 c1 c2 h0 .. hp (ncoef 6 + p)
//   model 2: c0 c1 c2 c3      (ncoef 6)      model 5: k                (ncoef 3)
// For HPC the correction order p follows from ncoef (p = ncoef - 6); the `order`
// argument of ms_calibration_fit chooses it (<= 0: by leave-one-out cross
// validation as the Python code does, > 0: that order, at most min(8, n - 8)).
// The least squares problems are solved with a Householder QR after scaling the
// columns (numpy uses an SVD; the solutions agree to rounding).
#include "common.h"
#include "e2_internal.h"   // np_sum: numpy's pairwise summation (the means and sums of hrms_calib)
#include <algorithm>
#include <limits>
#include <numeric>

namespace ms {
namespace {

const double NaN = std::numeric_limits<double>::quiet_NaN();

// least squares min |A c - b|, A column major (m rows, k columns); columns are
// scaled by their norms first (as np.polyfit / legfit do). false when singular.
bool lstsq(std::vector<double> A, std::vector<double> b, long m, int k, std::vector<double>& c) {
    if (m < k || k <= 0) return false;
    std::vector<double> scale(k);
    for (int j = 0; j < k; j++) {
        double s = 0;
        for (long i = 0; i < m; i++) s += A[j * m + i] * A[j * m + i];
        scale[j] = std::sqrt(s);
        if (!(scale[j] > 0)) scale[j] = 1.0;
        for (long i = 0; i < m; i++) A[j * m + i] /= scale[j];
    }
    std::vector<double> R(k * k, 0.0), v(m);
    for (int j = 0; j < k; j++) {
        double norm = 0;
        for (long i = j; i < m; i++) norm += A[j * m + i] * A[j * m + i];
        norm = std::sqrt(norm);
        if (!(norm > 0)) return false;
        const double alpha = A[j * m + j] > 0 ? -norm : norm;
        std::fill(v.begin(), v.end(), 0.0);
        for (long i = j; i < m; i++) v[i] = A[j * m + i];
        v[j] -= alpha;
        double vv = 0;
        for (long i = j; i < m; i++) vv += v[i] * v[i];
        if (vv > 0) {
            for (int q = j; q < k; q++) {
                double d = 0;
                for (long i = j; i < m; i++) d += v[i] * A[q * m + i];
                d = 2 * d / vv;
                for (long i = j; i < m; i++) A[q * m + i] -= d * v[i];
            }
            double d = 0;
            for (long i = j; i < m; i++) d += v[i] * b[i];
            d = 2 * d / vv;
            for (long i = j; i < m; i++) b[i] -= d * v[i];
        }
        for (int q = j; q < k; q++) R[j * k + q] = A[q * m + j];
    }
    c.assign(k, 0.0);
    for (int j = k - 1; j >= 0; j--) {
        double s = b[j];
        for (int q = j + 1; q < k; q++) s -= R[j * k + q] * c[q];
        if (!(std::fabs(R[j * k + j]) > 0)) return false;
        c[j] = s / R[j * k + j];
    }
    for (int j = 0; j < k; j++) {
        c[j] /= scale[j];
        if (!std::isfinite(c[j])) return false;
    }
    return true;
}

// numpy.polynomial.legendre.legval (Clenshaw recurrence, same operation order)
double legval(double x, const std::vector<double>& c) {
    const int n = (int)c.size();
    double c0, c1;
    if (n == 1) { c0 = c[0]; c1 = 0; }
    else if (n == 2) { c0 = c[0]; c1 = c[1]; }
    else {
        int nd = n;
        c0 = c[n - 2];
        c1 = c[n - 1];
        for (int i = 3; i <= n; i++) {
            const double tmp = c0;
            nd = nd - 1;
            c0 = c[n - i] - (c1 * (nd - 1)) / nd;
            c1 = tmp + (c1 * x * (2 * nd - 1)) / nd;
        }
    }
    return c0 + c1 * x;
}

// legendre.legvander column j at x
void legvander(double x, int deg, std::vector<double>& v) {
    v.assign(deg + 1, 1.0);
    if (deg > 0) {
        v[1] = x;
        for (int i = 2; i <= deg; i++) v[i] = (v[i - 1] * x * (2 * i - 1) - v[i - 2] * (i - 1)) / i;
    }
}

// x ** k of numpy for a whole number k: 1, x and x * x for k = 0, 1, 2 (its fast paths ones_like,
// positive and square), pow() above
inline double xpow(double x, int k) {
    switch (k) {
        case 0: return 1.0;
        case 1: return x;
        case 2: return x * x;
    }
    return std::pow(x, (double)k);
}

double hpc_t(double mz, double lo, double hi) {
    double t = 2.0 * (mz - lo) / std::max(hi - lo, 1e-9) - 1.0;
    return std::max(-1.0, std::min(1.0, t));
}

struct Cal {
    int model = 0;              // C numbering
    double x0 = 0, s = 1, lo = 0, hi = 0;
    std::vector<double> coef;   // model coefficients (without x0, s)
    std::vector<double> hpc;    // HPC correction
    int order = 0;              // HPC order

    double operator()(double mz) const {
        if (model == 5) return mz * coef[0];
        if (model == 3) {
            const double v = coef[0] + coef[1] * std::sqrt(std::max(mz, 0.0)) + coef[2] * mz;
            return v * v;
        }
        const double x = (mz - x0) / s;
        double corr = 0.0;
        for (size_t k = 0; k < coef.size(); k++) corr = corr + coef[k] * xpow(x, (int)k);
        const double m1 = mz + corr;
        if (model == 4) return m1 * (1.0 + legval(hpc_t(mz, lo, hi), hpc));
        return m1;
    }
};

int n_params(int model) {
    switch (model) {
        case 0: return 2; case 1: return 3; case 2: return 4; case 3: return 3; case 4: return 5; case 5: return 1;
    }
    throw std::invalid_argument("unknown calibration model");
}

const char* model_name(int model) {
    static const char* names[] = {"Linear", "Quadratic", "Cubic", "TOF (square root)", "HPC (quadratic + correction)",
                                  "Single point (ppm shift)"};
    return names[model];
}

Cal fit_plain(int model, const std::vector<double>& m, const std::vector<double>& r);

// _hpc_coef: Legendre fit of the relative error left by the quadratic
std::vector<double> hpc_coef(const std::vector<double>& m, const std::vector<double>& r, const Cal& base,
                             double lo, double hi, int p) {
    const long n = (long)m.size();
    std::vector<double> A((p + 1) * n), b(n), v, c;
    for (long i = 0; i < n; i++) {
        const double bm = base(m[i]);
        b[i] = (r[i] - bm) / bm;
        legvander(hpc_t(m[i], lo, hi), p, v);
        for (int j = 0; j <= p; j++) A[j * n + i] = v[j];
    }
    if (!lstsq(A, b, n, p + 1, c)) throw Error("HPC correction: singular fit");
    return c;
}

// hrms_calib.hpc_max_order: the quadratic and the correction leave at least 4
// calibrants of freedom (with n - 4 the fit went through every calibrant and the
// leave-one-out choice followed the noise on few calibrants)
int hpc_max_order(long n) { return (int)std::max<long>(0, std::min<long>(8, n - 8)); }

Cal fit_hpc(const std::vector<double>& m, const std::vector<double>& r, int order) {
    const long n = (long)m.size();
    Cal c = fit_plain(1, m, r);   // the quadratic base (model 1 while the correction is fitted)
    const int top = hpc_max_order(n);
    int best_p;
    if (order >= 0) {
        best_p = std::max(0, std::min(order, top));
    } else {
        best_p = 0;
        double best_e = NaN;
        for (int p = 0; p <= std::max(top, 0); p++) {
            std::vector<double> e2v(n);
            for (long i = 0; i < n; i++) {
                std::vector<double> mk, rk;
                for (long j = 0; j < n; j++) if (j != i) { mk.push_back(m[j]); rk.push_back(r[j]); }
                const Cal b = fit_plain(1, mk, rk);
                const double lo = *std::min_element(mk.begin(), mk.end()), hi = *std::max_element(mk.begin(), mk.end());
                const std::vector<double> cc = hpc_coef(mk, rk, b, lo, hi, p);
                const double pred = b(m[i]) * (1.0 + legval(hpc_t(m[i], lo, hi), cc));
                const double e = (pred - r[i]) / r[i];
                e2v[i] = e * e;
            }
            // np.sqrt(np.mean(np.square(errs))): pairwise sum as numpy
            const double e = std::sqrt(ms::e2::np_sum(e2v.data(), (size_t)n) / n);
            if (std::isnan(best_e) || e < 0.95 * best_e) { best_p = p; best_e = e; }
        }
    }
    c.order = best_p;
    c.hpc = hpc_coef(m, r, c, c.lo, c.hi, best_p);
    c.model = 4;
    return c;
}

Cal fit_plain(int model, const std::vector<double>& m, const std::vector<double>& r) {
    const long n = (long)m.size();
    const int need = n_params(model);
    if (n < need) {
        char buf[160];
        snprintf(buf, sizeof buf, "%s needs at least %d calibrant peak%s (found %ld)", model_name(model), need,
                 need > 1 ? "s" : "", n);
        throw Error(buf);
    }
    Cal c;
    c.model = model;
    c.lo = *std::min_element(m.begin(), m.end());
    c.hi = *std::max_element(m.begin(), m.end());
    c.x0 = ms::e2::np_sum(m.data(), (size_t)n) / n;   // m.mean() (numpy's pairwise sum)
    c.s = std::max(c.hi - c.lo, 1.0);
    if (model == 5) {
        std::vector<double> q(n), qq(n);
        for (long i = 0; i < n; i++) { q[i] = m[i] / r[i]; qq[i] = q[i] * q[i]; }
        c.coef = {ms::e2::np_sum(q.data(), (size_t)n) / ms::e2::np_sum(qq.data(), (size_t)n)};
    } else if (model == 3) {
        std::vector<double> A(3 * n), b(n, 1.0);
        for (long i = 0; i < n; i++) {
            const double sr = std::sqrt(r[i]);
            A[0 * n + i] = 1.0 / sr;
            A[1 * n + i] = std::sqrt(m[i]) / sr;
            A[2 * n + i] = m[i] / sr;
        }
        if (!lstsq(A, b, n, 3, c.coef)) throw Error("TOF calibration: singular fit");
    } else if (model == 4) {
        return fit_hpc(m, r, -1);
    } else {
        const int k = model + 2;   // deg + 1
        std::vector<double> A(k * n), b(n);
        for (long i = 0; i < n; i++) {
            const double x = (m[i] - c.x0) / c.s;
            for (int j = 0; j < k; j++) A[j * n + i] = xpow(x, j) / r[i];
            b[i] = (r[i] - m[i]) / r[i];
        }
        if (!lstsq(A, b, n, k, c.coef)) throw Error("calibration: singular fit");
    }
    return c;
}

// order: HPC correction order, < 0 = chosen by cross validation
Cal fit_model(int model, const std::vector<double>& m, const std::vector<double>& r, int order) {
    if (model == 4) {
        if ((long)m.size() < n_params(4)) fit_plain(4, m, r);   // throws the message
        return fit_hpc(m, r, order);
    }
    return fit_plain(model, m, r);
}

// a column of values every `s` doubles (the columns of an (n, 2) array without a copy: s = 2)
struct Col {
    const double* p;
    long s;
    double operator[](long i) const { return p[i * s]; }
};

template <class Arr> bool apex_at(const Arr& mz, const Arr& it, long n, long k, bool profile, double& pos, double& top);

// _peak_position: the most intense point in [lo, hi] (profile: Gaussian apex)
template <class Arr>
bool peak_position(const Arr& mz, const Arr& it, long n, double lo, double hi, bool profile, double& pos, double& top) {
    long k = -1;
    for (long i = 0; i < n; i++) {
        if (mz[i] >= lo && mz[i] <= hi && (k < 0 || it[i] > it[k])) k = i;
    }
    if (k < 0) return false;
    return apex_at(mz, it, n, k, profile, pos, top);
}

// _apex_at: position and height of the peak whose highest point is k
template <class Arr> bool apex_at(const Arr& mz, const Arr& it, long n, long k, bool profile, double& pos, double& top) {
    top = it[k];
    if (top <= 0) return false;
    if (!profile) { pos = mz[k]; return true; }
    long a = k, b = k;
    while (a > 0 && it[a - 1] >= 0.3 * top && it[a - 1] <= it[a]) a--;
    while (b < n - 1 && it[b + 1] >= 0.3 * top && it[b + 1] <= it[b]) b++;
    if (b - a >= 2) {
        const long m = b - a + 1;
        std::vector<double> A(3 * m), y(m), c;
        for (long i = 0; i < m; i++) {
            const double x = mz[a + i] - mz[k], w = it[a + i];   // np.polyfit(..., w=sqrt(it^2))
            A[0 * m + i] = x * x * w;
            A[1 * m + i] = x * w;
            A[2 * m + i] = w;
            y[i] = std::log(std::max(it[a + i], 1e-30)) * w;
        }
        if (lstsq(A, y, m, 3, c) && c[0] < 0) {
            const double p = -c[1] / (2 * c[0]);
            if (std::fabs(p) <= mz[b] - mz[a]) { pos = mz[k] + p; return true; }
        }
        std::vector<double> mw(m), ww(m);   // np.sum(mz * ww) / np.sum(ww): pairwise sums
        for (long i = a; i <= b; i++) { mw[i - a] = mz[i] * it[i]; ww[i - a] = it[i]; }
        pos = ms::e2::np_sum(mw.data(), (size_t)m) / ms::e2::np_sum(ww.data(), (size_t)m);
        return true;
    }
    if (0 < k && k < n - 1 && it[k - 1] > 0 && it[k + 1] > 0) {
        const double la = std::log(it[k - 1]), lb = std::log(it[k]), lc = std::log(it[k + 1]);
        const double den = la - 2 * lb + lc;
        if (den < 0) {
            const double off = 0.5 * (la - lc) / den;
            pos = mz[k] + off * (mz[k + 1] - mz[k]);
            return true;
        }
    }
    pos = mz[k];
    return true;
}


// ms_calibration_fit / ms_calibration_fit2: order < 0 = HPC order by cross validation, >= 0 = that order.
// The means and sums are numpy's (pairwise), so that the results are those of hrms_calib.Calibration
// up to the least squares solutions (QR here, an SVD in numpy).
void fit_full(int model, int order, const double* measured, const double* reference, long n, double* coef,
              int* ncoef, double* err_ppm, double* rms, double* cv_err, double* cv_rms) {
    if (!measured || !reference || !coef || !ncoef || n <= 0) throw std::invalid_argument("ms_calibration_fit: null argument");
    if (model < 0 || model > 5) throw std::invalid_argument("ms_calibration_fit: model must be 0..5");
    for (long i = 0; i < n; i++)
        if (!std::isfinite(measured[i]) || !std::isfinite(reference[i]) || !(reference[i] > 0) || !(measured[i] > 0))
            throw std::invalid_argument("ms_calibration_fit: the m/z values must be positive finite numbers");
    const std::vector<double> m(measured, measured + n), r(reference, reference + n);
    const Cal c = fit_model(model, m, r, order);
    std::vector<double> out = {c.x0, c.s};
    out.insert(out.end(), c.coef.begin(), c.coef.end());
    if (model == 4) out.insert(out.end(), c.hpc.begin(), c.hpc.end());
    if (out.size() > 16) throw Error("too many coefficients");
    std::copy(out.begin(), out.end(), coef);
    *ncoef = (int)out.size();
    std::vector<double> e2v(n);
    for (long i = 0; i < n; i++) {
        const double e = (c(m[i]) - r[i]) / r[i] * 1e6;
        if (err_ppm) err_ppm[i] = e;
        e2v[i] = e * e;
    }
    if (rms) *rms = std::sqrt(ms::e2::np_sum(e2v.data(), (size_t)n) / n);
    // cross validation: each calibrant predicted by the fit without it
    double cv = NaN;
    if (cv_err) for (long i = 0; i < n; i++) cv_err[i] = NaN;
    if (n > n_params(model)) {
        std::vector<double> errs(n);
        bool ok = true;
        for (long i = 0; i < n && ok; i++) {
            std::vector<double> mk, rk;
            for (long j = 0; j < n; j++) if (j != i) { mk.push_back(m[j]); rk.push_back(r[j]); }
            try {
                // the HPC order stays the one chosen for the full fit (hpc_order=self.hpc_order)
                const Cal ci = fit_model(model, mk, rk, c.order);
                const double e = (ci(m[i]) - r[i]) / r[i] * 1e6;
                if (!std::isfinite(e)) ok = false;
                errs[i] = e;
            } catch (const std::exception&) { ok = false; }
        }
        if (ok) {
            for (long i = 0; i < n; i++) e2v[i] = errs[i] * errs[i];
            cv = std::sqrt(ms::e2::np_sum(e2v.data(), (size_t)n) / n);
            if (cv_err) std::copy(errs.begin(), errs.end(), cv_err);
        }
    }
    if (cv_rms) *cv_rms = cv;
}

// hrms_calib.peak_near: of the peak maxima within target +- d (at least min_rel of the tallest
// there), the one nearest to target; its apex
bool peak_near(const double* mz, const double* it, long n, double target, double d, bool profile, double min_rel,
               double& pos, double& top) {
    const long a = std::lower_bound(mz, mz + n, target - d) - mz;    // searchsorted left
    const long b = std::upper_bound(mz, mz + n, target + d) - mz;    // searchsorted right
    if (b <= a) return false;
    std::vector<long> k;
    if (!profile) {
        for (long i = a; i < b; i++) k.push_back(i);
    } else {
        const long lo = std::max(a - 1, 0L), hi = std::min(b + 1, n);
        if (hi - lo < 3) {
            for (long i = a; i < b; i++) k.push_back(i);
        } else {
            for (long i = lo + 1; i < hi - 1; i++)
                if (it[i] > it[i - 1] && it[i] >= it[i + 1] && i >= a && i < b) k.push_back(i);
            if (k.empty()) {   // the window on the flank of a peak: its highest point (first maximum)
                long j = a;
                for (long i = a + 1; i < b; i++) if (it[i] > it[j]) j = i;
                k.push_back(j);
            }
        }
    }
    std::vector<long> kk;
    for (long i : k) if (it[i] > 0) kk.push_back(i);
    if (kk.empty()) return false;
    double mx = it[kk[0]];
    for (long i : kk) if (it[i] > mx) mx = it[i];
    long j = -1;
    double best = 0;
    for (long i : kk) {
        if (!(it[i] >= min_rel * mx)) continue;
        const double dd = std::fabs(mz[i] - target);
        if (j < 0 || dd < best) { j = i; best = dd; }
    }
    if (j < 0) return false;
    return apex_at(mz, it, n, j, profile, pos, top);
}

}  // namespace
}  // namespace ms

extern "C" {

MS_API int ms_calibration_fit(int model, int order, const double* measured, const double* reference, long n,
                              double* coef, int* ncoef, double* err_ppm, double* rms, double* cv_rms) {
    return ms::guarded([&] {
        ms::fit_full(model, order <= 0 ? -1 : order, measured, reference, n, coef, ncoef, err_ppm, rms, nullptr, cv_rms);
    });
}

MS_API int ms_calibration_apply(int model, int order, const double* coef, int ncoef, double m_lo, double m_hi,
                                const double* mz, double* out, long n) {
    (void)order;
    return ms::guarded([&] {
        if (!coef || (n > 0 && (!mz || !out))) throw std::invalid_argument("ms_calibration_apply: null argument");
        if (model < 0 || model > 5) throw std::invalid_argument("ms_calibration_apply: model must be 0..5");
        const int need = model == 4 ? 6 : (model == 3 ? 5 : (model == 5 ? 3 : model + 4));
        if (ncoef < need) throw std::invalid_argument("ms_calibration_apply: too few coefficients for this model");
        if (ncoef > 16) throw std::invalid_argument("ms_calibration_apply: more than 16 coefficients");
        ms::Cal c;
        c.model = model;
        c.x0 = coef[0];
        c.s = coef[1];
        c.lo = m_lo;
        c.hi = m_hi;
        const int nc = model == 4 ? 3 : need - 2;
        c.coef.assign(coef + 2, coef + 2 + nc);
        if (model == 4) {
            c.hpc.assign(coef + 5, coef + ncoef);
            c.order = ncoef - 6;
        }
        for (long i = 0; i < n; i++) out[i] = c(mz[i]);
    });
}

MS_API int ms_find_calibrants(const double* mz, const double* it, long n, int profile, const double* reference, long nref,
                              double tol_ppm, double min_rel, double* measured, double* relint, int* found) {
    return ms::guarded([&] {
        if (!reference || !measured || !relint || !found || (n > 0 && (!mz || !it)))
            throw std::invalid_argument("ms_find_calibrants: null argument");
        if (n < 0 || nref < 0) throw std::invalid_argument("ms_find_calibrants: negative length");
        double top = 0;
        for (long i = 0; i < n; i++) if (it[i] > top) top = it[i];
        for (long q = 0; q < nref; q++) {
            const double r = reference[q];
            measured[q] = 0.0;
            relint[q] = 0.0;
            found[q] = 0;
            if (n == 0 || r < mz[0] || r > mz[n - 1]) continue;
            const double d = r * tol_ppm * 1e-6;
            double pos = 0, inten = 0;
            const bool have = ms::peak_position(mz, it, n, r - d, r + d, profile != 0, pos, inten);
            const bool ok = have && top > 0 && inten >= min_rel * top;
            if (ok) {
                measured[q] = pos;
                relint[q] = 100.0 * inten / top;
                found[q] = 1;
            }
        }
        // two reference ions closer than the search window can find the same peak: it belongs
        // to the nearer one (hrms_calib.find_calibrants; the list is in increasing m/z there)
        for (long a = 0; a < nref; a++) {
            if (!found[a]) continue;
            const double ea = std::fabs(measured[a] - reference[a]) / reference[a];
            for (long b = 0; b < nref; b++) {
                if (b == a || !found[b] || std::fabs(measured[a] - measured[b]) >= 1e-9) continue;
                const double eb = std::fabs(measured[b] - reference[b]) / reference[b];
                if (eb < ea || (eb == ea && (reference[b] < reference[a] || (reference[b] == reference[a] && b < a)))) {
                    // a is the farther one (or a repeat of b in the list)
                    measured[a] = 0.0;
                    relint[a] = 0.0;
                    found[a] = 0;
                    break;
                }
            }
        }
    });
}

}

// ---- 3.1 additions (agent P) --------------------------------------------------------
extern "C" {

MS_API int ms_calibration_fit2(int model, int order, const double* measured, const double* reference, long n,
                               double* coef, int* ncoef, double* err_ppm, double* rms, double* cv_err, double* cv_rms) {
    return ms::guarded([&] {
        ms::fit_full(model, order < 0 ? -1 : order, measured, reference, n, coef, ncoef, err_ppm, rms, cv_err, cv_rms);
    });
}

MS_API int ms_peak_position(const double* mz, const double* it, long n, long stride, double lo, double hi, int profile,
                            double* pos, double* height) {
    bool none = false;
    const int rc = ms::guarded([&] {
        if (!pos || !height || n < 0 || stride < 1 || (n > 0 && (!mz || !it))) throw std::invalid_argument("ms_peak_position: bad arguments");
        *pos = 0.0;
        *height = 0.0;
        double p = 0, h = 0;
        const ms::Col cm{mz, stride}, ci{it, stride};
        if (!ms::peak_position(cm, ci, n, lo, hi, profile != 0, p, h)) { none = true; return; }
        *pos = p;
        *height = h;
    });
    return (rc == MS_OK && none) ? MS_NOT_FOUND : rc;
}

MS_API int ms_peak_near(const double* mz, const double* it, long n, double target, double d, int profile,
                        double min_rel, double* pos, double* height) {
    bool none = false;
    const int rc = ms::guarded([&] {
        if (!pos || !height || n < 0 || (n > 0 && (!mz || !it))) throw std::invalid_argument("ms_peak_near: bad arguments");
        if (std::isnan(target) || std::isnan(d)) throw std::invalid_argument("ms_peak_near: NaN target or window");
        *pos = 0.0;
        *height = 0.0;
        double p = 0, h = 0;
        if (n == 0 || !ms::peak_near(mz, it, n, target, d, profile != 0, min_rel, p, h)) { none = true; return; }
        *pos = p;
        *height = h;
    });
    return (rc == MS_OK && none) ? MS_NOT_FOUND : rc;
}

MS_API double ms_profile_apex(const double* mz, const double* it, long n, long i) {
    return ms::guarded_value(ms::NaN, [&] {
        if (!mz || !it || n <= 0 || i < 0 || i >= n) return ms::NaN;
        double p = 0, h = 0;
        if (!ms::apex_at(mz, it, n, i, true, p, h)) return mz[i];
        return p;
    });
}

MS_API int ms_find_calibrants2(const double* mz_p, const double* it_p, long n, long stride, int profile,
                               const double* reference, long nref, double tol_ppm, double min_rel, double* measured,
                               double* height, double* relint, int* found) {
    return ms::guarded([&] {
        if (!reference || !measured || !height || !relint || !found || n < 0 || nref < 0 || stride < 1 ||
            (n > 0 && (!mz_p || !it_p)))
            throw std::invalid_argument("ms_find_calibrants2: bad arguments");
        const ms::Col mz{mz_p, stride}, it{it_p, stride};
        double top = 0;   // it.max() (0 for an empty spectrum)
        for (long i = 0; i < n; i++) if (i == 0 || it[i] > top) top = it[i];
        for (long q = 0; q < nref; q++) {
            const double r = reference[q];
            measured[q] = 0.0;
            height[q] = 0.0;
            relint[q] = 0.0;
            found[q] = 0;
            if (n == 0 || r < mz[0] || r > mz[n - 1]) continue;
            const double d = r * tol_ppm * 1e-6;
            double pos = 0, inten = 0;
            const bool have = ms::peak_position(mz, it, n, r - d, r + d, profile != 0, pos, inten);
            if (!have) { inten = 0.0; continue; }
            const bool ok = top > 0 && inten >= min_rel * top;
            if (ok) {
                measured[q] = pos;
                height[q] = inten;
                relint[q] = 100.0 * inten / top;
                found[q] = 1;
            }
        }
    });
}

}

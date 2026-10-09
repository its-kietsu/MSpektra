// PDA (UV/Vis) computations of the LCMS window: port of lcms_pda.PDAData (chromatogram at a
// wavelength and band, max plot, spectrum at a time or over a range minus a background).
// The absorbance matrix is the float32 matrix of the Python object (n_times x n_wl, row major,
// mAU); the means are numpy's in float32: the band of a chromatogram and the spectrum over a time
// range are both sums of whole columns / rows one after the other (the layouts numpy gives the
// fancy indexed copies), divided as numpy does (in double, rounded to float32).
#include "common.h"
#include <algorithm>
#include <cstdint>
#include <limits>
#if defined(__SSE2__)
#include <emmintrin.h>
#endif

namespace ms {
namespace {

// first index of the smallest |v - x| (np.argmin; a NaN wins)
long nearest(const double* v, long n, double x) {
    long k = 0;
    double best = std::fabs(v[0] - x);
    for (long i = 1; i < n && !std::isnan(best); i++) {
        const double d = std::fabs(v[i] - x);
        if (std::isnan(d) || d < best) { best = d; k = i; }
    }
    return k;
}

// _mean_spectrum(t0, t1) as double (t1 NaN: a single spectrum)
void mean_spectrum(const float* a, long nt, long nw, const double* times, double interval_s, double t0, double t1,
                   double* out) {
    if (std::isnan(t1) || std::fabs(t1 - t0) < interval_s / 60.0) {
        const long i = nearest(times, nt, t0);
        for (long j = 0; j < nw; j++) out[j] = (double)a[i * nw + j];
        return;
    }
    const double lo = std::min(t0, t1), hi = std::max(t0, t1);
    std::vector<float> s(nw, 0.0f);
    long k = 0;
    for (long i = 0; i < nt; i++) {
        if (!(times[i] >= lo && times[i] <= hi)) continue;
        const float* row = a + i * nw;
        for (long j = 0; j < nw; j++) s[j] += row[j];   // add.reduce over axis 0: row after row, float32
        k++;
    }
    if (k == 0) {
        const long i = nearest(times, nt, (lo + hi) / 2);
        for (long j = 0; j < nw; j++) out[j] = (double)a[i * nw + j];
        return;
    }
    for (long j = 0; j < nw; j++) out[j] = (double)(float)((double)s[j] / (double)k);
}

void check(const float* a, long nt, long nw, const char* fn) {
    if (!a || nt <= 0 || nw <= 0) throw std::invalid_argument(std::string(fn) + ": empty or missing absorbance matrix");
    if ((double)nt * (double)nw > 4e9) throw std::invalid_argument(std::string(fn) + ": matrix too large");
}

}  // namespace
}  // namespace ms

extern "C" {

// ---- 3.1 additions (agent P) --------------------------------------------------------

MS_API int ms_pda_chromatogram(const float* a, long nt, long nw, const double* wl, double wl0, double bw, float* out) {
    return ms::guarded([&] {
        ms::check(a, nt, nw, "ms_pda_chromatogram");
        if (!wl || !out) throw std::invalid_argument("ms_pda_chromatogram: null argument");
        if (std::isnan(wl0) || std::isnan(bw)) throw std::invalid_argument("ms_pda_chromatogram: NaN wavelength or band");
        // _wl_index: the wavelengths within wl0 +- bw / 2, else the nearest one
        const double half = std::max(bw, 0.0) / 2.0;
        std::vector<long> sel;
        for (long j = 0; j < nw; j++) if (std::fabs(wl[j] - wl0) <= half + 1e-9) sel.push_back(j);
        if (sel.empty()) sel.push_back(ms::nearest(wl, nw, wl0));
        const long ns = (long)sel.size();
        // A[:, sel].mean(axis=1): the fancy index gives a Fortran ordered copy, so numpy adds its
        // columns one after the other (not the pairwise sum of a contiguous row)
        for (long i = 0; i < nt; i++) {
            const float* r = a + i * nw;
            float s = 0.0f;
            for (long k = 0; k < ns; k++) s += r[sel[k]];
            out[i] = (float)((double)s / (double)ns);
        }
    });
}

MS_API int ms_pda_max_plot(const float* a, long nt, long nw, float* out) {
    return ms::guarded([&] {
        ms::check(a, nt, nw, "ms_pda_max_plot");
        if (!out) throw std::invalid_argument("ms_pda_max_plot: null argument");
        for (long i = 0; i < nt; i++) {   // A.max(axis=1): a NaN wins
            const float* r = a + i * nw;
            float m = r[0];
            int nan = 0;
            long j = 0;
#if defined(__SSE2__)
            // four lanes; max is exact, so the order does not matter (NaN: flagged apart)
            __m128 mv = _mm_set1_ps(r[0]), nv = _mm_setzero_ps();
            for (; j + 4 <= nw; j += 4) {
                const __m128 v = _mm_loadu_ps(r + j);
                nv = _mm_or_ps(nv, _mm_cmpunord_ps(v, v));
                mv = _mm_max_ps(v, mv);
            }
            float lanes[4];
            _mm_storeu_ps(lanes, mv);
            nan = _mm_movemask_ps(nv) != 0;
            for (int q = 0; q < 4; q++) m = lanes[q] > m ? lanes[q] : m;
#endif
            for (; j < nw; j++) {
                const float v = r[j];
                nan |= (v != v);
                m = v > m ? v : m;
            }
            out[i] = nan ? std::numeric_limits<float>::quiet_NaN() : m;
        }
    });
}

// spectrum(t0, t1, bg): t1 NaN = the spectrum nearest to t0; bg_mode 0 none, 1 the spectrum
// nearest to b0, 2 the mean over b0 .. b1
MS_API int ms_pda_spectrum(const float* a, long nt, long nw, const double* times, double interval_s,
                           double t0, double t1, int bg_mode, double b0, double b1, double* out) {
    return ms::guarded([&] {
        ms::check(a, nt, nw, "ms_pda_spectrum");
        if (!times || !out) throw std::invalid_argument("ms_pda_spectrum: null argument");
        if (std::isnan(t0) || (bg_mode > 0 && std::isnan(b0)) || (bg_mode == 2 && std::isnan(b1)))
            throw std::invalid_argument("ms_pda_spectrum: NaN time");
        ms::mean_spectrum(a, nt, nw, times, interval_s, t0, t1, out);
        if (bg_mode == 1 || bg_mode == 2) {
            std::vector<double> bg(nw);
            ms::mean_spectrum(a, nt, nw, times, interval_s, b0, bg_mode == 2 ? b1 : std::numeric_limits<double>::quiet_NaN(),
                              bg.data());
            for (long j = 0; j < nw; j++) out[j] = out[j] - bg[j];
        }
    });
}

}

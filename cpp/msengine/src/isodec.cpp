// IsoDec through UniDec's isodeclib (dynamic library): port of
// ms_deconv.isodec_deconvolute and of the Python runtime around the library
// (IsoDec/runtime.py batch_process_spectrum, c_interface.py process_spectrum,
// match.py MatchedCollection / MatchedMass).
//
// Flow: the spectrum is restricted (minimum intensity and baseline of
// run_method), centroided as IsoDec does (local maxima within +-window
// points above threshold x the top), handed to process_spectrum() of the
// library with the IsoSettings structure of c_interface.py (the default
// model built into the library, analyte type "Peptide"), and the matched
// peaks are merged into masses as MatchedCollection.add_pk_to_masses does.
// The result lists the monoisotopic masses with their total intensity, charge
// states and apex m/z, and a stick spectrum of narrow Gaussians.
#include "e2_internal.h"
#include "f2_internal.h"
#include "fileio.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <memory>
#include <string>
#include <vector>
#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

using ms::e2::Spec;
using namespace ms::unidec;

namespace {

// ============================================================ the library's structures (c_interface.py)
struct MPStruct {
    float mz;
    int z;
    float monoiso;
    float peakmass;
    float avgmass;
    float area;
    float peakint;
    int matchedindsiso[64];
    int matchedindsexp[64];
    float isomz[64];
    float isodist[64];
    float isomass[64];
    float monoisos[16];
    int startindex;
    int endindex;
    float score;
    int realisolength;
};

struct IDSettings {
    int phaseres;
    int verbose;
    int peakwindow;
    float peakthresh;
    int minpeaks;
    float css_thresh;
    float matchtol;
    int maxshift;
    float mzwindow[2];
    float plusoneintwindow[2];
    int knockdown_rounds;
    float min_score_diff;
    float minareacovered;
    int isolength;
    double mass_diff_c;
    float adductmass;
    int minusoneaszero;
    float isotopethreshold;
    float datathreshold;
    float zscore_threshold;
};

typedef int (*process_spectrum_fn)(const double* cmz, const float* cint, int n, const char* fname, MPStruct* out,
                                   IDSettings settings, const char* type);

// IsoDecConfig defaults (unidecstructure.IsoDecConfig)
struct IsoConfig {
    double meanpeakspacing_thresh = 0.01;
    int background_subtraction = 0;
    double mass_diff_c = 1.0033;
    int peakwindow = 80;
    int phaseres = 8;
    double matchtol = 5;
    int minpeaks = 3;
    double peakthresh = 0.0001;
    double css_thresh = 0.7;
    int maxshift = 3;
    double mzwindowlb = -1.05, mzwindowub = 4.05;
    double plusoneintwindowlb = 0.1, plusoneintwindowub = 0.6;
    int knockdown_rounds = 5;
    double min_score_diff = 0.1;
    double minareacovered = 0.20;
    int minusoneaszero = 1;
    double isotopethreshold = 0.01;
    double datathreshold = 0.05;
    double zscore_threshold = 0.95;
    double adductmass = 1.007276467;
    int verbose = 0;
};

IDSettings to_settings(const IsoConfig& c) {
    IDSettings s;
    memset(&s, 0, sizeof s);
    s.phaseres = c.phaseres;
    s.verbose = c.verbose;
    s.peakwindow = c.peakwindow;
    s.peakthresh = (float)c.peakthresh;
    s.minpeaks = c.minpeaks;
    s.css_thresh = (float)c.css_thresh;
    s.matchtol = (float)c.matchtol;
    s.maxshift = c.maxshift;
    s.mzwindow[0] = (float)c.mzwindowlb; s.mzwindow[1] = (float)c.mzwindowub;
    s.plusoneintwindow[0] = (float)c.plusoneintwindowlb; s.plusoneintwindow[1] = (float)c.plusoneintwindowub;
    s.knockdown_rounds = c.knockdown_rounds;
    s.min_score_diff = (float)c.min_score_diff;
    s.minareacovered = (float)c.minareacovered;
    s.isolength = 64;
    s.mass_diff_c = c.mass_diff_c;
    s.adductmass = (float)c.adductmass;
    s.minusoneaszero = c.minusoneaszero;
    s.isotopethreshold = (float)c.isotopethreshold;
    s.datathreshold = (float)c.datathreshold;
    s.zscore_threshold = (float)c.zscore_threshold;
    return s;
}

// ============================================================ loading the library
struct IsoLib {
    void* handle = nullptr;
    process_spectrum_fn process_spectrum = nullptr;
    std::string path;
    ~IsoLib() {
#ifdef _WIN32
        if (handle) FreeLibrary((HMODULE)handle);
#else
        if (handle) dlclose(handle);
#endif
    }
};

std::unique_ptr<IsoLib> load_isodec(const std::string& lib_dir) {
    std::unique_ptr<IsoLib> lib(new IsoLib);
#ifdef _WIN32
    const char* name = "isodeclib.dll";
#elif defined(__APPLE__)
    const char* name = "isodeclib.dylib";
#else
    const char* name = "isodeclib.so";
#endif
    std::vector<std::string> candidates;
    if (!lib_dir.empty()) {
        candidates.push_back(join_path(lib_dir, name));
        candidates.push_back(join_path(join_path(lib_dir, "IsoDec"), name));
        candidates.push_back(join_path(join_path(join_path(lib_dir, "unidec"), "IsoDec"), name));
    }
    candidates.push_back(name);
    std::string tried;
    for (const std::string& c : candidates) {
#ifdef _WIN32
        // the altered search path finds isogen.dll next to it; the path is UTF-8 (wide API:
        // the ANSI one cannot open a folder with characters outside the code page)
        const std::wstring wc = ms::utf8_to_wide(c);
        HMODULE h = LoadLibraryExW(wc.c_str(), nullptr, LOAD_WITH_ALTERED_SEARCH_PATH);
        if (!h) h = LoadLibraryW(wc.c_str());
        if (h) {
            lib->handle = h;
            lib->process_spectrum = reinterpret_cast<process_spectrum_fn>(reinterpret_cast<void*>(GetProcAddress(h, "process_spectrum")));
        }
#else
        void* h = dlopen(c.c_str(), RTLD_NOW | RTLD_LOCAL);
        if (h) {
            lib->handle = h;
            lib->process_spectrum = (process_spectrum_fn)dlsym(h, "process_spectrum");
        }
#endif
        if (lib->handle) {
            lib->path = c;
            if (!lib->process_spectrum) throw std::runtime_error("isodeclib has no process_spectrum: " + c);
            return lib;
        }
        tried += (tried.empty() ? "" : ", ") + c;
    }
    throw std::runtime_error("IsoDec library not found (tried " + tried + ")");
}

// ============================================================ datatools
// fastpeakdetect: local maxima within +- window points above threshold x the top
Spec fast_peakdetect(const Spec& d, int window, double threshold) {
    Spec out;
    long n = (long)d.x.size();
    if (n == 0) return out;
    double maxval = *std::max_element(d.y.begin(), d.y.end());
    if (maxval == 0) return out;
    for (long i = 0; i < n; i++) {
        if (!(d.y[i] > maxval * threshold)) continue;
        long start = std::max(0L, i - window), end = std::min(n, i + window + 1);
        if (start >= end) continue;
        double testmax = d.y[start];
        for (long j = start + 1; j < end; j++) testmax = std::max(testmax, d.y[j]);
        bool first = true;
        for (long j = start; j < i; j++) if (d.y[j] == d.y[i]) { first = false; break; }
        if (d.y[i] == testmax && first) { out.x.push_back(d.x[i]); out.y.push_back(d.y[i]); }
    }
    return out;
}

// datatools.datacompsub (IsoDec version): windowed minimum, Gaussian smoothed
void iso_datacompsub(Spec& d, double buff) {
    long n = (long)d.x.size();
    std::vector<double> mins(n);
    for (long i = 0; i < n; i++) {
        long a = std::max(0L, i - (long)std::fabs(buff)), b = std::min(n, i + (long)std::fabs(buff));
        double m = std::numeric_limits<double>::infinity();
        for (long j = a; j < b; j++) m = std::min(m, d.y[j]);
        mins[i] = b > a ? m : 0.0;
    }
    gaussian_filter_reflect(mins, std::fabs(buff) * 2);
    for (long i = 0; i < n; i++) d.y[i] -= mins[i];
}

// remove_noise_cdata(data, 100, 1.5, "median")
Spec remove_noise_cdata(const Spec& d, long localmin, double factor) {
    long n = (long)d.x.size();
    if (n < localmin) return d;
    std::vector<double> lm;
    for (long i = 0; i < n - localmin; i++) {
        std::vector<double> w(d.y.begin() + i, d.y.begin() + i + localmin);
        lm.push_back(ms::e2::median(std::move(w)));
    }
    double last = lm.empty() ? 0.0 : lm.back();
    while ((long)lm.size() < n) lm.push_back(last);
    // np.convolve(lm, ones(localmin) / localmin, mode="same")
    std::vector<double> noise(n, 0.0);
    long k = localmin;
    long off = (k - 1) / 2;   // 'same' output j = full[j + off]
    for (long j = 0; j < n; j++) {
        long f = j + off;     // full index: sum lm[i] for i in [f - k + 1, f]
        long a = std::max(0L, f - k + 1), b = std::min(n - 1, f);
        double s = 0;
        for (long i = a; i <= b; i++) s += lm[i];
        noise[j] = s / (double)k;
    }
    Spec out;
    for (long i = 0; i < n; i++) if (d.y[i] - noise[i] * factor > 0) { out.x.push_back(d.x[i]); out.y.push_back(d.y[i]); }
    return out;
}

// datatools.fastnearest on a sorted array
long fastnearest(const std::vector<double>& a, double t) {
    if (a.empty()) return 0;
    long start = 0, length = (long)a.size() - 1;
    long diff = length - start;
    while (diff > 1) {
        long mid = start + (length - start) / 2;
        if (t < a[mid]) length = mid;
        else if (t == a[mid]) return mid;
        else start = mid;
        diff = length - start;
    }
    return std::fabs(t - a[start]) >= std::fabs(t - a[length]) ? length : start;
}

bool within_ppm(double theo, double exp, double tol) { return std::fabs((theo - exp) / theo * 1e6) <= tol; }

// ============================================================ match.py
struct Dist { std::vector<double> x, y; };   // (mass, intensity) rows

struct MPeak {
    int z = 0;
    float mz = 0, monoiso = 0, peakmass = 0, avgmass = 0, peakint = 0, matchedintensity = 0;
    std::vector<float> monoisos;
    Dist massdist;            // (isomass, isodist) above 0.1 % of the top
    Dist centroids;           // the centroids of the cluster
    Dist decon_centroids;     // centroids moved to the mass axis
};

// numpy pairwise summation of a float32 array
float np_sum_f32(const float* a, size_t n) {
    if (n < 8) { float r = 0; for (size_t i = 0; i < n; i++) r += a[i]; return r; }
    if (n <= 128) {
        float r[8];
        for (int j = 0; j < 8; j++) r[j] = a[j];
        size_t i = 8;
        for (; i + 8 <= n; i += 8) for (int j = 0; j < 8; j++) r[j] += a[i + j];
        float res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        for (; i < n; i++) res += a[i];
        return res;
    }
    size_t n2 = n / 2; n2 -= n2 % 8;
    return np_sum_f32(a, n2) + np_sum_f32(a + n2, n - n2);
}

// find_matched_intensities with max_shift 0, z 1: the tallest centroid within tol ppm (of peakmz) of each isotope
std::vector<double> matched_intensities(const Dist& cent, const Dist& iso, double tol) {
    size_t m = iso.x.size();
    std::vector<double> out(m, 0.0);
    double diff = iso.x[0] * tol * 1e-6;
    for (size_t i = 0; i < m; i++) {
        double lo = iso.x[i] - diff, hi = iso.x[i] + diff;
        for (size_t j = 0; j < cent.x.size(); j++) {
            double mz2 = cent.x[j];
            if (mz2 > hi) break;
            if (mz2 < lo) continue;
            if (cent.y[j] > out[i]) out[i] = cent.y[j];
        }
    }
    return out;
}

// calc_css_from_data (cosine similarity, with the "-1 isotope as zero" term of
// calculate_cosinesimilarity that reads the last matched intensity)
double calc_css(const Dist& cent, const Dist& iso0, double tol) {
    if (cent.x.empty() || iso0.x.empty()) return 0.0;
    std::vector<double> ci = matched_intensities(cent, iso0, tol);
    double cmax = *std::max_element(ci.begin(), ci.end());
    double imax = *std::max_element(iso0.y.begin(), iso0.y.end());
    double norm = cmax / imax;
    double ab = 0, a2 = 0, b2 = 0;
    {   // minusoneaszero: cent_intensities[max_shift + shift - 1] = cent_intensities[-1]
        double a = ci.back();
        a2 += a * a;
    }
    for (size_t i = 0; i < iso0.y.size(); i++) {
        double a = ci[i], b = iso0.y[i] * norm;
        ab += a * b; a2 += a * a; b2 += b * b;
    }
    if (ab == 0 || a2 == 0 || b2 == 0) return 0.0;
    return ab / (std::sqrt(a2) * std::sqrt(b2));
}

// merge_decon_centroids: sorted, neighbours within ppm_tol merged (intensity weighted)
Dist merge_decon_centroids(Dist c, double ppm_tol) {
    size_t n = c.x.size();
    if (n == 0) return c;
    std::vector<size_t> idx(n);
    for (size_t i = 0; i < n; i++) idx[i] = i;
    std::stable_sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return c.x[a] < c.x[b]; });
    Dist out;
    double cm = c.x[idx[0]], ci = c.y[idx[0]];
    double ws = cm * ci, ti = ci;
    for (size_t k = 1; k < n; k++) {
        double nm = c.x[idx[k]], ni = c.y[idx[k]];
        double ppm = std::fabs(cm - nm) / cm * 1e6;
        if (ppm <= ppm_tol) {
            ws += nm * ni; ti += ni;
            cm = ws / ti; ci = ti;
        } else {
            out.x.push_back(cm); out.y.push_back(ci);
            cm = nm; ci = ni; ws = cm * ci; ti = ci;
        }
    }
    out.x.push_back(cm); out.y.push_back(ci);
    return out;
}

// merge_massdist: the distribution shifted and scaled onto the deconvolved centroids (float32 arrays)
void merge_massdist(Dist& md, const Dist& dc, double tol) {
    if (md.x.empty() || dc.x.empty()) return;   // (md.x[0] below; fastnearest of an empty list gives 0)
    double max_intensity = 0, diff_at_max = 0, normfactor = 1;
    for (size_t i = 0; i < md.x.size(); i++) {
        if (md.y[i] > max_intensity) {
            long idx = fastnearest(dc.x, md.x[i]);
            if (within_ppm(dc.x[idx], md.x[i], tol)) {
                max_intensity = md.y[i];
                diff_at_max = std::fabs(dc.x[idx] - md.x[i]);
                normfactor = md.y[i] != 0 ? dc.y[idx] / md.y[i] : 1;
            }
        }
    }
    if (diff_at_max < tol * md.x[0] / 1e6) for (double& v : md.x) v = (double)(float)(v + diff_at_max);
    if (max_intensity != 0) for (double& v : md.y) v = (double)(float)(v * normfactor);
}

struct MMass {
    float monoiso = 0;
    std::vector<float> monoisos;
    float apexintensity = 0, totalintensity = 0, avgmass = 0;
    std::vector<int> zs;
    std::vector<float> mzs, mzints;
    Dist decon_centroids, massdist;
    int totalpeaks = 1;

    explicit MMass(const MPeak& pk) {
        monoiso = pk.monoiso; monoisos = pk.monoisos; apexintensity = pk.peakint;
        mzs.push_back(pk.mz); zs.push_back(pk.z); totalintensity = pk.matchedintensity; mzints.push_back(pk.peakint);
        avgmass = pk.avgmass; decon_centroids = pk.decon_centroids; massdist = pk.massdist;
    }

    bool check_if_match(const MPeak& pk, const IsoConfig& c) const {
        if (!(std::fabs((double)pk.monoiso - (double)monoiso) <= c.maxshift * c.mass_diff_c * 1.1)) return false;
        bool ppm_ok = false;
        for (int mm = -c.maxshift; mm <= c.maxshift; mm++)
            if (within_ppm(monoiso, (double)pk.monoiso + mm * c.mass_diff_c, c.matchtol)) { ppm_ok = true; break; }
        if (!ppm_ok) return false;
        if (decon_centroids.x.empty() || pk.massdist.x.empty()) return false;
        return calc_css(decon_centroids, pk.massdist, c.matchtol) >= c.css_thresh;
    }

    void merge_in_pk(const MPeak& pk, const IsoConfig& c) {
        if (pk.peakint > apexintensity) apexintensity = pk.peakint;
        if (std::find(zs.begin(), zs.end(), pk.z) == zs.end()) { zs.push_back(pk.z); mzs.push_back(pk.mz); mzints.push_back(pk.peakint); }
        Dist all = decon_centroids;
        all.x.insert(all.x.end(), pk.decon_centroids.x.begin(), pk.decon_centroids.x.end());
        all.y.insert(all.y.end(), pk.decon_centroids.y.begin(), pk.decon_centroids.y.end());
        decon_centroids = merge_decon_centroids(all, c.matchtol);
        // monoisotopic mass candidates (float32 arithmetic as in the Python)
        for (float m : pk.monoisos) {
            if (monoisos.empty()) { monoisos.push_back(m); continue; }   // (monoisos[0] of an empty list below)
            std::vector<double> mi(monoisos.begin(), monoisos.end());
            long ci = fastnearest(mi, m);
            if (within_ppm(m, monoisos[ci], c.matchtol)) {
                float inten = pk.matchedintensity;
                monoisos[ci] = monoisos[ci] * totalintensity + m * inten / (totalintensity + inten);
            } else {
                monoisos.push_back(m);
            }
        }
        if (within_ppm(monoiso, pk.monoiso, c.matchtol)) {
            float inten = pk.matchedintensity;
            monoiso = (monoiso * totalintensity + pk.monoiso * inten) / (totalintensity + inten);
        } else {
            double css_new = calc_css(decon_centroids, pk.massdist, c.matchtol);
            double css_old = calc_css(decon_centroids, massdist, c.matchtol);
            if (css_new > css_old) { monoiso = pk.monoiso; massdist = pk.massdist; }
        }
        merge_massdist(massdist, decon_centroids, c.matchtol);
        totalintensity += pk.matchedintensity;
        totalpeaks++;
    }
};

struct Collection {
    std::vector<MMass> masses;
    std::vector<double> monoisos;   // kept sorted, parallel to masses

    void add(const MPeak& pk, const IsoConfig& c) {
        if (masses.empty()) { monoisos.push_back(pk.monoiso); masses.emplace_back(pk); return; }
        // fastwithin_abstol_withnearest
        long idx = fastnearest(monoisos, pk.monoiso);
        double tol = c.maxshift * 1.25;
        std::vector<long> indices;
        if (std::fabs(monoisos[idx] - pk.monoiso) <= tol) {
            indices.push_back(idx);
            long up = idx + 1, lo = idx - 1;
            bool au = true, al = true;
            while (au || al) {
                if (au) {
                    if (up >= (long)monoisos.size()) au = false;
                    else if (std::fabs(monoisos[up] - pk.monoiso) <= tol) { indices.push_back(up); up++; }
                    else au = false;
                }
                if (al) {
                    if (lo < 0) al = false;
                    else if (std::fabs(monoisos[lo] - pk.monoiso) <= tol) { indices.push_back(lo); lo--; }
                    else al = false;
                }
            }
        }
        double nearest_mass = monoisos[idx];
        std::vector<long> matched;
        for (long i : indices) if (masses[i].check_if_match(pk, c)) matched.push_back(i);
        if (!matched.empty()) {
            // several matches: the closest retention time, all equal here -> the first
            masses[matched[0]].merge_in_pk(pk, c);
            return;
        }
        if (pk.monoiso > nearest_mass) idx++;
        if (idx >= (long)monoisos.size()) { monoisos.push_back(pk.monoiso); masses.emplace_back(pk); }
        else { monoisos.insert(monoisos.begin() + idx, pk.monoiso); masses.insert(masses.begin() + idx, MMass(pk)); }
    }
};

}  // namespace

// ============================================================ result object
struct ms_isodec_result {
    std::vector<double> pk_mono, pk_height, pk_mz;
    std::vector<int> pk_z;
    std::vector<std::string> pk_charges;
    std::vector<const char*> pk_charges_c;
    std::vector<double> mass_x, mass_y;
    std::string notes;
};

namespace {

struct IPeak { double mass, height, mz, avg; std::vector<int> zs; };

void run_isodec(const double* mz, const double* it, long n0, const std::string& lib_dir, const Params& p, ms_isodec_result& out) {
    using namespace ms::e2;
    double zlo_d = p.num("z_lo", 1), zhi_d = p.num("z_hi", 50);
    if (!std::isfinite(zlo_d) || !std::isfinite(zhi_d) || zlo_d < 1 || zhi_d < zlo_d)
        throw std::invalid_argument("Use whole-number charges in increasing order");
    if (zhi_d > 10000) throw std::invalid_argument("Use charges from 1 to 10000");
    int zlo = (int)zlo_d, zhi = (int)zhi_d;
    double lo_m = p.num("mass_lo", 0), hi_m = p.num("mass_hi", std::numeric_limits<double>::infinity());
    // ---- run_method: restriction, minimum intensity, baseline
    Spec spec = restrict_spec(mz, it, n0, p.num("mz_lo", 0), p.num("mz_hi", 0));
    size_t n = spec.x.size();
    double base = n ? *std::max_element(spec.y.begin(), spec.y.end()) : 0.0;
    double min_int = p.num("min_intensity", 0);
    bool pct = p.inum("min_intensity_pct", 0) != 0 || p.str("min_intensity_unit") == "pct";
    double thr = min_int > 0 ? (pct ? min_int / 100.0 * base : min_int) : 0.0;
    std::vector<char> keep;
    long npk = 0;
    if (thr > 0 && n) { auto pm = peak_mask(spec, thr); keep = pm.first; npk = pm.second; }
    else {
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
            double pc = base > 0 ? 100.0 * thr / base : 0.0;
            snprintf(b, sizeof b, ", peaks from %s up (%.3g %% of the base peak %s): %ld peaks kept whole, %ld of %ld points",
                     fmt_int(thr).c_str(), pc, fmt_int(base).c_str(), npk, points_used, (long)n);
            summary += b;
        } else {
            snprintf(b, sizeof b, ", no minimum intensity (base peak %s)", fmt_int(base).c_str());
            summary += b;
        }
    }
    if (p.inum("baseline", 0)) subtract_baseline(spec, p.num("baseline_width", 15.0) > 0 ? p.num("baseline_width", 15.0) : 15.0);
    if (thr > 0) {
        bool anypos = false;
        for (size_t i = 0; i < n; i++) { if (!keep[i]) spec.y[i] = 0.0; else if (spec.y[i] > 0) anypos = true; }
        if (npk < 1 || points_used < 5 || !anypos)
            throw std::invalid_argument("No data above the minimum intensity (" + fmt_int(thr) + "; the base peak in the m/z range is " +
                                        fmt_int(base) + "). Lower the minimum intensity.");
    }
    // ---- isodec_deconvolute
    {
        Spec pos;
        for (size_t i = 0; i < n; i++) if (spec.y[i] > 0) { pos.x.push_back(spec.x[i]); pos.y.push_back(spec.y[i]); }
        spec = pos;
        n = spec.x.size();
    }
    if (n < 5) throw std::invalid_argument("Too few data points in the m/z range");
    IsoConfig c;
    c.phaseres = (int)p.inum("phaseres", p.inum("iso_phaseres", 8));
    if (c.phaseres <= 0) c.phaseres = 8;
    {
        double a = p.num("adduct_mass", 0);
        double sign = p.num("sign", 1);
        if (a == 0 || !std::isfinite(a)) a = PROTON * (sign < 0 ? -1 : 1);
        c.adductmass = a;
    }
    auto getp = [&](const char* key, bool& found) -> double {
        found = p.has(std::string("iso_") + key) || p.has(key);
        return p.has(std::string("iso_") + key) ? p.num(std::string("iso_") + key, 0) : p.num(key, 0);
    };
    bool f;
    double v;
    v = getp("css_thresh", f); if (f) c.css_thresh = v;
    v = getp("matchtol", f); if (f) c.matchtol = v;
    auto to_int = [](double x) { return std::isfinite(x) ? (int)std::max(-1e9, std::min(1e9, x)) : 0; };   // ((int)NaN is undefined)
    v = getp("minpeaks", f); if (f) c.minpeaks = to_int(v);
    v = getp("maxshift", f); if (f) c.maxshift = to_int(v);
    v = getp("knockdown_rounds", f); if (f) c.knockdown_rounds = to_int(v);
    v = getp("min_score_diff", f); if (f) c.min_score_diff = v;
    v = getp("minareacovered", f); if (f) c.minareacovered = v;
    v = getp("isotopethreshold", f); if (f) c.isotopethreshold = v;
    v = getp("datathreshold", f); if (f) c.datathreshold = v;
    v = getp("zscore_threshold", f); if (f) c.zscore_threshold = v;
    v = getp("background_subtraction", f); if (f) c.background_subtraction = to_int(v);
    int window = (int)p.inum("iso_window", p.inum("window", 5));
    double threshold = p.num("iso_thresh", p.num("thresh", 0.0001));
    // ---- batch_process_spectrum: centroids
    Spec newdat = spec;
    if (c.background_subtraction > 0) iso_datacompsub(newdat, c.background_subtraction);
    Spec cent = n >= 3 ? fast_peakdetect(newdat, window, threshold * 0.1) : Spec();
    if (cent.x.size() >= 2) {
        double med = median_diff(cent.x.data(), cent.x.size());
        if (med <= c.meanpeakspacing_thresh) cent = remove_noise_cdata(cent, 100, 1.5);
    }
    // ---- the library
    std::unique_ptr<IsoLib> lib = load_isodec(lib_dir);
    size_t nc = cent.x.size();
    std::vector<MPeak> peaks_found;
    if (nc > 0) {
        std::vector<float> cint(nc);
        for (size_t i = 0; i < nc; i++) cint[i] = (float)cent.y[i];
        std::vector<MPStruct> elems(nc);
        memset(elems.data(), 0, nc * sizeof(MPStruct));
        IDSettings settings = to_settings(c);
        int nmatched = lib->process_spectrum(cent.x.data(), cint.data(), (int)nc, nullptr, elems.data(), settings, "Peptide");
        if (nmatched < 0 || nmatched > (int)nc) throw std::runtime_error("IsoDec returned an invalid number of peaks");
        for (int k = 0; k < nmatched; k++) {
            const MPStruct& q = elems[k];
            if (q.z == 0) continue;
            if (q.startindex < 0 || q.endindex < 0 || q.endindex <= q.startindex)
                throw std::runtime_error("Invalid start or end index in matched peak.");
            int iend = std::min(q.endindex, (int)nc - 1);   // Python slices centroids[start:end + 1]
            if (getenv("MS_DEBUG")) fprintf(stderr, "isodec peak z %d mz %.4f mono %.4f start %d end %d peakint %g\n", q.z, q.mz, q.monoiso, q.startindex, q.endindex, q.peakint);
            MPeak pk;
            pk.z = q.z; pk.mz = q.mz; pk.monoiso = q.monoiso; pk.peakmass = q.peakmass; pk.avgmass = q.avgmass; pk.peakint = q.peakint;
            for (int i = 0; i < 16; i++) if (q.monoisos[i] > 0) pk.monoisos.push_back(q.monoisos[i]);
            float dmax = q.isodist[0];
            for (int i = 1; i < 64; i++) dmax = std::max(dmax, q.isodist[i]);
            std::vector<float> kept_d;
            for (int i = 0; i < 64; i++) if (q.isodist[i] > dmax * 0.001f) {
                kept_d.push_back(q.isodist[i]);
                pk.massdist.x.push_back(q.isomass[i]); pk.massdist.y.push_back(q.isodist[i]);
            }
            pk.matchedintensity = np_sum_f32(kept_d.data(), kept_d.size());
            for (int i = q.startindex; i <= iend; i++) {
                pk.centroids.x.push_back(cent.x[i]); pk.centroids.y.push_back(cent.y[i]);
                pk.decon_centroids.x.push_back(cent.x[i] * q.z - c.adductmass * q.z);
                pk.decon_centroids.y.push_back(cent.y[i]);
            }
            peaks_found.push_back(pk);
        }
    }
    Collection coll;
    for (const MPeak& pk : peaks_found) coll.add(pk, c);
    // ---- peaks of the result
    std::vector<IPeak> peaks;
    for (const MMass& mm : coll.masses) {
        std::vector<int> zs;
        for (int z : mm.zs) if (zlo <= z && z <= zhi && std::find(zs.begin(), zs.end(), z) == zs.end()) zs.push_back(z);
        std::sort(zs.begin(), zs.end());
        if (zs.empty() || !(lo_m <= mm.monoiso && mm.monoiso <= hi_m)) continue;
        double inten = mm.totalintensity != 0 ? mm.totalintensity : mm.apexintensity;
        peaks.push_back({(double)mm.monoiso, inten, (double)mm.mzs[0], (double)mm.avgmass, zs});
    }
    std::stable_sort(peaks.begin(), peaks.end(), [](const IPeak& a, const IPeak& b) { return a.mass < b.mass; });
    double top = 1e-12;
    for (const IPeak& q : peaks) top = std::max(top, q.height);
    double pthr = p.num("peak_threshold", p.num("peak_thresh", 0));
    std::vector<IPeak> kept;
    for (const IPeak& q : peaks) if (q.height >= pthr * top) kept.push_back(q);
    peaks = kept;
    // ---- stick spectrum: narrow Gaussians around each mass
    double step = p.num("mass_step", 0.01);
    if (!(step > 0)) step = 0.01;
    step = std::min(step, 0.05);
    if (!peaks.empty()) {
        double span = peaks.back().mass - peaks.front().mass + 15;
        step = std::max(step, span / 400000.0);
        std::vector<double> grid = np_arange(peaks.front().mass - 5, peaks.back().mass + 10, step);
        std::vector<double> y(grid.size(), 0.0);
        double sig = std::max(3 * step, 0.02);
        for (const IPeak& q : peaks) {
            size_t a = std::lower_bound(grid.begin(), grid.end(), q.mass - 6 * sig) - grid.begin();
            size_t b = std::lower_bound(grid.begin(), grid.end(), q.mass + 6 * sig) - grid.begin();
            for (size_t i = a; i < b; i++) {
                double d = (grid[i] - q.mass) / sig;
                y[i] = std::max(y[i], q.height * std::exp(-0.5 * d * d));
            }
        }
        out.mass_x = grid; out.mass_y = y;
    }
    for (const IPeak& q : peaks) {
        out.pk_mono.push_back(q.mass); out.pk_height.push_back(q.height); out.pk_mz.push_back(q.mz);
        out.pk_z.push_back(q.zs.front());
        std::string zs;
        for (size_t i = 0; i < q.zs.size(); i++) zs += (i ? ", " : "") + std::to_string(q.zs[i]);
        out.pk_charges.push_back(zs);
    }
    for (const std::string& s : out.pk_charges) out.pk_charges_c.push_back(s.c_str());
    char b[256];
    snprintf(b, sizeof b, "IsoDec: %ld monoisotopic mass%s from isotope distributions, charge %d to %d", (long)peaks.size(),
             peaks.size() != 1 ? "es" : "", zlo, zhi);
    out.notes = b;
    if (!summary.empty()) out.notes += "; " + summary;
}

}  // namespace

extern "C" {

MS_API int ms_isodec(const double* mz, const double* it, long n, const char* lib_dir, const char* params, ms_isodec_result** out) {
    if (out) *out = nullptr;
    return ms::guarded([&] {
        if (!mz || !it || n < 0 || !out) throw std::invalid_argument("ms_isodec: bad arguments");
        std::unique_ptr<ms_isodec_result> r(new ms_isodec_result);
        Params p = parse_params(params);
        run_isodec(mz, it, n, lib_dir ? lib_dir : "", p, *r);
        *out = r.release();
    });
}

MS_API void ms_isodec_free(ms_isodec_result* r) { delete r; }

MS_API long ms_isodec_peaks(const ms_isodec_result* r, const double** mono, const double** height,
                            const double** mz_apex, const int** z, const char*** charges) {
    if (!r) return 0;
    if (mono) *mono = r->pk_mono.data();
    if (height) *height = r->pk_height.data();
    if (mz_apex) *mz_apex = r->pk_mz.data();
    if (z) *z = r->pk_z.data();
    if (charges) *charges = const_cast<const char**>(r->pk_charges_c.data());
    return (long)r->pk_mono.size();
}
MS_API long ms_isodec_mass(const ms_isodec_result* r, const double** mass, const double** it) {
    if (!r) return 0;
    if (mass) *mass = r->mass_x.data();
    if (it) *it = r->mass_y.data();
    return (long)r->mass_x.size();
}
MS_API const char* ms_isodec_notes(const ms_isodec_result* r) { return r ? r->notes.c_str() : ""; }

}

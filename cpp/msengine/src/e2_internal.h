// Internal helpers of the peak picking and maximum entropy code (agent E2).
// Everything here is private to peaks.cpp / maxent.cpp; the spectrum
// helpers (peak width, resolution) are local versions, independent of
// spectrum.cpp, so the files can be merged by copying.
#pragma once
#include "common.h"
#include <vector>
#include <cstddef>
#include <string>
#include <utility>

namespace ms {
namespace e2 {

constexpr double PROTON = 1.007276467;

struct Peak {
    double mass = 0, apex = 0, height = 0, area = 0;
    int n_iso = 0;
};

// ---- small numerics (peaks.cpp)
double np_sum(const double* a, size_t n);                 // numpy pairwise summation order
double median(std::vector<double> v);                     // numpy median (mean of the middle two)
double median_diff(const double* x, size_t n);            // median of the point spacing
double py_round(double v);                                // Python round(): half to even

// Local maxima of a mass spectrum within +-window (Da) above threshold x the
// top: mass = centroid of the top half, apex = grid point, area = sum over
// +-window/2 (ms_deconv.pick_mass_peaks). Sorted by mass.
std::vector<Peak> pick_mass_peaks(const double* x, const double* y, size_t n, double window, double threshold);

// Isotope resolved spectrum -> one entry per species (ms_deconv.group_isotopes).
std::vector<Peak> group_isotopes(const double* x, const double* y, size_t n, double threshold, double spacing = 1.00235);

// Apex of the peak near the grid point m0 (ms_deconv.refine_peak_mass).
double refine_peak_mass(const double* x, const double* y, size_t n, double m0);

// Weighted least squares of a small dense system (Householder QR); A is m x n
// row major, returns x (n). Columns of A are used as they are.
void lstsq_small(const std::vector<double>& A, const std::vector<double>& b, int m, int n, std::vector<double>& x);

// FWHM of peak i with interpolated half height crossings; <= 0 when undefined
double fwhm_at(const double* x, const double* y, size_t n, size_t i);
double estimate_peak_width(const double* x, const double* y, size_t n);     // FWHM (m/z) of the tallest peak
double estimate_resolution(const double* x, const double* y, size_t n, int npk = 5);

// ---- spectrum preparation of run_method (maxent.cpp; shared with unidec.cpp, isodec.cpp)
struct Spec { std::vector<double> x, y; };
// finite, positive m/z, sorted, restricted to lo..hi when hi > lo (ms_deconv._restrict)
Spec restrict_spec(const double* mz, const double* it, long n, double lo, double hi);
// rolling baseline (morphological opening) of `factor` peak widths (ms_deconv.subtract_baseline)
void subtract_baseline(Spec& s, double factor);
// points of the peaks whose top reaches thr, kept whole; number of such peaks (ms_deconv.peak_mask)
std::pair<std::vector<char>, long> peak_mask(const Spec& s, double thr);
// intensity for the notes: 850, 42000, 1.25e6 (ms_deconv._fmt_int)
std::string fmt_int(double v);

}  // namespace e2
}  // namespace ms

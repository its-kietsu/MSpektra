// Find mass (ms_find_mass.py of LCMS Analysis): the arithmetic of _profile, bit for bit (findmass.cpp).
#pragma once
#include <vector>

namespace ms {

// The window ends of _profile, one row per charge state, one column per grid point (numpy's order):
// centre = ((M + grid[g]) + z * c) / z, lo = centre - hw / z, hi = centre + hw / z.
void mass_profile_keys(double M, double hw, const int* zs, int nz, double c, const double* grid, long ng,
                       std::vector<double>& lo, std::vector<double>& hi);

// The profile of one spectrum (mz, it: n points, n > 0) as _profile computes it: cs = [0, cumsum(it)],
// P[g] = sum over the states (in order) of cs[searchsorted(mz, hi)] - cs[searchsorted(mz, lo)], searchsorted
// 'left' as numpy does it (its binary search, also for an m/z axis that is not sorted). cs, ilo, ihi: scratch.
void mass_profile_row(const double* mz, const double* it, long n, const std::vector<double>& lo,
                      const std::vector<double>& hi, int nz, long ng, std::vector<double>& cs, std::vector<long>& ilo,
                      std::vector<long>& ihi, double* P);

// The m/z calibration of ms_calibration_apply on n values (calib.cpp; out may be mz): ms_mass_profile_cal
// calibrates each spectrum so, as the window does (hrms_calib.Calibration calls ms_calibration_apply).
void calibration_apply_n(int model, const double* coef, int ncoef, double m_lo, double m_hi, const double* mz,
                         double* out, long n);

}  // namespace ms

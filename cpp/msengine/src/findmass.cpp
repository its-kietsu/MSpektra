// Find mass (ms_find_mass.py of LCMS Analysis): the arithmetic of _profile, in numpy's order of operations
// (compiled without fused multiply add, see CMakeLists.txt), so that ms_mass_profile gives the values of the
// Python code bit for bit. The file access is in file.cpp (ms_mass_profile).
#include "findmass.h"

namespace ms {

namespace {

// numpy's "less" for float64 (NaN sorted last)
inline bool npy_less(double a, double b) { return a < b || (b != b && a == a); }

// numpy.searchsorted(arr, key, side='left') (npysort binsearch): the keys in their order, the search range
// of each key narrowed by the previous one
void searchsorted_left(const double* arr, long n, const double* key, long nk, long* ret) {
    if (nk <= 0) return;
    long min_idx = 0, max_idx = n;
    double last = key[0];
    for (long i = 0; i < nk; i++) {
        const double kv = key[i];
        if (npy_less(last, kv)) {
            max_idx = n;
        } else {
            min_idx = 0;
            max_idx = (max_idx < n) ? (max_idx + 1) : n;
        }
        last = kv;
        while (min_idx < max_idx) {
            const long mid = min_idx + ((max_idx - min_idx) >> 1);
            if (npy_less(arr[mid], kv)) min_idx = mid + 1;
            else max_idx = mid;
        }
        ret[i] = min_idx;
    }
}

}  // namespace

void mass_profile_keys(double M, double hw, const int* zs, int nz, double c, const double* grid, long ng,
                       std::vector<double>& lo, std::vector<double>& hi) {
    lo.resize((size_t)nz * (size_t)ng);
    hi.resize((size_t)nz * (size_t)ng);
    for (int k = 0; k < nz; k++) {
        const double z = (double)zs[k];
        const double zc = z * c;    // zz * c
        const double wz = hw / z;   // w / zz
        for (long g = 0; g < ng; g++) {
            const double cen = ((M + grid[g]) + zc) / z;
            lo[(size_t)k * ng + g] = cen - wz;
            hi[(size_t)k * ng + g] = cen + wz;
        }
    }
}

void mass_profile_row(const double* mz, const double* it, long n, const std::vector<double>& lo,
                      const std::vector<double>& hi, int nz, long ng, std::vector<double>& cs, std::vector<long>& ilo,
                      std::vector<long>& ihi, double* P) {
    // np.concatenate([[0.0], np.cumsum(it)])
    cs.resize((size_t)n + 1);
    cs[0] = 0.0;
    if (n > 0) cs[1] = it[0];
    for (long i = 1; i < n; i++) cs[i + 1] = cs[i] + it[i];
    const long nk = (long)nz * ng;
    ihi.resize((size_t)nk);
    ilo.resize((size_t)nk);
    searchsorted_left(mz, n, hi.data(), nk, ihi.data());
    searchsorted_left(mz, n, lo.data(), nk, ilo.data());
    for (long g = 0; g < ng; g++) P[g] = 0.0;
    for (int k = 0; k < nz; k++) {
        const long* a = ihi.data() + (size_t)k * ng;
        const long* b = ilo.data() + (size_t)k * ng;
        for (long g = 0; g < ng; g++) P[g] += cs[a[g]] - cs[b[g]];
    }
}

}  // namespace ms

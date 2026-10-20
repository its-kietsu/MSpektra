// Scans kept as float32 arrays (ShimadzuReader, ArrayReader): several mass chromatograms in one pass.
#pragma once
#include "common.h"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace ms {

// out[w * idx.size() + k] = the xic() of these readers for window w (|mz32 - (float)mz| <= (float)tol in float32,
// the intensities summed in double in the order of the points) in scan idx[k], bit for bit: the points of a
// scan in increasing m/z are found by binary search on the same float32 test (fl(x - m) does not decrease
// with x, so the points that pass are one run), the others are tested one by one; each window adds the same
// points in the same order. Parallel over the scans (each writes its own outputs); polls *stop.
template <class Off>
void xic_multi_f32(const std::vector<float>& mz_, const std::vector<float>& it_, const std::vector<Off>& offsets_,
                   const std::vector<long>& idx, const double* mz, const double* tol, int nw, double* out,
                   const std::atomic<bool>* stop) {
    const long n = (long)idx.size();
    for (long i : idx)
        if (i < 0 || (size_t)i + 1 >= offsets_.size()) throw std::invalid_argument("scan index out of range");
    if (n == 0 || nw <= 0) return;
    std::vector<float> m(nw), t(nw);
    std::vector<char> fin(nw);
    for (int w = 0; w < nw; w++) {
        m[w] = (float)mz[w];
        t[w] = (float)tol[w];
        fin[w] = std::isfinite(m[w]) && std::isfinite(t[w]);
    }
    const int nt = threads();
    (void)nt;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 16) num_threads(nt) if (n > 32)
#endif
    for (long k = 0; k < n; k++) {
        if (stop && stop->load(std::memory_order_relaxed)) continue;
        const size_t a = (size_t)offsets_[idx[k]], b = (size_t)offsets_[idx[k] + 1];
        const float* x = mz_.data() + a;
        const float* y = it_.data() + a;
        const size_t len = b - a;
        bool sorted = true;   // (a NaN m/z makes it false)
        for (size_t j = 1; j < len && sorted; j++) sorted = x[j - 1] <= x[j];
        for (int w = 0; w < nw; w++) {
            const float mw = m[w], tw = t[w];
            size_t j0 = 0, j1 = len;
            if (sorted && fin[w]) {
                j0 = (size_t)(std::partition_point(x, x + len, [mw, tw](float v) { return (float)(v - mw) < -tw; }) - x);
                j1 = (size_t)(std::partition_point(x + j0, x + len, [mw, tw](float v) { return (float)(v - mw) <= tw; }) - x);
            }
            double s = 0;
            for (size_t j = j0; j < j1; j++)
                if (std::fabs(x[j] - mw) <= tw) s += (double)y[j];
            out[(size_t)w * (size_t)n + (size_t)k] = s;
        }
    }
}

}  // namespace ms

// UniDec engine for one dimensional mass spectra, built into msengine: a C++ port of
// udmain.c (run_unidec, SetupInputs, SetupDeconvolution, RunIteration, CheckConvergence,
// SetupOutputs), udcore.c, udtools.c, udstruct.c (SetDefaultConfig, PostImport, LoadConfig)
// and udio.c (SetupZtab) of UniDec 8.2.1. Scores (UD_score.c) and peak detection
// (UD_analysis.c) are in udscore.cpp.
//
// UniDec engine by Michael T. Marty (Marty et al., Anal. Chem. 2015, 87 (8), 4370-4376,
// DOI: 10.1021/acs.analchem.5b00140), written with Erik Marklund and Andrew Baldwin.
// Copyright (c) 2016, University of Oxford
//               2017-2025, Arizona Board of Regents on behalf of the University of Arizona
// All rights reserved. Distributed under the UniDec License (BSD with a citation requirement),
// whose full text with the conditions and the disclaimer is in udengine.h and in
// LICENSES/UniDec_LICENSE.txt of the package.
//
// What is different from unidec.exe, and why (the arithmetic of every step is the C source's,
// in the same order, in single precision; see the comments marked "port:"):
//  * no files: the input arrays come from the caller, the outputs go into Output (the caller
//    turns them into what it used to read from _mass.txt, _error.txt and the .bin files);
//  * no exit(), abort() or printf: every stop of the engine is an exception with a plain message
//    (the caller turns it into an error code), including "Setup is bad. No points are allowed"
//    after which unidec.exe wrote its outputs from a null array and Windows ended it (0xC0000409);
//  * no shared state: all arrays belong to the run; two runs in one process are independent;
//  * cancellation: a flag checked between iterations and inside the long loops;
//  * two data races of the C code made deterministic: MakeSparseBlur removes points from the
//    allowed set while other threads read it (here: evaluated in the order of the single threaded
//    engine, so the result no longer depends on the threads), and Reconvolve's maximum was
//    updated by all threads without synchronisation (here: the true maximum);
//  * powf(x, 2) is written x * x (what the compilers make of it anyway, exact);
//  * is_peak (peak detection) scans only the points within the window instead of the whole mass
//    axis for every candidate (the same points: the axis is sorted), arrays that the outputs no
//    longer need are freed before the mass grid is allocated, and temporary arrays of the
//    iterations are allocated once instead of every iteration;
//  * not ported (unidec.exe features MS Analysis does not use; refused with a message): ion
//    mobility, CD-MS, MetaUniDec/HDF5, mass list (mfile), manual assignment, DoubleDec (the only
//    user of FFTW in the one dimensional engine; no FFT library is linked).
#include "udengine.h"
#include "udmath.h"
#include "../common.h"
#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace ms {
namespace udengine {

namespace {

// ---------------------------------------------------------------- udtools.c
inline long idx2(int ncols, long r, int c) { return r * (long)ncols + c; }
inline long idx3(int ncols, int nrows, long r, int c, int d) { return r * (long)ncols * nrows + (long)c * nrows + d; }
inline int mod_(int a, int b) { const int r = a % b; return r < 0 ? r + b : r; }
inline int indexmod(int length, int r, int c) { return mod_(c - r, length); }
inline int fixk(int k, int lengthmz) {
    k = std::abs(k);
    if (k >= lengthmz) k = 2 * lengthmz - k - 2;
    return k;
}
inline float sq(float x) { return x * x; }   // port: powf(x, 2.f)
inline float clip0(float x) { return x > 0 ? x : 0; }   // clip(x, 0)
inline float calcmz(float mass, int z, float adductmass) { return (mass + (float)z * adductmass) / (float)z; }
inline float calcmass(float mz, int z, float adductmass) { return (float)z * (mz - adductmass); }
inline float Gaus(float x, float x0, float sig) { return udmath::expf(-sq(x - x0) / (2.f * sig * sig)); }
inline float nativecharge(float mass, float fudge) { return 0.0467f * (float)std::pow((double)mass, (double)0.533f) + fudge; }
inline float LinearInterpolate(float y1, float y2, float mu) { return y1 * (1 - mu) + y2 * mu; }
inline float LinearInterpolatePosition(float x1, float x2, float x) {
    if (x2 - x1 == 0) return 0;
    return (x - x1) / (x2 - x1);
}
inline float CubicInterpolate(float y0, float y1, float y2, float y3, float mu) {
    const float mu2 = mu * mu;
    const float a0 = y3 - y2 - y0 + y1;
    const float a1 = y0 - y1 - a0;
    const float a2 = y2 - y0;
    const float a3 = y1;
    return a0 * mu * mu2 + a1 * mu2 + a2 * mu + a3;
}

// mzpeakshape (sig 0 and an unknown shape are refused before the run instead of exit(103) / exit(14))
inline float mzpeakshape(float x, float y, float sig, int psfun) {
    if (psfun == 0) return udmath::expf(-sq(x - y) / (2.f * sig * sig));
    if (psfun == 1) return sq(sig / 2.f) / (sq(x - y) + sq(sig / 2.f));
    if (y < x) return udmath::expf(-sq(x - y) / (2.f * sig * sig * 0.180337f));
    return (sig / 2.f) * (sig / 2.f) / (sq(x - y) + sq(sig / 2));
}

int nearfast(const float* d, float point, int numdat) {
    int start = 0, length = numdat - 1;
    int diff = length - start;
    while (diff > 1) {
        const int mid = start + (length - start) / 2;
        if (point < d[mid]) length = mid;
        else if (point == d[mid]) return mid;
        else if (point > d[mid]) start = mid;
        diff = length - start;
    }
    return std::fabs(point - d[start]) >= std::fabs(point - d[length]) ? length : start;
}

// nearfast() for queries that do not decrease: the same index, found by moving on from the bracket of
// the previous query instead of a new binary search. d must be strictly increasing (the engine's m/z,
// SetupZtab refuses equal neighbours). nearfast returns the index of a point equal to the query, else
// of the nearer of the two points around it (of the upper one when |p - d[s]| >= |p - d[s + 1]|, the
// distances in float), with the pair (0, 1) below the data and (n - 2, n - 1) above.
struct NearSweep {
    const float* d = nullptr;
    int n = 0, s = -1;   // s: the largest index with d[s] <= p, kept within 0 .. n - 2 (-1: not started)
    void reset(const float* data, int count) { d = data; n = count; s = -1; }
    int operator()(float p) {
        if (n < 2) return 0;
        if (s < 0) {
            int lo = 0, hi = n - 1;   // d[lo] <= p < d[hi] when p is inside
            if (!(p >= d[0])) s = 0;
            else if (p >= d[n - 1]) s = n - 2;
            else {
                while (hi - lo > 1) {
                    const int mid = lo + (hi - lo) / 2;
                    if (d[mid] <= p) lo = mid;
                    else hi = mid;
                }
                s = lo;
            }
        }
        while (s < n - 2 && d[s + 1] <= p) s++;
        if (d[s] == p) return s;
        if (d[s + 1] == p) return s + 1;
        return std::fabs(p - d[s]) >= std::fabs(p - d[s + 1]) ? s + 1 : s;
    }
};

float vmax(const float* a, long n) {
    if (n <= 0) return 0;
    float m = a[0];
    for (long i = 0; i < n; i++) if (a[i] > m) m = a[i];
    return m;
}
float vsum(const float* a, long n) {
    float s = 0;
    for (long i = 0; i < n; i++) s += a[i];
    return s;
}
void norm_sum(float* a, int n) {
    const float norm = vsum(a, n);
    if (norm > 0) for (int i = 0; i < n; i++) a[i] = a[i] / norm;
}
// port: unidec.exe was built with Intel's compiler, which turns the plain sums of errfunspeedy,
// Average and CheckConvergence into SSE code with four interleaved partial sums (lane l adds the
// elements l, l + 4, ...; then (lane 0 + lane 2) + (lane 1 + lane 3); then the last n % 4 elements).
// The same order is used here: with the serial order of the C source the error (and so R squared
// and the UniScore) came out different from unidec.exe by up to 1e-4 (relative).
struct Sum4 {
    float acc[4] = {0, 0, 0, 0};
    void add(long i, float v) { acc[i & 3] += v; }
    // the total of elements 0 .. n - 1 added with add(i, ...) for i < n - n % 4, plus the tail
    template <class F> float total(long n, F&& tail) const {
        float s = (acc[0] + acc[2]) + (acc[1] + acc[3]);
        for (long i = n - n % 4; i < n; i++) s += tail(i);
        return s;
    }
};

float average(int n, const float* a) {
    Sum4 s;
    const long n4 = n - n % 4;
    for (long i = 0; i < n4; i++) s.add(i, a[i]);
    const float t1 = s.total(n, [&](long i) { return a[i]; });
    const float t2 = (float)n;
    if (t2 == 0) return 0;
    return t1 / t2;
}

// ---------------------------------------------------------------- the state of one run
struct Run {
    const Config& c;
    const Control& ctl;
    int L, Z;
    long ln;
    const float* mz;          // inp.dataMZ
    const float* in;          // inp.dataInt
    std::vector<int> nztab;
    std::vector<float> mtab;
    std::vector<char> ibarr;  // inp.barr
    // Decon
    std::vector<float> fitdat, baseline, noise, massgrid, massaxis, massaxisval, blur, newblur, mzdist, rmzdist;
    std::vector<int> starttab, endtab;
    float error = 0, rsquared = 0, uniscore = 0, conv = 0;
    int iterations = 0, mlen = 0, maxlength = 0;
    // IntraDecon
    int mlength = 0, zlength = 0, numclose = 0;
    float betafactor = 1.0f;
    int off = 0;
    std::vector<char> barr;
    std::vector<int> mind, zind, closemind, closezind, closeind;
    std::vector<float> mdist, zdist, oldblur, closeval, closearray, dataInt2;
    // scratch of the iterations (port: allocated once)
    std::vector<float> deltas, denom, scratch;
    // the engine's peaks (score)
    std::vector<float> peakx_out, peaky_out, dscores_out;

    Run(const Config& cfg, const Control& k, const float* x, const float* y, int n)
        : c(cfg), ctl(k), L(n), Z(cfg.numz), ln((long)n * cfg.numz), mz(x), in(y) {}

    void check() const { ctl.check(); }
};

// ---------------------------------------------------------------- inputs (SetupZtab, SetupInputs)
void setup_ztab(Run& r) {
    r.nztab.resize(r.Z);
    for (int i = 0; i < r.Z; i++) r.nztab[i] = i + r.c.startz;
    for (int j = 0; j < r.Z; j++)
        if (r.nztab[j] == 0)
            throw std::invalid_argument("The UniDec engine cannot use a charge of 0 (it is in the charge range " +
                                        std::to_string(r.c.startz) + " to " + std::to_string(r.c.endz) + ").");
    for (int i = 0; i < r.L - 1; i++)
        if (r.mz[i] == r.mz[i + 1]) {
            char b[200];
            snprintf(b, sizeof b, "The UniDec engine cannot use two data points with the same m/z (%f).", (double)r.mz[i]);
            throw std::invalid_argument(b);
        }
}

void setup_inputs(Run& r) {
    const Config& c = r.c;
    const int L = r.L, Z = r.Z;
    r.mtab.assign(r.ln, 0.0f);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < L; i++)
        for (int j = 0; j < Z; j++) r.mtab[idx2(Z, i, j)] = calcmass(r.mz[i], r.nztab[j], c.adductmass);
    // port: ignorezeros() is left out: SetLimits (TestMass) sets every element afterwards
    r.ibarr.assign(r.ln, 0);
    // port: nativecharge (a power function, slow) only where its value can change the outcome. For masses
    // between masslb and massub it lies between the two bounds below (with a margin for its rounding);
    // a charge that passes both tests at either bound with a margin of one passes them for every mass.
    double nl_lo = 0, nl_hi = 0;
    {
        const double lb = std::max((double)c.masslb, 0.0), ub = std::max((double)c.massub, 0.0);
        nl_lo = 0.0467 * std::pow(lb, (double)0.533f) * (1 - 1e-5);
        nl_hi = 0.0467 * std::pow(ub, (double)0.533f) * (1 + 1e-5);
    }
#pragma omp parallel for schedule(static)
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < Z; j++) {
            const float testmass = r.mtab[idx2(Z, i, j)];
            char ok = 0;
            if (testmass < c.massub && testmass > c.masslb) {
                const double z = (double)r.nztab[j];
                if (z < nl_lo + (double)c.nativezub - 1 && z > nl_hi + (double)c.nativezlb + 1) {
                    ok = 1;
                } else {
                    const float nativelimit = nativecharge(testmass, 0);
                    ok = (float)r.nztab[j] < nativelimit + c.nativezub && (float)r.nztab[j] > nativelimit + c.nativezlb;
                }
            }
            r.ibarr[idx2(Z, i, j)] = ok;
        }
    }
}

// ---------------------------------------------------------------- peak shape (SetStartsEnds, MakePeakShape1D/2D, SetUpPeakShape)
int set_starts_ends(Run& r) {
    const Config& c = r.c;
    int maxlength = 1;
    NearSweep ss, se;   // port: the points rise with i (see NearSweep)
    ss.reset(r.mz, r.L);
    se.reset(r.mz, r.L);
    for (int i = 0; i < r.L; i++) {
        const float window = c.psmzthresh;   // variablepw 0
        float point = r.mz[i] - window;
        int start, end;
        if (point < r.mz[0] && c.speedyflag == 0) start = 0 - nearfast(r.mz, 2 * r.mz[0] - point, r.L);
        else start = ss(point);
        r.starttab[i] = start;
        point = r.mz[i] + window;
        if (point > r.mz[r.L - 1] && c.speedyflag == 0) end = r.L - 1 + nearfast(r.mz, 2 * r.mz[0] - point, r.L);
        else end = se(point);
        r.endtab[i] = end;
        if (end - start > maxlength) maxlength = end - start;
    }
    return maxlength;
}

void make_peak_shape_2d(Run& r, int makereverse, int inflateflag) {
    const Config& c = r.c;
    float mzsig = std::fabs(c.mzsig);
    if (inflateflag == 1) mzsig *= c.peakshapeinflate;
    const int maxlength = r.maxlength;
#pragma omp parallel for schedule(static)
    for (int i = 0; i < r.L; i++) {
        const int start = r.starttab[i], end = r.endtab[i];
        for (int j = start; j <= end; j++) {
            const int j2 = fixk(j, r.L);
            r.mzdist[idx2(maxlength, i, j2 - start)] = mzpeakshape(r.mz[i], r.mz[j2], mzsig, c.psfun);
            if (makereverse == 1) r.rmzdist[idx2(maxlength, i, j2 - start)] = mzpeakshape(r.mz[j2], r.mz[i], mzsig, c.psfun);
        }
    }
}

void make_peak_shape_1d(Run& r, int makereverse, int inflateflag) {
    const Config& c = r.c;
    const float binsize = r.mz[1] - r.mz[0];
    const float newrange = c.psmzthresh / binsize;
    float mzsig = std::fabs(c.mzsig);
    if (inflateflag == 1) mzsig *= c.peakshapeinflate;
    if (!(std::fabs(newrange) < 2e9f))
        throw std::invalid_argument("The UniDec engine's peak shape is too wide for the bins of its input (peak width and bin size).");
    for (int n = (int)-newrange; n < (int)newrange; n++) {
        r.mzdist[indexmod(r.L, 0, n)] = mzpeakshape(0, (float)n * binsize, mzsig, c.psfun);
        if (makereverse == 1) r.rmzdist[indexmod(r.L, 0, n)] = mzpeakshape((float)n * binsize, 0, mzsig, c.psfun);
    }
}

void setup_peak_shape(Run& r) {
    const Config& c = r.c;
    r.starttab.assign(r.L, 0);
    r.endtab.assign(r.L, 0);
    if (c.mzsig != 0) {
        r.maxlength = set_starts_ends(r);
        long pslen = r.L;
        if (c.speedyflag == 0) pslen = (long)r.L * r.maxlength;
        if (pslen > INT_MAX)
            throw std::invalid_argument("The UniDec engine's peak shape table would exceed its 32 bit index (use Data reduction: Automatic or a smaller peak width).");
        r.mzdist.assign(pslen, 0.0f);
        int makereverse = 0;
        if (c.mzsig < 0 || c.beta < 0) {
            makereverse = 1;
            r.rmzdist.assign(pslen, 0.0f);
        }
        if (c.speedyflag == 0) make_peak_shape_2d(r, makereverse, 1);
        else make_peak_shape_1d(r, makereverse, 0);
    } else {
        r.maxlength = 0;
    }
}

// ---------------------------------------------------------------- neighbourhood blur (SetUpBlur, MakeSparseBlur)
// port: MakeSparseBlur ran in parallel and removed points (barr = 0) while other threads read barr, so
// which points survived near the edges depended on thread timing. Here the neighbours are found in
// parallel (that does not depend on barr) and the removal is applied in one pass in the order of the
// single threaded engine (i ascending, then j), reading barr as that engine would have.
void make_sparse_blur(Run& r) {
    const Config& c = r.c;
    const int L = r.L, Z = r.Z, numclose = r.numclose;
    const float molig = c.molig;
    char* barr = r.barr.data();
    int* closeind = r.closeind.data();
    float* closearray = r.closearray.data();
    // 1: candidate neighbour (before the barr test) and its value, for every allowed point. port: for one
    // charge index and neighbour the predicted m/z rise with i, so the nearest point is found with one
    // NearSweep per (j, k) and thread (the iterations of a thread are consecutive) instead of a binary search
#pragma omp parallel
    {
    std::vector<NearSweep> sweep((size_t)Z * numclose);
    for (NearSweep& w : sweep) w.reset(r.mz, L);
#pragma omp for schedule(static)
    for (int i = 0; i < L; i++) {
        if (r.ctl.stop()) continue;
        for (int j = 0; j < Z; j++) {
            if (barr[idx2(Z, i, j)] == 1) {
                float mzsig = c.mzsig;
                if (mzsig == 0) {
                    int i1 = i - 1, i2 = i + 1;
                    if (i >= L - 1) i2 = i;
                    if (i == 0) i1 = i;
                    mzsig = 2.0f * std::fabs(r.mz[i2] - r.mz[i1]);
                    if (mzsig > c.massbins || mzsig == 0) mzsig = c.massbins * 2;
                }
                const float newthreshold = mzsig * 2;
                for (int k = 0; k < numclose; k++) {
                    const long ck = idx3(Z, numclose, i, j, k);
                    const int indz = j + r.closezind[k];
                    if (indz < 0 || indz >= Z || (r.nztab[j] + r.closezind[k]) == 0) {
                        closeind[ck] = -1;
                        closearray[ck] = 0;
                        continue;
                    }
                    const float point = calcmz(r.mtab[idx2(Z, i, j)] + (float)r.closemind[k] * molig, r.nztab[j] + r.closezind[k], c.adductmass);
                    if (point < r.mz[0] - newthreshold || point > r.mz[L - 1] + newthreshold) {
                        closeind[ck] = -1;
                        closearray[ck] = 0;
                        continue;
                    }
                    // the predicted m/z rise with i only for a positive charge (and a molig term that does not change with i)
                    const int ind = (r.nztab[j] + r.closezind[k]) > 0 ? sweep[(size_t)j * numclose + k](point) : nearfast(r.mz, point, L);
                    const float closepoint = r.mz[ind];
                    if (std::fabs(point - closepoint) < newthreshold) {
                        closeind[ck] = (int)idx2(Z, ind, indz);
                        closearray[ck] = r.closeval[k] * mzpeakshape(point, closepoint, mzsig, c.psfun);
                    } else {
                        closeind[ck] = -1;
                        closearray[ck] = 0;
                    }
                }
            } else {
                for (int k = 0; k < numclose; k++) {
                    closeind[idx3(Z, numclose, i, j, k)] = -1;
                    closearray[idx3(Z, numclose, i, j, k)] = 0;
                }
            }
        }
    }
    }
    r.check();
    // 2: the barr test and the removal of points with fewer than two allowed neighbours, in order
    for (long cell = 0; cell < r.ln; cell++) {
        if (barr[cell] != 1) continue;
        int num = 0;
        int* ci = closeind + cell * numclose;
        float* ca = closearray + cell * numclose;
        for (int k = 0; k < numclose; k++) {
            if (ci[k] < 0) continue;
            if (barr[ci[k]] == 1) num += 1;
            else { ci[k] = -1; ca[k] = 0; }
        }
        if (num < 2 && c.manualflag == 0) barr[cell] = 0;
        if ((cell & 0xfffff) == 0) r.check();
    }
}

void setup_blur(Run& r) {
    const Config& c = r.c;
    if (c.zsig >= 0 && c.msig >= 0) {
        r.zlength = 1 + 2 * (int)c.zsig;
        r.mlength = 1 + 2 * (int)c.msig;
    } else {
        r.zlength = c.zsig != 0 ? 1 + 2 * (int)((double)(3 * std::fabs(c.zsig)) + 0.5) : 1;
        r.mlength = c.msig != 0 ? 1 + 2 * (int)((double)(3 * std::fabs(c.msig)) + 0.5) : 1;
    }
    if (r.zlength < 1 || r.mlength < 1 || r.zlength > 100000 || r.mlength > 100000)
        throw std::invalid_argument("The UniDec engine cannot use this charge or mass smoothing width.");
    r.numclose = r.mlength * r.zlength;
    r.mind.assign(r.mlength, 0);
    r.mdist.assign(r.mlength, 0.0f);
    for (int i = 0; i < r.mlength; i++) {
        r.mind[i] = i - (r.mlength - 1) / 2;
        r.mdist[i] = c.msig != 0 ? Gaus((float)i, ((float)r.mlength - 1) / 2.f, c.msig) : 1;
    }
    r.zind.assign(r.zlength, 0);
    r.zdist.assign(r.zlength, 0.0f);
    for (int i = 0; i < r.zlength; i++) {
        r.zind[i] = i - (r.zlength - 1) / 2;
        r.zdist[i] = c.zsig != 0 ? Gaus((float)i, ((float)r.zlength - 1) / 2.f, c.zsig) : 1;
    }
    r.closemind.assign(r.numclose, 0);
    r.closezind.assign(r.numclose, 0);
    r.closeval.assign(r.numclose, 0.0f);
    if ((double)r.numclose * (double)r.ln > (double)INT_MAX)
        throw std::invalid_argument("The UniDec engine's neighbour table (" + std::to_string((long long)r.numclose * r.ln) +
                                    " entries) would exceed its 32 bit index. Use Data reduction: Automatic, a narrower m/z or charge range, or smaller smoothing widths.");
    r.closeind.assign((size_t)r.numclose * r.ln, 0);
    r.closearray.assign((size_t)r.numclose * r.ln, 0.0f);
    for (int k = 0; k < r.numclose; k++) {
        r.closemind[k] = r.mind[k % r.mlength];
        r.closezind[k] = r.zind[k / r.mlength];
        r.closeval[k] = r.zdist[k / r.mlength] * r.mdist[k % r.mlength];
    }
    norm_sum(r.mdist.data(), r.mlength);
    norm_sum(r.zdist.data(), r.zlength);
    norm_sum(r.closeval.data(), r.numclose);
    bool any = false;
    for (long i = 0; i < r.ln && !any; i++) any = r.barr[i] == 1;
    if (!any)   // port: exit(10)
        throw std::runtime_error("The UniDec engine found no data point that gives a mass inside the mass range at the chosen charges. Check the m/z, mass and charge ranges.");
    make_sparse_blur(r);
}

// ---------------------------------------------------------------- baseline (aggressiveflag 1 and 2; MS Analysis does not use it)
void midblur_baseline(Run& r, std::vector<float>& base, int mult) {
    const int L = r.L;
    if (mult == 0) mult = L / 400;
    std::vector<float> temp(base);
#pragma omp parallel
    {
        const int window = 25, len = window * 2;
        std::vector<float> med(len);
#pragma omp for schedule(static)
        for (int i = 0; i < L; i++) {
            int index = 0;
            for (int j = -window; j < window; j++) {
                const int k = i + j * mult;
                float newval = 0;
                if (k >= 0 && k < L) newval = temp[k];
                if (k < 0 && -k < L) newval = temp[-k];   // port: bounds checked (the C code read outside)
                if (k >= L && 2 * L - k >= 0 && 2 * L - k < L) newval = temp[2 * L - k];
                med[index++] = newval;
            }
            std::sort(med.begin(), med.end());
            float val = 0;
            index = 0;
            for (int j = 0; j < window; j++) { val += med[j]; index++; }
            base[i] = index != 0 ? val / (float)index : 0;
        }
    }
}

void blur_baseline(Run& r, std::vector<float>& base, float mzsig, int mult, int filterwidth) {
    const int L = r.L;
    const int mulin = mult;
    std::vector<float> temp(base);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < L; i++) {
        // port: mult is private to each point (the C code shared it between the threads)
        int m = mult;
        const float mzdiff = i > 0 ? r.mz[i] - r.mz[i - 1] : r.mz[i + 1] - r.mz[i];
        if (mulin == 0 && mzdiff > 0) m = (int)(2 * mzsig / mzdiff);
        if (m < 1) m = 1;
        float val = 0;
        const int window = filterwidth;
        for (int j = -window; j < window; j++) {
            const int k = i + j * m;
            float newval = 0;
            if (k >= 0 && k < L) newval = temp[k];
            if (k < 0 && -k < L) newval = temp[-k];
            if (k >= L && 2 * L - k >= 0 && 2 * L - k < L) newval = temp[2 * L - k];
            val += newval;
        }
        base[i] = val / ((float)window * 2 + 1);
    }
}

void deconvolve_baseline(Run& r, std::vector<float>& base, float mzsig) {
    midblur_baseline(r, base, 0);
    midblur_baseline(r, base, 5);
    std::vector<float> dn(base);
    for (int i = 0; i < r.L; i++)
        if (dn[i] != 0 && r.in[i] >= 0) dn[i] = r.in[i] / dn[i];
    midblur_baseline(r, dn, 0);
    midblur_baseline(r, dn, 5);
    for (int i = 0; i < r.L; i++) base[i] = base[i] * dn[i];
    (void)mzsig;
}

// ---------------------------------------------------------------- SetupDeconvolution
// returns true when no point is left (the "badness" of the C code)
bool setup_deconvolution(Run& r) {
    const Config& c = r.c;
    const int L = r.L, Z = r.Z;
    r.barr = r.ibarr;
    std::vector<char>().swap(r.ibarr);   // port: inp.barr is not used again
    setup_peak_shape(r);
    r.check();
    setup_blur(r);
    r.check();
    if (c.intthresh != -1) {   // KillB
        for (int i = 0; i < L; i++)
            if (r.in[i] <= c.intthresh)
                for (int j = 0; j < Z; j++) r.barr[idx2(Z, i, j)] = 0;
    }
    r.blur.assign(r.ln, 0.0f);
    r.newblur.assign(r.ln, 0.0f);
    r.oldblur.assign(r.ln, 0.0f);
    if (c.baselineflag == 1) {
        r.baseline.assign(L, 0.0f);
        r.noise.assign(L, 0.0f);
    }
    for (int i = 0; i < L; i++) {
        const float val = r.in[i] / ((float)(Z + 2));
        if (c.baselineflag == 1) {
            r.baseline[i] = val;
            r.noise[i] = val;
        }
        for (int j = 0; j < Z; j++) r.blur[idx2(Z, i, j)] = r.barr[idx2(Z, i, j)] == 1 ? val : 0;
    }
    r.oldblur = r.blur;
    r.newblur = r.blur;
    r.dataInt2.assign(r.in, r.in + L);
    if (c.baselineflag == 1 && c.mzsig != 0) {
        deconvolve_baseline(r, r.baseline, std::fabs(c.mzsig));
        if (c.aggressiveflag == 2) {
            for (int i = 0; i < 10; i++) deconvolve_baseline(r, r.baseline, std::fabs(c.mzsig));
            for (int i = 0; i < L; i++)
                if (r.baseline[i] > 0) r.dataInt2[i] -= r.baseline[i];
        }
    }
    const float dmax = vmax(r.dataInt2.data(), L);
    if (dmax > 1) r.betafactor = dmax;
    for (long i = 0; i < r.ln; i++)
        if (r.barr[i] == 1) return false;
    return true;
}

// ---------------------------------------------------------------- blurs
void blur_it(Run& r) {
    if (r.numclose == 1) { r.newblur = r.blur; return; }
    const int nc = r.numclose;
    const float* blur = r.blur.data();
    float* nb = r.newblur.data();
#pragma omp parallel for schedule(static)
    for (long i = 0; i < r.ln; i++) {
        if ((i & 0xffff) == 0 && r.ctl.stop()) continue;
        float temp = 0;
        if (r.barr[i] == 1) {
            for (int k = 0; k < nc; k++) {
                const int ci = r.closeind[i * nc + k];
                if (ci != -1) temp += r.closearray[i * nc + k] * blur[ci];
            }
        }
        nb[i] = temp;
    }
}

void blur_it_mean(Run& r) {
    if (r.numclose == 1) { r.newblur = r.blur; return; }
    const int nc = r.numclose;
    const float zerolog = r.c.zerolog;
    const float* blur = r.blur.data();
    const char* barr = r.barr.data();
    const int* closeind = r.closeind.data();
    const float* closearray = r.closearray.data();
    float* nb = r.newblur.data();
    // port: the logarithms and exponentials of a block of points are taken together (udmath log_batch and
    // exp_batch: the values of udmath::logf / expf, two at a time); the sum of each point is formed in the
    // order of the C code
    const long B = 512;
#pragma omp parallel
    {
        std::vector<float> arg((size_t)B * nc), lg((size_t)B * nc), earg(B), ex(B);
        std::vector<char> pos((size_t)B * nc);
        std::vector<long> act(B);
#pragma omp for schedule(static)
        for (long b0 = 0; b0 < r.ln; b0 += B) {
            if (r.ctl.stop()) continue;
            const long b1 = std::min(r.ln, b0 + B);
            long na = 0;
            for (long i = b0; i < b1; i++) {
                if (barr[i] == 1) act[na++] = i;
                else nb[i] = 0;
            }
            for (long a = 0; a < na; a++) {
                const long i = act[a];
                for (int k = 0; k < nc; k++) {
                    float temp2 = 0;
                    const int ci = closeind[i * nc + k];
                    if (ci != -1) temp2 = blur[ci] * closearray[i * nc + k];
                    pos[a * nc + k] = temp2 > 0;
                    arg[a * nc + k] = temp2 > 0 ? temp2 : 1.0f;
                }
            }
            udmath::log_batch(arg.data(), lg.data(), na * nc);
            for (long a = 0; a < na; a++) {
                float temp = 0;
                for (int k = 0; k < nc; k++) temp += pos[a * nc + k] ? lg[a * nc + k] : zerolog;
                earg[a] = temp / (float)nc;
            }
            udmath::exp_batch(earg.data(), ex.data(), na);
            for (long a = 0; a < na; a++) nb[act[a]] = ex[a];
        }
    }
}

void blur_it_hybrid1(Run& r) {
    if (r.numclose == 1) { r.newblur = r.blur; return; }
    const int L = r.L, Z = r.Z, nc = r.numclose;
    const float zerolog = r.c.zerolog;
#pragma omp parallel for schedule(static)
    for (int i = 0; i < L; i++) {
        if (r.ctl.stop()) continue;
        for (int j = 0; j < Z; j++) {
            float temp = 0;
            if (r.barr[idx2(Z, i, j)] == 1) {
                for (int n = 0; n < r.mlength; n++) {
                    float temp2 = 0;
                    for (int k = 0; k < r.zlength; k++) {
                        const int m = (int)idx2(r.mlength, k, n);
                        float temp3 = 0;
                        const long q = idx3(Z, nc, i, j, m);
                        if (r.closeind[q] != -1) temp3 = r.blur[r.closeind[q]] * r.closearray[q];
                        if (temp3 > 0) temp2 += udmath::logf(temp3);
                        else temp2 += zerolog;
                    }
                    temp += udmath::expf(temp2 / (float)r.zlength) * r.mdist[n];
                }
            }
            r.newblur[idx2(Z, i, j)] = temp;
        }
    }
}

void blur_it_hybrid2(Run& r) {
    if (r.numclose == 1) { r.newblur = r.blur; return; }
    const int L = r.L, Z = r.Z, nc = r.numclose;
    const float zerolog = r.c.zerolog;
#pragma omp parallel for schedule(static)
    for (int i = 0; i < L; i++) {
        if (r.ctl.stop()) continue;
        for (int j = 0; j < Z; j++) {
            float temp = 0;
            if (r.barr[idx2(Z, i, j)] == 1) {
                for (int n = 0; n < r.mlength; n++) {
                    float temp2 = 0;
                    for (int k = 0; k < r.zlength; k++) {
                        const int m = (int)idx2(r.mlength, k, n);
                        const long q = idx3(Z, nc, i, j, m);
                        if (r.closeind[q] != -1) temp2 += r.blur[r.closeind[q]] * r.zdist[k] * r.closearray[q];
                    }
                    if (temp2 > 0) temp += udmath::logf(temp2);
                    else temp += zerolog;
                }
                temp = udmath::expf(temp / (float)r.mlength);
            }
            r.newblur[idx2(Z, i, j)] = temp;
        }
    }
}

// ---------------------------------------------------------------- iteration
void convolve_simp(Run& r, const float* mzdist, const float* deltas, float* denom) {
    const int L = r.L, maxlength = r.maxlength;
    const int* starttab = r.starttab.data();
    const int* endtab = r.endtab.data();
    if (r.c.speedyflag == 0) {
#pragma omp parallel for schedule(static)
        for (int i = 0; i < L; i++) {
            float cv = 0;
            for (int k = starttab[i]; k <= endtab[i]; k++) {
                const int k2 = fixk(k, L);
                const int start = starttab[k2];
                cv += deltas[k2] * mzdist[idx2(maxlength, k2, i - start)];
            }
            denom[i] = cv;
        }
    } else {
        // port: indexmod(L, k, i) = (i - k) mod L without the division: i - k lies in (-L, L) (start and
        // end are indexes of the data), so the k below and above i are two runs (same order of the sum)
#pragma omp parallel for schedule(static)
        for (int i = 0; i < L; i++) {
            float cv = 0;
            const int s0 = starttab[i], s1 = endtab[i];
            const int mid = std::min(s1, i);
            int k = s0;
            for (; k <= mid; k++) cv += deltas[k] * mzdist[i - k];
            for (; k <= s1; k++) cv += deltas[k] * mzdist[i - k + L];
            denom[i] = cv;
        }
    }
}

void sum_deltas(Run& r, const float* blur, float* deltas) {
    const int L = r.L, Z = r.Z;
    const char* barr = r.barr.data();
#pragma omp parallel for schedule(static)
    for (int i = 0; i < L; i++) {
        float temp = 0;
        for (int j = 0; j < Z; j++)
            if (barr[idx2(Z, i, j)] == 1) temp += blur[idx2(Z, i, j)];
        deltas[i] = temp;
    }
}

void apply_ratios(Run& r, const float* blur, const float* dn, float* blur2) {
    const int L = r.L, Z = r.Z;
    const char* barr = r.barr.data();
#pragma omp parallel for schedule(static)
    for (int i = 0; i < L; i++)
        for (int j = 0; j < Z; j++) {
            const long q = idx2(Z, i, j);
            blur2[q] = barr[q] == 1 ? dn[i] * blur[q] : 0;
        }
}

void deconvolve_iteration_speedy(Run& r) {
    const Config& c = r.c;
    const int L = r.L;
    float* deltas = r.deltas.data();
    float* denom = r.denom.data();
    if (c.aggressiveflag == 1 && c.mzsig != 0) blur_baseline(r, r.baseline, std::fabs(c.mzsig), 0, c.filterwidth);
    sum_deltas(r, r.newblur.data(), deltas);
    if (c.mzsig != 0 && c.psig >= 0) convolve_simp(r, r.mzdist.data(), deltas, denom);
    else std::memcpy(denom, deltas, sizeof(float) * L);
    if (c.aggressiveflag == 1)
        for (int i = 0; i < L; i++) denom[i] += r.baseline[i];
#pragma omp parallel for schedule(static)
    for (int i = 0; i < L; i++)
        if (denom[i] != 0 && r.dataInt2[i] >= 0) denom[i] = r.dataInt2[i] / denom[i];
    if (c.mzsig < 0) {
        convolve_simp(r, r.rmzdist.data(), denom, deltas);
        std::memcpy(denom, deltas, sizeof(float) * L);
    }
    apply_ratios(r, r.newblur.data(), denom, r.blur.data());
    if (c.aggressiveflag == 1) {
        std::vector<float> d(denom, denom + L);
        blur_baseline(r, d, std::fabs(c.mzsig), 0, c.filterwidth);
        for (int i = 0; i < L; i++) r.baseline[i] = r.baseline[i] * d[i];
    }
}

void softargmax(Run& r, float beta) {
    const int L = r.L, Z = r.Z;
    r.scratch = r.blur;
    const float* nb = r.scratch.data();
    float* blur = r.blur.data();
#pragma omp parallel for schedule(static)
    for (int i = 0; i < L; i++) {
        float sum2 = 0, sum1 = 0, factor = 0, min2 = 1.0f;
        for (int j = 0; j < Z; j++) {
            const float d = nb[idx2(Z, i, j)];
            sum1 += d;
            const float e = udmath::expf(beta * d);
            if (e < min2) min2 = e;
            blur[idx2(Z, i, j)] = e;
            sum2 += e;
        }
        const float dn = sum2 - min2 * (float)Z;
        if (dn != 0) factor = sum1 / dn;
        if (factor > 0) {
            for (int j = 0; j < Z; j++) {
                blur[idx2(Z, i, j)] -= min2;
                blur[idx2(Z, i, j)] *= factor;
            }
        } else {
            for (int j = 0; j < Z; j++) blur[idx2(Z, i, j)] = 0;
        }
    }
}

void point_smoothing(Run& r, int width) {
    const int L = r.L, Z = r.Z;
    const float fwidth = (float)width;
    r.scratch = r.blur;
    const float* nb = r.scratch.data();
    float* blur = r.blur.data();
    const char* barr = r.barr.data();
#pragma omp parallel for schedule(static)
    for (int i = 0; i < L; i++) {
        for (int j = 0; j < Z; j++) {
            if (barr[idx2(Z, i, j)] == 1) {
                int low = i - width;
                if (low < 0) low = 0;
                int high = i + width + 1;
                if (high > L) high = L;
                float sum = 0;
                for (int k = low; k < high; k++) sum += nb[idx2(Z, k, j)];
                blur[idx2(Z, i, j)] = sum / (1.0f + 2.0f * fwidth);
            }
        }
    }
}

void suppression_harmonic(Run& r) {
    const int L = r.L, Z = r.Z;
    r.scratch = r.blur;
    const float* input = r.scratch.data();
    float* blur = r.blur.data();
    const int zmax = r.nztab[Z - 1], zmin = r.nztab[0];
    const int end = (zmax / 2) - zmin;
    for (int i = 0; i < L; i++) {   // port: serial (the C code wrote cells of other rows' j2 in parallel; same result)
        for (int j = 0; j < end; j++) {
            const int j2 = 2 * j + zmin;
            if (j2 < 0 || j2 >= Z) continue;   // port: the C code read outside the row here
            const float val = input[idx2(Z, i, j)], val2 = input[idx2(Z, i, j2)];
            if (val < val2) blur[idx2(Z, i, j)] = 0;
            if (val > val2) blur[idx2(Z, i, j2)] = 0;
        }
    }
}

void suppression_satellite(Run& r, int n) {
    const int L = r.L, Z = r.Z;
    if (n <= 0 || n >= Z / 2) return;
    r.scratch = r.blur;
    const float* input = r.scratch.data();
    float* blur = r.blur.data();
#pragma omp parallel for schedule(static)
    for (int i = 0; i < L; i++) {
        for (int j = n; j < Z - n - 1; j++) {
            float maxval = 0;
            for (int k = 0; k < 2 * n + 1; k++) {
                const int index = j - n + k;
                if (index < 0 || index >= Z) continue;
                if (input[idx2(Z, i, index)] > maxval) maxval = input[idx2(Z, i, index)];
            }
            for (int k = 0; k < 2 * n + 1; k++) {
                const int index = j - n + k;
                if (index < 0 || index >= Z) continue;
                if (input[idx2(Z, i, index)] < maxval) blur[idx2(Z, i, index)] = 0;
            }
        }
    }
}

void highest_n_chargestates(Run& r, int n, float zcutpercent) {
    const int L = r.L, Z = r.Z;
    if (n <= 0 || n >= Z) return;
    float* blur = r.blur.data();
#pragma omp parallel
    {
        std::vector<float> tmp(Z);
#pragma omp for schedule(static)
        for (int i = 0; i < L; i++) {
            for (int j = 0; j < Z; j++) tmp[j] = blur[idx2(Z, i, j)];
            for (int a = 0; a < n; a++) {
                int maxidx = a;
                for (int b = a + 1; b < Z; b++) if (tmp[b] > tmp[maxidx]) maxidx = b;
                const float t = tmp[a]; tmp[a] = tmp[maxidx]; tmp[maxidx] = t;
            }
            const float threshold = tmp[n - 1];
            for (int j = 0; j < Z; j++)
                if (blur[idx2(Z, i, j)] < threshold) blur[idx2(Z, i, j)] = blur[idx2(Z, i, j)] * zcutpercent;
        }
    }
}

void clip_minor_chargestates(Run& r, float zcutoff, float zcutpercent) {
    const int L = r.L, Z = r.Z;
    if (zcutoff <= 0.0f || zcutoff >= 1.0f) return;
    float* blur = r.blur.data();
#pragma omp parallel for schedule(static)
    for (int i = 0; i < L; i++) {
        float maxval = 0.0f;
        for (int j = 0; j < Z; j++) maxval = std::max(maxval, blur[idx2(Z, i, j)]);
        const float threshold = zcutoff * maxval;
        for (int j = 0; j < Z; j++)
            if (blur[idx2(Z, i, j)] < threshold) blur[idx2(Z, i, j)] = blur[idx2(Z, i, j)] * zcutpercent;
    }
}

void run_iteration(Run& r, int iterations) {
    const Config& c = r.c;
    r.iterations = iterations;
    if (c.beta > 0 && iterations > 0) softargmax(r, c.beta / r.betafactor);
    if (c.psig >= 1 && iterations > 0) point_smoothing(r, std::abs((int)c.psig));
    if (c.suppression_satellite > 0 && iterations > c.suppression_startit) suppression_satellite(r, c.suppression_satellite);
    if (c.suppression_harmonic > 0 && iterations > c.suppression_startit) suppression_harmonic(r);
    if (c.suppression_topn > 0 && iterations > c.suppression_startit) highest_n_chargestates(r, c.suppression_topn, c.suppression_percent);
    if (c.suppression_topx > 0 && iterations > c.suppression_startit) clip_minor_chargestates(r, c.suppression_topx, c.suppression_percent);
    if (c.zsig >= 0 && c.msig >= 0) blur_it_mean(r);
    else if (c.zsig > 0 && c.msig < 0) blur_it_hybrid1(r);
    else if (c.zsig < 0 && c.msig > 0) blur_it_hybrid2(r);
    else blur_it(r);
    r.check();
    deconvolve_iteration_speedy(r);
}

int check_convergence(Run& r, int iterations) {
    const Config& c = r.c;
    const float* blur = r.blur.data();
    const float* old = r.oldblur.data();
    const char* barr = r.barr.data();
    Sum4 sd, st;   // port: the order of unidec.exe (see Sum4)
    const long n4 = r.ln - r.ln % 4;
    for (long i = 0; i < n4; i++) {
        if (barr[i] == 1) {
            sd.add(i, sq(blur[i] - old[i]));
            st.add(i, blur[i]);
        }
    }
    const float diff = sd.total(r.ln, [&](long i) { return barr[i] == 1 ? sq(blur[i] - old[i]) : 0.0f; });
    const float tot = st.total(r.ln, [&](long i) { return barr[i] == 1 ? blur[i] : 0.0f; });
    if (tot != 0) r.conv = diff / tot;
    else if (r.conv == 12345678) return 1;
    if (r.conv < 0.000001) {
        if (r.off == 1 && c.numit > 0) return 1;
        r.off = 1;
    }
    (void)iterations;
    std::memcpy(r.oldblur.data(), r.blur.data(), (size_t)r.ln * sizeof(float));
    return 0;
}

// ---------------------------------------------------------------- outputs
void getfitdatspeedy(Run& r, float* fitdat, const float* blur, float maxint) {
    const int L = r.L, Z = r.Z;
    std::vector<float> deltas(L);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < L; i++) {
        float temp = 0;
        for (int j = 0; j < Z; j++) temp += blur[idx2(Z, i, j)];
        deltas[i] = temp;
    }
    if (r.maxlength != 0) convolve_simp(r, r.mzdist.data(), deltas.data(), fitdat);
    else std::memcpy(fitdat, deltas.data(), sizeof(float) * L);
    float fitmax = 0;
    for (int i = 0; i < L; i++) if (fitdat[i] > fitmax) fitmax = fitdat[i];
    if (fitmax != 0)
        for (int i = 0; i < L; i++) fitdat[i] = fitdat[i] < 0 ? 0 : fitdat[i] * maxint / fitmax;
}

float errfunspeedy(Run& r) {
    const Config& c = r.c;
    const int L = r.L;
    float maxint = 0;
    for (int i = 0; i < L; i++) if (r.in[i] > maxint) maxint = r.in[i];
    getfitdatspeedy(r, r.fitdat.data(), r.blur.data(), maxint);
    if (c.baselineflag == 1)
        for (int i = 0; i < L; i++) r.fitdat[i] += r.baseline[i];
    for (int i = 0; i < L; i++) if (r.fitdat[i] < 0) r.fitdat[i] = 0;   // ApplyCutoff(fitdat, 0)
    const float fitmean = average(L, r.in);
    Sum4 se, ss;   // port: the order of unidec.exe (see Sum4)
    const long n4 = L - L % 4;
    for (long i = 0; i < n4; i++) {
        se.add(i, sq(r.fitdat[i] - r.in[i]));
        ss.add(i, sq(r.in[i] - fitmean));
    }
    const float error = se.total(L, [&](long i) { return sq(r.fitdat[i] - r.in[i]); });
    const float sstot = ss.total(L, [&](long i) { return sq(r.in[i] - fitmean); });
    if (sstot != 0) r.rsquared = 1 - (error / sstot);
    return error;
}

// port: the maximum is the true one (the C code updated one shared variable from all threads)
float reconvolve(Run& r) {
    const Config& c = r.c;
    const int L = r.L, Z = r.Z, maxlength = r.maxlength;
    const char* barr = r.barr.data();
    const float* blur = r.blur.data();
    float* nb = r.newblur.data();
    float newblurmax = 0;
#pragma omp parallel
    {
        float tmax = 0;
#pragma omp for schedule(static)
        for (int i = 0; i < L; i++) {
            if (r.ctl.stop()) continue;
            for (int j = 0; j < Z; j++) {
                float cv = 0;
                if (barr[idx2(Z, i, j)] == 1) {
                    for (int k = r.starttab[i]; k <= r.endtab[i]; k++) {
                        if (c.speedyflag == 0) {
                            const int k2 = fixk(k, L);
                            if (blur[idx2(Z, k2, j)] != 0) {
                                const int start = r.starttab[k2];
                                cv += blur[idx2(Z, k2, j)] * r.mzdist[idx2(maxlength, k2, i - start)];
                            }
                        } else if (blur[idx2(Z, k, j)] != 0) {
                            cv += blur[idx2(Z, k, j)] * r.mzdist[k <= i ? i - k : i - k + L];   // indexmod(L, k, i)
                        }
                    }
                }
                nb[idx2(Z, i, j)] = cv;
                if (cv > tmax) tmax = cv;
            }
        }
#pragma omp critical(udengine_reconvolve)
        if (tmax > newblurmax) newblurmax = tmax;
    }
    r.check();
    return newblurmax;
}

void integrate_transform(Run& r, const float* blur, float massmax, float massmin) {
    const int L = r.L, Z = r.Z;
    for (int i = 0; i < L; i++) {
        if ((i & 0xfff) == 0) r.check();
        for (int j = 0; j < Z; j++) {
            const float testmass = r.mtab[idx2(Z, i, j)];
            if (testmass < massmax && testmass > massmin) {
                const int index = nearfast(r.massaxis.data(), testmass, r.mlen);
                const float newval = blur[idx2(Z, i, j)];
                if (r.massaxis[index] == testmass) {
                    r.massaxisval[index] += newval;
                    r.massgrid[idx2(Z, index, j)] += newval;
                }
                if (r.massaxis[index] < testmass && index < r.mlen - 2) {
                    const int index2 = index + 1;
                    const float interpos = LinearInterpolatePosition(r.massaxis[index], r.massaxis[index2], testmass);
                    r.massaxisval[index] += (1.0f - interpos) * newval;
                    r.massgrid[idx2(Z, index, j)] += (1.0f - interpos) * newval;
                    r.massaxisval[index2] += interpos * newval;
                    r.massgrid[idx2(Z, index2, j)] += interpos * newval;
                }
                if (r.massaxis[index] > testmass && index > 0) {
                    const int index2 = index - 1;
                    const float interpos = LinearInterpolatePosition(r.massaxis[index], r.massaxis[index2], testmass);
                    r.massaxisval[index] += (1 - interpos) * newval;
                    r.massgrid[idx2(Z, index, j)] += (1 - interpos) * newval;
                    r.massaxisval[index2] += interpos * newval;
                    r.massgrid[idx2(Z, index2, j)] += interpos * newval;
                }
            }
        }
    }
}

void interpolate_transform(Run& r, const float* blur) {
    const int L = r.L, Z = r.Z;
    const float startmzval = r.mz[0], endmzval = r.mz[L - 1];
#pragma omp parallel for schedule(static)
    for (int i = 0; i < r.mlen; i++) {
        if ((i & 0xfff) == 0 && r.ctl.stop()) continue;
        float val = 0;
        for (int j = 0; j < Z; j++) {
            const float mztest = calcmz(r.massaxis[i], r.nztab[j], r.c.adductmass);
            if (mztest > startmzval && mztest < endmzval) {
                int index = nearfast(r.mz, mztest, L);
                int index2 = index;
                float newval = 0;
                if (r.mz[index] == mztest) {
                    newval = blur[idx2(Z, index, j)];
                    val += newval;
                    r.massgrid[idx2(Z, i, j)] = newval;
                } else {
                    if (r.mz[index] > mztest && index > 1 && index < L - 1) {
                        index2 = index;
                        index = index - 1;
                    } else if (r.mz[index] < mztest && index < L - 2 && index > 0) {
                        index2 = index + 1;
                    }
                    if (index2 > index && (r.mz[index2] - r.mz[index]) != 0) {
                        const float mu = (mztest - r.mz[index]) / (r.mz[index2] - r.mz[index]);
                        const float y0 = blur[idx2(Z, index - 1, j)];
                        const float y1 = blur[idx2(Z, index, j)];
                        const float y2 = blur[idx2(Z, index2, j)];
                        const float y3 = blur[idx2(Z, index2 + 1, j)];
                        newval = clip0(CubicInterpolate(y0, y1, y2, y3, mu));
                        val += newval;
                        r.massgrid[idx2(Z, i, j)] = newval;
                    }
                }
            }
        }
        r.massaxisval[i] = val;
    }
    r.check();
}

void smart_transform(Run& r, const float* blur) {
    const Config& c = r.c;
    const int L = r.L, Z = r.Z;
    const float startmzval = r.mz[0], endmzval = r.mz[L - 1];
    const float* massaxis = r.massaxis.data();
    const float* mz = r.mz;
    // port: for one charge the m/z of the mass axis rise with i (for z > 0): NearSweep per charge and thread
#pragma omp parallel
    {
    std::vector<NearSweep> s0(Z), s1(Z), s2(Z);
    for (int j = 0; j < Z; j++) { s0[j].reset(mz, L); s1[j].reset(mz, L); s2[j].reset(mz, L); }
#pragma omp for schedule(static)
    for (int i = 0; i < r.mlen; i++) {
        if ((i & 0xfff) == 0 && r.ctl.stop()) continue;
        float val = 0;
        for (int j = 0; j < Z; j++) {
            const int z = r.nztab[j];
            const float mtest = massaxis[i];
            auto cmz = [&](float m) { return calcmz(m, z, c.adductmass); };
            const float mztest = cmz(mtest);
            float mzlower, mlower;
            if (i > 0) {
                mlower = massaxis[i - 1];
                mzlower = cmz(mlower);
            } else {
                mzlower = mztest;
                mlower = mtest;
            }
            float mzupper, mupper;
            if (i < r.mlen - 1) {
                mupper = massaxis[i + 1];
                mzupper = cmz(mupper);
            } else {
                mzupper = mztest;
                mupper = mtest;
            }
            if (mzupper > startmzval && mzlower < endmzval) {
                int index = z > 0 ? s0[j](mztest) : nearfast(mz, mztest, L);
                const int index1 = z > 0 ? s1[j](mzlower) : nearfast(mz, mzlower, L);
                int index2 = z > 0 ? s2[j](mzupper) : nearfast(mz, mzupper, L);
                const float imz = mz[index];
                float newval = 0;
                if (index2 - index1 < 5) {
                    if (imz == mztest) {
                        newval = clip0(blur[idx2(Z, index, j)]);
                        val += newval;
                        r.massgrid[idx2(Z, i, j)] = newval;
                    } else {
                        int edge = 0;
                        index2 = index;
                        if (imz > mztest) index = index - 1;
                        else if (imz < mztest) index2 = index + 1;
                        if (index < 1 || index2 >= L - 1) edge = 1;
                        if (index < 0 || index2 >= L) edge = 2;
                        // port: edge tested first (the C code read mz[-1] or mz[L] here before finding edge != 0)
                        if (edge == 0 && index2 > index && (mz[index2] - mz[index]) != 0) {
                            const float mu = (mztest - mz[index]) / (mz[index2] - mz[index]);
                            const float y0 = blur[idx2(Z, index - 1, j)];
                            const float y1 = blur[idx2(Z, index, j)];
                            const float y2 = blur[idx2(Z, index2, j)];
                            const float y3 = blur[idx2(Z, index2 + 1, j)];
                            newval = clip0(CubicInterpolate(y0, y1, y2, y3, mu));
                            val += newval;
                            r.massgrid[idx2(Z, i, j)] = newval;
                        } else if (edge == 1 && (mz[index2] - mz[index]) != 0) {
                            const float mu = (mztest - mz[index]) / (mz[index2] - mz[index]);
                            const float y1 = blur[idx2(Z, index, j)];
                            const float y2 = blur[idx2(Z, index2, j)];
                            newval = clip0(LinearInterpolate(y1, y2, mu));
                            val += newval;
                            r.massgrid[idx2(Z, i, j)] = newval;
                        } else if (edge == 2) {
                            // port: the C code read mz[-1] or mz[L] in this test before fixing the indexes
                            if (index2 == 0) { index = 0; index2 = 1; }
                            if (index == L - 1) { index = L - 1; index2 = L - 2; }
                            if (index >= 0 && index < L && index2 >= 0 && index2 < L && (mz[index2] - mz[index]) != 0) {
                                const float mu = (mztest - mz[index]) / (mz[index] - mz[index2]);
                                const float y1 = blur[idx2(Z, index, j)];
                                newval = clip0(LinearInterpolate(y1, 0, mu));
                                val += newval;
                                r.massgrid[idx2(Z, i, j)] = newval;
                            }
                        }
                    }
                } else {
                    float num = 0;
                    for (int k = index1; k < index2 + 1; k++) {
                        const float kmz = mz[k];
                        const float km = (kmz - c.adductmass) * (float)z;
                        float scale;
                        if (mztest < kmz && km < mupper) scale = LinearInterpolatePosition(mupper, mtest, km);
                        else if (kmz < mztest && km > mlower) scale = LinearInterpolatePosition(mlower, mtest, km);
                        else if (kmz == mztest) scale = 1;
                        else scale = 0;
                        newval += scale * blur[idx2(Z, k, j)];
                        num += scale;
                    }
                    if (num != 0) newval /= num;
                    newval = clip0(newval);
                    val += newval;
                    r.massgrid[idx2(Z, i, j)] = newval;
                }
            }
        }
        r.massaxisval[i] = val;
    }
    }
    r.check();
}

}  // namespace

// in udscore.cpp: score() of UD_score.c (peak detection, DScores, UniScore)
float score_run(const Config& c, const float* mz, const float* in, int L, const int* nztab, const float* newblur,
                const std::vector<float>& massaxis, const std::vector<float>& massaxisval, const std::vector<float>& massgrid,
                int mlen, float rsquared, std::vector<float>& peakx, std::vector<float>& peaky, std::vector<float>& dscores,
                const Control& ctl);

namespace {

void setup_outputs(Run& r) {
    const Config& c = r.c;
    const int L = r.L, Z = r.Z;
    // port: the neighbour table and the previous iterate are not used again: freed before the outputs
    std::vector<int>().swap(r.closeind);
    std::vector<float>().swap(r.closearray);
    std::vector<float>().swap(r.oldblur);
    std::vector<float>().swap(r.scratch);
    if (c.peakshapeinflate != 1 && c.mzsig != 0) {
        if (c.speedyflag == 0) make_peak_shape_2d(r, 0, 0);
        else make_peak_shape_1d(r, 0, 0);
    }
    const float blurmax = vmax(r.blur.data(), r.ln);
    float cutoff = 0;
    if (blurmax != 0) cutoff = 0.000001f;
    {
        const float cut = blurmax * cutoff;
        for (long i = 0; i < r.ln; i++) if (r.blur[i] < cut) r.blur[i] = 0;
    }
    r.fitdat.assign(L, 0.0f);
    r.error = errfunspeedy(r);
    if (c.intthresh != -1)
        for (int i = 0; i < L - 1; i++)
            if (r.in[i] == 0 && r.in[i + 1] == 0) { r.fitdat[i] = 0; r.fitdat[i + 1] = 0; }
    if (c.orbimode == 1) {
        for (int i = 0; i < L; i++)
            for (int j = 0; j < Z; j++) {
                const float z = (float)r.nztab[j];
                if (z != 0) r.blur[idx2(Z, i, j)] /= z;
            }
    }
    float newblurmax = blurmax;
    if (c.rawflag == 0 || c.rawflag == 2) {
        if (c.mzsig != 0) newblurmax = reconvolve(r);
        else r.newblur = r.blur;
    }
    float massmax = c.masslb, massmin = c.massub;
    if (c.fixedmassaxis == 0) {
        const float thr = newblurmax * cutoff;
#pragma omp parallel
        {
            float tmax = c.masslb, tmin = c.massub;
#pragma omp for schedule(static)
            for (int i = 0; i < L; i++) {
                for (int j = 0; j < Z; j++) {
                    const long q = idx2(Z, i, j);
                    if (r.newblur[q] * (float)r.barr[q] > thr) {
                        float testmax = r.mtab[q] + c.psmzthresh * (float)r.nztab[j] + c.massbins;
                        float testmin = r.mtab[q] - c.psmzthresh * (float)r.nztab[j];
                        testmin = std::round(testmin / c.massbins) * c.massbins;
                        testmax = std::round(testmax / c.massbins) * c.massbins;
                        if (testmax > tmax) tmax = testmax;
                        if (testmin < tmin) tmin = testmin;
                    }
                }
            }
#pragma omp critical(udengine_massrange)
            {
                if (tmax > massmax) massmax = tmax;
                if (tmin < massmin) massmin = tmin;
            }
        }
    } else {
        massmax = c.massub;
        massmin = c.masslb;
    }
    const double mlen_d = (double)(int)((massmax - massmin) / c.massbins);
    r.mlen = (int)((massmax - massmin) / c.massbins);
    if (r.mlen < 1) {
        // "ERROR: No masses detected": a dummy axis over the mass range. port: the C code allocated the
        // mass grid for the old length (0 or less) and then wrote the new length from it
        massmax = c.massub;
        massmin = c.masslb;
        r.mlen = (int)((massmax - massmin) / c.massbins);
        if (r.mlen < 1) r.mlen = 1;
        if ((double)r.mlen * Z > (double)INT_MAX)
            throw std::invalid_argument("The UniDec engine's mass x charge grid would exceed its 32 bit index. Use a larger mass step or narrower mass and charge ranges.");
        r.massaxis.assign(r.mlen, 0.0f);
        r.massaxisval.assign(r.mlen, 0.0f);
        r.massgrid.assign((size_t)r.mlen * Z, 0.0f);
        for (int i = 0; i < r.mlen; i++) r.massaxis[i] = massmin + (float)i * c.massbins;
        r.uniscore = 0;
        return;
    }
    if (mlen_d * Z > (double)INT_MAX)
        throw std::invalid_argument("The UniDec engine's mass x charge grid would exceed its 32 bit index. Use a larger mass step or narrower mass and charge ranges.");
    r.massaxis.assign(r.mlen, 0.0f);
    r.massaxisval.assign(r.mlen, 0.0f);
    r.massgrid.assign((size_t)r.mlen * Z, 0.0f);
    for (int i = 0; i < r.mlen; i++) r.massaxis[i] = massmin + (float)i * c.massbins;
    const float* blur = (c.rawflag == 1 || c.rawflag == 3) ? r.blur.data() : r.newblur.data();
    if (c.rawflag >= 0 && c.rawflag <= 3) {
        if (c.poolflag == 0) integrate_transform(r, blur, massmax, massmin);
        else if (c.poolflag == 1) interpolate_transform(r, blur);
        else smart_transform(r, blur);   // poolflag 2 (others refused before the run)
    }
    if (c.silent == 0) {
        r.uniscore = score_run(c, r.mz, r.in, L, r.nztab.data(), r.newblur.data(), r.massaxis, r.massaxisval, r.massgrid, r.mlen,
                               r.rsquared, r.peakx_out, r.peaky_out, r.dscores_out, r.ctl);
    }
}

}  // namespace

// ---------------------------------------------------------------- PostImport, LoadConfig
void Config::post_import() {
    if (psfun == 0) mzsig = mzsig / 2.35482f;
    if (zpsfun == 0) csig = csig / 2.35482f;
    dtsig = dtsig / 2.35482f;
    numz = endz - startz + 1;
    if (linflag != -1) speedyflag = linflag != 2 ? 1 : 0;
    if (massub < 0 || masslb < 0) {
        fixedmassaxis = 1;
        massub = std::fabs(massub);
        masslb = std::fabs(masslb);
    }
    if (massbins == 0) massbins = 1;
    baselineflag = (aggressiveflag == 1 || aggressiveflag == 2) ? 1 : 0;
    if (psig < 0) mzsig /= 3;
    psmzthresh = psthresh * std::fabs(mzsig) * peakshapeinflate;
}

void Config::load_line(const std::string& key, const std::string& value) {
    const char* x = key.c_str();
    const char* y = value.c_str();
    auto has = [&](const char* k) { return std::strstr(x, k) != nullptr; };
    auto I = [&]() { return (int)std::strtol(y, nullptr, 10); };
    auto F = [&]() { return std::strtof(y, nullptr); };
    if (has("mfile")) mflag = 1;
    if (has("numit")) numit = I();
    if (has("startz")) startz = I();
    if (has("endz")) endz = I();
    if (has("zzsig")) zsig = F();
    if (has("psig")) psig = F();
    if (has("beta")) beta = F();
    if (has("mzsig")) mzsig = F();
    if (has("msig")) msig = F();
    if (has("molig")) molig = F();
    if (has("massub")) massub = F();
    if (has("masslb")) masslb = F();
    if (has("psfun")) psfun = I();
    if (has("zpsfn")) zpsfun = I();
    if (has("mtabsig")) mtabsig = F();
    if (has("massbins")) massbins = F();
    if (has("psthresh")) psthresh = F();
    if (has("speedy")) speedyflag = I();
    if (has("aggressive")) aggressiveflag = I();
    if (has("adductmass")) adductmass = F();
    if (has("rawflag")) rawflag = I();
    if (has("nativezub")) nativezub = F();
    if (has("nativezlb")) nativezlb = F();
    if (has("poolflag")) poolflag = I();
    if (has("manualfile")) manualflag = 1;
    if (has("intthresh")) intthresh = F();
    if (has("peakshapeinflate")) peakshapeinflate = F();
    if (has("isotopemode")) isotopemode = I();
    if (has("orbimode")) orbimode = I();
    if (has("imflag")) imflag = I();
    if (has("cdmsflag")) cdmsflag = I();
    if (has("linflag")) linflag = I();
    if (has("csig")) csig = F();
    if (has("dtsig")) dtsig = F();
    if (has("baselineflag")) baselineflag = I();
    if (has("noiseflag")) noiseflag = I();
    if (has("filterwidth")) filterwidth = I();
    if (has("zerolog")) zerolog = F();
    if (has("peakwindow")) peakwin = F();
    if (has("peakthresh")) peakthresh = F();
    if (has("peaknorm")) peaknorm = I();
    if (has("doubledec")) doubledec = I();
    if (has("suppression_topn")) suppression_topn = I();
    if (has("suppression_topx")) suppression_topx = F();
    if (has("suppression_percent")) suppression_percent = F();
    if (has("suppression_startit")) suppression_startit = I();
    if (has("suppression_harmonic")) suppression_harmonic = I();
    if (has("suppression_satellite")) suppression_satellite = I();
}

void Control::check() const {
    if (stop()) throw ms::Cancelled();
}

// ---------------------------------------------------------------- the run
void run(const Config& c, const float* mz, const float* intensity, int lengthmz, Output& out, const Control& ctl) {
    const auto t0 = std::chrono::steady_clock::now();
    // what this port does not have (unidec.exe would run another program part or read more files)
    if (c.imflag == 1 || c.cdmsflag == 1 || c.metamode != -2)
        throw std::invalid_argument("The UniDec engine in the library runs one dimensional mass spectra only (no ion mobility, CD-MS or MetaUniDec).");
    if (c.mflag == 1 || c.manualflag == 1)
        throw std::invalid_argument("The UniDec engine in the library has no mass list or manual assignment.");
    if (c.doubledec)
        throw std::invalid_argument("The UniDec engine in the library has no DoubleDec.");
    if (lengthmz < 2) throw std::invalid_argument("The UniDec engine needs at least two data points.");
    if (c.numz < 1) throw std::invalid_argument("The UniDec engine needs at least one charge state.");
    if (c.psfun < 0 || c.psfun > 2)   // port: exit(14) in mzpeakshape
        throw std::invalid_argument("The UniDec engine stopped because the peak shape is not valid (peak shape " + std::to_string(c.psfun) + ").");
    if (!(c.mzsig == c.mzsig) || !(c.massbins > 0) || !std::isfinite(c.massbins) || !std::isfinite(c.masslb) || !std::isfinite(c.massub))
        throw std::invalid_argument("The UniDec engine cannot use these settings (peak width, mass range or mass step not valid).");
    if (c.poolflag < 0 || c.poolflag > 2)   // port: exit(1987)
        throw std::invalid_argument("The UniDec engine cannot use this way of turning m/z into mass (poolflag " + std::to_string(c.poolflag) + ").");
    if ((double)lengthmz * (double)c.numz > (double)INT_MAX)
        throw std::invalid_argument("The UniDec engine's m/z x charge grid would exceed its 32 bit index. Use Data reduction: Automatic or a narrower m/z or charge range.");
    for (int i = 0; i < lengthmz; i++)
        if (!std::isfinite(mz[i]) || !std::isfinite(intensity[i]))
            throw std::invalid_argument("The UniDec engine cannot use data that are not finite numbers.");
#ifdef _OPENMP
    if (ctl.threads > 0) omp_set_num_threads(ctl.threads);
#endif
    Run r(c, ctl, mz, intensity, lengthmz);
    setup_ztab(r);
    setup_inputs(r);
    r.check();
    const bool bad = setup_deconvolution(r);
    if (bad)   // port: "Setup is bad. No points are allowed", then a crash while writing the outputs
        throw std::runtime_error("The UniDec engine found no point of the m/z x charge grid it can use: every data point with a mass in the range lacks the same mass at a neighbouring charge state (or has zero intensity). Check the m/z, mass and charge ranges and the minimum intensity, or use Data reduction: Automatic.");
    r.deltas.assign(lengthmz, 0.0f);
    r.denom.assign(lengthmz, 0.0f);
    for (int it = 0; it < std::abs(c.numit); it++) {
        r.check();
        run_iteration(r, it);
        if (ctl.iteration) ctl.iteration->store(it + 1, std::memory_order_relaxed);
        if (c.numit < 10 || it % 10 == 0 || it % 10 == 1 || it > 0.9 * c.numit)
            if (check_convergence(r, it) == 1) break;
    }
    r.check();
    std::vector<float>().swap(r.deltas);
    std::vector<float>().swap(r.denom);
    setup_outputs(r);
    out.fitdat.swap(r.fitdat);
    if (c.baselineflag == 1) out.baseline.swap(r.baseline);
    if (ctl.keep_grid && (c.rawflag == 0 || c.rawflag == 1)) out.grid.swap(c.rawflag == 0 ? r.newblur : r.blur);
    out.massaxis.swap(r.massaxis);
    out.massaxisval.swap(r.massaxisval);
    out.massgrid.swap(r.massgrid);
    out.peakx.swap(r.peakx_out);
    out.peaky.swap(r.peaky_out);
    out.dscores.swap(r.dscores_out);
    out.error = r.error;
    out.rsquared = r.rsquared;
    out.uniscore = r.uniscore;
    out.iterations = r.iterations;
    out.mlen = r.mlen;
    out.lengthmz = lengthmz;
    out.numz = r.Z;
    out.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

// ---------------------------------------------------------------- "%f" of a float, read back
double printed6(float v, bool ties_even) {
    if (!std::isfinite(v)) return (double)v;
    const double a = std::fabs((double)v);
    if (a >= 8388608.0) return (double)v;   // 2^23 and up: whole numbers, printed exactly
    // v has at most 24 significant bits and 10^6 = 2^6 * 15625: the product is exact in a double
    const double t = a * 1e6;
    double k = std::floor(t);
    const double frac = t - k;
    if (frac > 0.5 || (frac == 0.5 && (!ties_even || std::fmod(k, 2.0) != 0))) k += 1;
    const double res = k / 1e6;   // correctly rounded, as strtod reads the printed digits
    return std::signbit(v) ? -res : res;
}

}  // namespace udengine
}  // namespace ms

// UniDec engine for one dimensional mass spectra, built into msengine: the scores of the
// deconvolution (UD_score.c of UniDec 8.2.1: get_fwhms, uscore, mscore, csscore, fscore,
// score_from_peaks, peaks_no_score, score) and the peak detection they use (UD_analysis.c:
// peak_detect, is_peak, peak_norm). Scores: Kostelic and Marty, Methods Mol. Biol. 2022, 2500, 159-180.
//
// UniDec engine by Michael T. Marty (Marty et al., Anal. Chem. 2015, 87 (8), 4370-4376,
// DOI: 10.1021/acs.analchem.5b00140), written with Erik Marklund and Andrew Baldwin.
// Copyright (c) 2016, University of Oxford
//               2017-2025, Arizona Board of Regents on behalf of the University of Arizona
// All rights reserved. Distributed under the UniDec License (BSD with a citation requirement),
// whose full text with the conditions and the disclaimer is in udengine.h and in
// LICENSES/UniDec_LICENSE.txt of the package.
//
// port: the arithmetic is the C code's, in the same order. is_peak compared every point of the
// mass axis with every candidate (minutes on fine mass axes); here only the points within the
// window around the candidate are compared, which on the sorted axis are the same points. The
// "Peak: Mass: ... DScore: ..." lines unidec.exe printed are the arrays peakx, peaky and dscores.
#include "udengine.h"
#include "../common.h"
#include <algorithm>
#include <cmath>
#include <vector>

namespace ms {
namespace udengine {

namespace {

inline long idx2(int ncols, long r, int c) { return r * (long)ncols + c; }

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

struct In {
    const Config& c;
    const float* mz;
    const float* in;
    int L;
    const int* nztab;
    const float* newblur;
    const float* massaxis;
    const float* massaxisval;
    const float* massgrid;
    int mlen;
};

void get_fwhms(const In& d, int plen, const float* peakx, float* fwhmlow, float* fwhmhigh, int* badfwhm) {
    const float* massaxis = d.massaxis;
    const float* masssum = d.massaxisval;
    const int mlen = d.mlen;
    for (int i = 0; i < plen; i++) {
        const float peak = peakx[i];
        const int index = nearfast(massaxis, peak, mlen);
        const float max = masssum[index];
        const float halfmax = max / 2.0f;
        int lindex = index, hindex = index;
        for (;;) {
            lindex -= 1;
            if (lindex < 0) break;
            if (masssum[lindex] <= halfmax) break;
        }
        for (;;) {
            hindex += 1;
            if (hindex >= mlen) break;
            if (masssum[hindex] <= halfmax) break;
        }
        // port: the C code read massaxis[-1] or massaxis[mlen] when the peak reached the end of the axis
        float mlow = massaxis[std::max(lindex, 0)];
        float mhigh = massaxis[std::min(hindex, mlen - 1)];
        const float highdiff = mhigh - peak, lowdiff = peak - mlow, fwhm = mhigh - mlow;
        const float threshold = 0.75f;
        const float mult = threshold / (1 - threshold);
        if (fwhm != 0) {
            const float rh = highdiff / fwhm, rl = lowdiff / fwhm;
            if (rh > threshold) {
                mhigh = peak + lowdiff * mult;
                badfwhm[i] = 1;
            } else if (rl > threshold) {
                mlow = peak - highdiff * mult;
                badfwhm[i] = 1;
            }
        }
        fwhmlow[i] = mlow;
        fwhmhigh[i] = mhigh;
    }
}

float uscore(const In& d, float mlow, float mhigh) {
    const Config& c = d.c;
    float numerator = 0, denominator = 0;
    const float datamin = d.mz[0], datamax = d.mz[d.L - 1];
    for (int i = 0; i < c.numz; i++) {
        const float z = (float)d.nztab[i];
        // port: unidec.exe divides by z once (two divisions by the same value became one reciprocal)
        const float rz = 1.0f / z;
        const float lmz = (mlow + z * c.adductmass) * rz;
        const float hmz = (mhigh + z * c.adductmass) * rz;
        if (hmz > datamin && lmz < datamax) {
            int lindex = nearfast(d.mz, lmz, d.L);
            int hindex = nearfast(d.mz, hmz, d.L);
            if (d.mz[lindex] < lmz) lindex += 1;
            if (d.mz[hindex] > hmz) hindex -= 1;
            // port: the three sums in four interleaved partial sums, as unidec.exe adds them (see Sum4 in udengine.cpp)
            float se[4] = {0, 0, 0, 0}, sd[4] = {0, 0, 0, 0}, sc[4] = {0, 0, 0, 0};
            const int n = hindex - lindex + 1;
            const int n4 = n > 0 ? n - n % 4 : 0;
            auto term = [&](int j, float& e, float& dt, float& dc) {
                const float data = d.in[j];
                float decon = d.newblur[idx2(c.numz, j, i)];
                if (c.orbimode == 1) decon *= z;
                e = std::fabs(data - decon);
                dt = data;
                dc = decon;
            };
            for (int q = 0; q < n4; q++) {
                float e, dt, dc;
                term(lindex + q, e, dt, dc);
                se[q & 3] += e; sd[q & 3] += dt; sc[q & 3] += dc;
            }
            float sumerrors = (se[0] + se[2]) + (se[1] + se[3]);
            float sumdata = (sd[0] + sd[2]) + (sd[1] + sd[3]);
            float sumdecon = (sc[0] + sc[2]) + (sc[1] + sc[3]);
            for (int q = n4; q < n; q++) {
                float e, dt, dc;
                term(lindex + q, e, dt, dc);
                sumerrors += e; sumdata += dt; sumdecon += dc;
            }
            float per = 0;
            if (sumdata > 0) per = 1 - (sumerrors / sumdata);
            sumdecon = sumdecon * sumdecon;   // powf(sumdecon, 2)
            numerator += sumdecon * per;
            denominator += sumdecon;
        }
    }
    return denominator != 0 ? numerator / denominator : 0;
}

float mscore(const In& d, float mlow, float mhigh) {
    const Config& c = d.c;
    const float* massaxis = d.massaxis;
    const float* masssum = d.massaxisval;
    float numerator = 0, denominator = 0;
    const float datamin = massaxis[0], datamax = massaxis[d.mlen - 1];
    if (mhigh > datamin && mlow < datamax) {
        int lindex = nearfast(massaxis, mlow, d.mlen);
        int hindex = nearfast(massaxis, mhigh, d.mlen);
        if (massaxis[lindex] < mlow) lindex += 1;
        if (massaxis[hindex] > mhigh) hindex -= 1;
        float mmax = 0, msum = 0;
        for (int j = lindex; j <= hindex; j++) {
            const float data = masssum[j];
            if (data > mmax) mmax = data;
            msum += data;
        }
        for (int i = 0; i < c.numz; i++) {
            float sumdecon = 0;
            for (int j = lindex; j <= hindex; j++) sumdecon += d.massgrid[idx2(c.numz, j, i)];
            if (sumdecon > 0) {
                float sumerrors = 0, sumdata = 0;
                for (int j = lindex; j <= hindex; j++) {
                    const float data = masssum[j];
                    const float decon = d.massgrid[idx2(c.numz, j, i)] / sumdecon * msum;
                    sumerrors += std::fabs(data - decon);
                    sumdata += data;
                }
                float per = 0;
                if (sumdata > 0) per = 1 - (sumerrors / sumdata);
                sumdecon = sumdecon * sumdecon;
                numerator += sumdecon * per;
                denominator += sumdecon;
            }
        }
    }
    return denominator != 0 ? numerator / denominator : 0;
}

float csscore(const In& d, float mlow, float mhigh, std::vector<float>& zvals) {
    const Config& c = d.c;
    const float* massaxis = d.massaxis;
    float score = 0;
    std::fill(zvals.begin(), zvals.end(), 0.0f);
    const float datamin = massaxis[0], datamax = massaxis[d.mlen - 1];
    if (mhigh > datamin && mlow < datamax) {
        int lindex = nearfast(massaxis, mlow, d.mlen);
        int hindex = nearfast(massaxis, mhigh, d.mlen);
        if (massaxis[lindex] < mlow) lindex += 1;
        if (massaxis[hindex] > mhigh) hindex -= 1;
        float zmax = 0, zsum = 0;
        int zmaxindex = -1;
        for (int i = 0; i < c.numz; i++) {
            float zval = 0;
            for (int j = lindex; j <= hindex; j++) zval += d.massgrid[idx2(c.numz, j, i)];
            zvals[i] = zval;
            zsum += zval;
            if (zval > zmax) { zmax = zval; zmaxindex = i; }
        }
        float badarea = 0;
        int index = zmaxindex;
        float lowval = zmax;
        while (index < c.numz - 1) {
            index += 1;
            const float v = zvals[index];
            if (v < lowval) lowval = v;
            else badarea += v - lowval;
        }
        index = zmaxindex;
        lowval = zmax;
        while (index > 0) {
            index -= 1;
            const float v = zvals[index];
            if (v < lowval) lowval = v;
            else badarea += v - lowval;
        }
        if (zsum != 0) score = 1 - (badarea / zsum);
    }
    return score;
}

float find_minimum(const In& d, float lowpt, float highpt) {
    const int lindex = nearfast(d.massaxis, lowpt, d.mlen);
    const int hindex = nearfast(d.massaxis, highpt, d.mlen);
    float minval = d.massaxisval[hindex];
    for (int i = lindex; i < hindex; i++)
        if (d.massaxisval[i] < minval) minval = d.massaxisval[i];
    return minval;
}

float score_minimum(float height, float min) {
    const float hh = height / 2;
    if (min > hh) return 1 - ((min - hh) / (height - hh));
    return 1;
}

float fscore(const In& d, int plen, const float* peakx, float height, float mlow, float mhigh, float peak, int badfwhm) {
    const Config& c = d.c;
    float fs = 1;
    if (badfwhm == 1) {
        const float highdiff = mhigh - peak, lowdiff = peak - mlow;
        if (lowdiff > highdiff) fs *= score_minimum(height, find_minimum(d, mlow - c.massbins, peak));
        else fs *= score_minimum(height, find_minimum(d, peak, mhigh + c.massbins));
    }
    for (int i = 0; i < plen; i++) {
        const float peak2 = peakx[i];
        if (peak != peak2) {
            if (peak2 < peak && peak2 > mlow) fs *= score_minimum(height, find_minimum(d, peak2, peak));
            if (peak2 > peak && peak2 < mhigh) fs *= score_minimum(height, find_minimum(d, peak, peak2));
        }
    }
    return fs;
}

// peak_detect with is_peak: a point is a peak when it is not below the threshold, no point within
// +- window (in mass units) is higher, and no earlier point within it is equally high
int peak_detect(const float* x, const float* y, int n, float window, float thresh, int normthresh, float* peakx, float* peaky,
                const Control& ctl) {
    int plen = 0;
    float max = n > 0 ? y[0] : 0;
    for (int i = 0; i < n; i++) if (y[i] > max) max = y[i];
    float threshval = 0;
    if (normthresh == 1) {
        if (thresh <= 1) {
            threshval = thresh * max;
        } else {
            std::vector<float> sorted(y, y + n);
            std::sort(sorted.begin(), sorted.end());
            const float percentage = 0.9f;
            int threshold_index = 0;
            for (int i = 0; i < n; i++) {
                if (sorted[i] > 0) {
                    threshold_index = i + (int)((n - i) * percentage);
                    break;
                }
            }
            if (threshold_index < n) threshval = sorted[threshold_index] * thresh;
        }
    } else {
        threshval = thresh;
    }
    for (int index = 0; index < n; index++) {
        const float xval = x[index], yval = y[index];
        if (yval < threshval) continue;
        if ((index & 0x3ff) == 0) ctl.check();
        bool peak = true;
        // port: the axis is sorted, so |x[i] - xval| <= window holds on one run of points around index
        if (std::fabs(x[index] - xval) <= window) {
            for (int i = index; i >= 0 && std::fabs(x[i] - xval) <= window; i--) {
                if (y[i] > yval || (y[i] == yval && i < index)) { peak = false; break; }
            }
            if (peak)
                for (int i = index + 1; i < n && std::fabs(x[i] - xval) <= window; i++)
                    if (y[i] > yval) { peak = false; break; }
        }
        if (peak) {
            peakx[plen] = xval;
            peaky[plen] = yval;
            plen++;
        }
    }
    return plen;
}

void peak_norm(float* peaky, int plen, int peaknorm) {
    float norm = 0;
    if (peaknorm == 1) { norm = plen > 0 ? peaky[0] : 0; for (int i = 0; i < plen; i++) if (peaky[i] > norm) norm = peaky[i]; }
    if (peaknorm == 2) { for (int i = 0; i < plen; i++) norm += peaky[i]; }
    if (norm != 0)
        for (int i = 0; i < plen; i++) peaky[i] /= norm;
}

}  // namespace

float score_run(const Config& c, const float* mz, const float* in, int L, const int* nztab, const float* newblur,
                const std::vector<float>& massaxis, const std::vector<float>& massaxisval, const std::vector<float>& massgrid,
                int mlen, float rsquared, std::vector<float>& peakx, std::vector<float>& peaky, std::vector<float>& dscores,
                const Control& ctl) {
    peakx.assign(mlen, 0.0f);
    peaky.assign(mlen, 0.0f);
    const int plen = peak_detect(massaxis.data(), massaxisval.data(), mlen, c.peakwin, c.peakthresh, c.normthresh, peakx.data(),
                                 peaky.data(), ctl);
    peakx.resize(plen);
    peaky.resize(plen);
    peakx.shrink_to_fit();
    peaky.shrink_to_fit();
    dscores.assign(plen, 0.0f);
    if (plen == 0) return 0;
    peak_norm(peaky.data(), plen, c.peaknorm);
    In d{c, mz, in, L, nztab, newblur, massaxis.data(), massaxisval.data(), massgrid.data(), mlen};
    std::vector<float> fwhmlow(plen), fwhmhigh(plen);
    std::vector<int> badfwhm(plen, 0);
    get_fwhms(d, plen, peakx.data(), fwhmlow.data(), fwhmhigh.data(), badfwhm.data());
    float numerator = 0, denominator = 0, uniscore = 0;
    const float threshold = 0;
    std::vector<float> zvals(c.numz);
    for (int i = 0; i < plen; i++) {
        if ((i & 0x3f) == 0) ctl.check();
        const float xfwhm = 2;
        const float m = peakx[i];
        const float ival = peaky[i];
        const float l = m - (m - fwhmlow[i]) * xfwhm;
        const float h = m + (fwhmhigh[i] - m) * xfwhm;
        const int index = nearfast(massaxis.data(), m, mlen);
        const float height = massaxisval[index];
        const float usc = uscore(d, l, h);
        const float msc = mscore(d, l, h);
        const float cssc = csscore(d, l, h, zvals);
        const float fsc = fscore(d, plen, peakx.data(), height, fwhmlow[i], fwhmhigh[i], m, badfwhm[i]);
        const float dsc = usc * msc * cssc * fsc;
        dscores[i] = dsc;
        if (dsc > threshold) {
            numerator += ival * ival * dsc;
            denominator += ival * ival;
        }
    }
    if (denominator != 0) uniscore = rsquared * numerator / denominator;
    return uniscore;
}

}  // namespace udengine
}  // namespace ms

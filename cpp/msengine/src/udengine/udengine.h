// UniDec engine for one dimensional mass spectra, built into msengine.
//
// A C++ port of the engine of UniDec 8.2.1 (unidec/src: udmain.c, udcore.c, udtools.c,
// udstruct.c, udio.c, UD_score.c, UD_analysis.c): what unidec.exe does for a conf.dat
// without ion mobility, CD-MS, mass list, manual assignment or DoubleDec, run on arrays
// in memory instead of files. See udengine.cpp for what was changed and why.
//
// UniDec engine by Michael T. Marty (Marty et al., Anal. Chem. 2015, 87 (8), 4370-4376,
// DOI: 10.1021/acs.analchem.5b00140), written with Erik Marklund and Andrew Baldwin.
// If you publish results obtained with it, please cite that paper.
//
// UniDec License:
//
// Copyright (c) 2016, University of Oxford
//               2017-2025, Arizona Board of Regents on behalf of the University of Arizona
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 1. Redistributions of source code must retain the above copyright
//    notice, this list of conditions, and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions, and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
// 3. Neither the name of the copyright holders nor the
//    names of its contributors may be used to endorse or promote products
//    derived from this software without specific prior written permission.
// 4. Any publications that result from use of the software should cite Marty et al. Anal. Chem. 2015.
//    DOI: 10.1021/acs.analchem.5b00140. If UniDec is redistributed or incorporated into other software,
//    it must be clearly indicated to the end user that UniDec is being used, and the request to cite
//    Marty et al. Anal. Chem. 2015. DOI: 10.1021/acs.analchem.5b00140 must be passed on to the end user.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ''AS IS'' AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#pragma once
#include <atomic>
#include <string>
#include <vector>

namespace ms {
namespace udengine {

// The settings of the engine (struct Config of udstruct.h, the fields a one dimensional run reads),
// with the defaults of SetDefaultConfig. Fill it as LoadConfig would from conf.dat (floats as
// strtof reads the text), then call post_import() (PostImport).
struct Config {
    int numit = 50, numz = 0, endz = 100, startz = 1;
    float zsig = 1, psig = 1, beta = 0, mzsig = 1, msig = 0, molig = 0, massub = 5000000, masslb = 100;
    int psfun = 0, zpsfun = 0;
    float psmzthresh = 0, mtabsig = 0, massbins = 100, psthresh = 6;
    int speedyflag = 0, linflag = -1, aggressiveflag = 0;
    float adductmass = 1.007276467f;
    int rawflag = 1;
    float nativezub = 100, nativezlb = -200;
    int poolflag = 1;
    float intthresh = 0, peakshapeinflate = 1;
    int fixedmassaxis = 0, isotopemode = 0, baselineflag = 1, noiseflag = 0, orbimode = 0, datanorm = 1;
    int filterwidth = 20;
    float zerolog = -12;
    float peakwin = 500, peakthresh = 0.1f;
    int peaknorm = 1, normthresh = 1, doubledec = 0, silent = 0;
    int suppression_topn = 0;
    float suppression_topx = 0, suppression_percent = 0;
    int suppression_startit = 3, suppression_harmonic = 0, suppression_satellite = 0;
    // features of the engine this port does not have (refused with a plain message when set)
    int mflag = 0, manualflag = 0, imflag = 0, cdmsflag = 0, metamode = -2;
    float csig = 1, dtsig = 0.2f;

    void post_import();                       // PostImport of udstruct.c
    // one "key value" line of conf.dat as LoadConfig reads it (substring match of the key, as the engine does)
    void load_line(const std::string& key, const std::string& value);
};

// What the engine computes (the arrays unidec.exe writes as _fitdat.bin, _baseline.bin, _grid.bin,
// _massgrid.bin and _mass.txt, and the values of _error.txt), as floats.
struct Output {
    std::vector<float> fitdat;               // lengthmz
    std::vector<float> baseline;             // lengthmz (baselineflag 1 only)
    std::vector<float> grid;                 // lengthmz x numz (only when Control::keep_grid; rawflag 0: newblur, 1: blur)
    std::vector<float> massaxis, massaxisval;// mlen
    std::vector<float> massgrid;             // mlen x numz
    std::vector<float> peakx, peaky, dscores;// the engine's peaks with their DScores (UD_score.c)
    float error = 0, rsquared = 0, uniscore = 0;
    int iterations = 0, mlen = 0, lengthmz = 0, numz = 0;
    double seconds = 0;                      // wall time of the deconvolution
};

// The caller's handles: a cancel flag the engine checks between iterations and inside its long
// loops, the iterations done (for the progress display), and the thread count.
struct Control {
    const std::atomic<int>* cancel = nullptr;
    std::atomic<long>* iteration = nullptr;
    int threads = 0;                          // 0: OpenMP default
    bool keep_grid = false;                   // fill Output::grid
    bool stop() const { return cancel && cancel->load(std::memory_order_relaxed) != 0; }
    void check() const;                       // throws ms::Cancelled when cancelled
};

// Runs the engine (run_unidec of udmain.c without the files) on lengthmz points (m/z ascending,
// as floats the way the engine reads its input). Throws std::invalid_argument for settings the
// engine cannot use, std::runtime_error with a plain message where the engine would stop or crash,
// ms::Cancelled when cancelled, std::bad_alloc when memory runs out. No exit(), no global state.
void run(const Config& config, const float* mz, const float* intensity, int lengthmz, Output& out, const Control& ctl);

// The text "%f" makes of a float, read back as a double (what the readers of _mass.txt, _error.txt
// and the engine's peak lines get). ties_even: how an exact tie at the 7th decimal is rounded.
double printed6(float v, bool ties_even = true);

}  // namespace udengine
}  // namespace ms

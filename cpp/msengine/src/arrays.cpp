// Scans handed over as arrays (ms_open_arrays): LC-MS files whose format is read
// in Python (Agilent ChemStation .D and Waters .raw through rainbow, Thermo .raw
// through RawFileReader, mzXML, ANDI netCDF). The Python reference is
// lcms_sources.ArrayLCDFile (lcms_data.LCDFile on the same arrays).
//
// The spectra are kept as float32 and summed on m/z bins, exactly as the
// ShimadzuReader does it (LCDFile._binned + _fill, LCDFile.chromatogram "xic"),
// so that both give the same values as the Python reference bit for bit.
#include "file.h"
#include <algorithm>
#include <cmath>
#include <cstdio>

namespace ms {
namespace {

const long long kMaxBins = 20000000;   // lcms_data.MAX_BINS

class ArrayReader : public Reader {
public:
    ArrayReader(const std::string& kind, const std::string& instrument, long n, const double* rt, const int* event,
                const long long* offsets, const double* mz, const double* it, const double* tic, const double* bpc,
                int n_events, const int* ev_pol, const double* ev_lo, const double* ev_hi, long n_msms)
        : kind_(kind), instrument_(instrument), n_msms_(n_msms) {
        if (n <= 0 || !rt || !event || !offsets || !tic || !bpc) throw std::invalid_argument("ms_open_arrays: no scans");
        if (n_events <= 0 || !ev_pol) throw std::invalid_argument("ms_open_arrays: no scan events");
        if (offsets[0] != 0) throw std::invalid_argument("ms_open_arrays: offsets must start at 0");
        for (long k = 0; k < n; k++)
            if (offsets[k + 1] < offsets[k]) throw std::invalid_argument("ms_open_arrays: offsets must not decrease");
        const long long total = offsets[n];
        if (total > 0 && (!mz || !it)) throw std::invalid_argument("ms_open_arrays: null spectra");
        for (int e = 0; e < n_events; e++) {
            ScanInfo ev;
            ev.polarity = ev_pol[e] < 0 ? -1 : 1;
            ev.mz_lo = ev_lo ? ev_lo[e] : 0;
            ev.mz_hi = ev_hi ? ev_hi[e] : 0;
            ev.profile = false;
            events_.push_back(ev);
        }
        mz_.resize((size_t)total);
        it_.resize((size_t)total);
        for (long long j = 0; j < total; j++) {
            mz_[j] = (float)mz[j];
            it_[j] = (float)it[j];
        }
        offsets_.assign(offsets, offsets + n + 1);
        scans_.resize((size_t)n);
        for (long k = 0; k < n; k++) {
            if (event[k] < 0 || event[k] >= n_events) throw std::invalid_argument("ms_open_arrays: scan event out of range");
            if (k && rt[k] < rt[k - 1]) throw std::invalid_argument("ms_open_arrays: the scans must be in time order");
            ScanInfo& s = scans_[(size_t)k];
            s.rt = rt[k];
            s.event = event[k];
            s.polarity = events_[event[k]].polarity;
            s.profile = false;
            s.mz_lo = events_[event[k]].mz_lo;
            s.mz_hi = events_[event[k]].mz_hi;
            s.tic = tic[k];
            s.bpc = bpc[k];
        }
    }

    const char* kind() const override { return kind_.c_str(); }
    const std::string& instrument() const override { return instrument_; }
    const std::string& recalibration() const override { return empty_; }
    long n_msms() const override { return n_msms_; }
    const std::vector<ScanInfo>& scans() const override { return scans_; }
    const std::vector<ScanInfo>* events() const override { return &events_; }

    void spectrum(long i, Spectrum& out) override {
        out.profile = false;
        centroids(i, out.mz, out.it);
    }
    void centroids(long i, std::vector<double>& mz, std::vector<double>& it) override {
        if (i < 0 || (size_t)i + 1 >= offsets_.size()) throw std::invalid_argument("scan index out of range");
        const size_t a = (size_t)offsets_[i], b = (size_t)offsets_[i + 1];
        mz.resize(b - a);
        it.resize(b - a);
        for (size_t j = a; j < b; j++) { mz[j - a] = mz_[j]; it[j - a] = it_[j]; }
    }

    // LCDFile.chromatogram(kind="xic"): |mz32 - mz| <= tol in float32, intensities summed
    bool xic(int, double mz, double tol, std::vector<double>& y) override {
        const size_t n = scans_.size();
        y.assign(n, 0.0);
        const float m = (float)mz, t = (float)tol;
        for (size_t k = 0; k < n; k++) {
            double s = 0;
            for (size_t j = (size_t)offsets_[k]; j < (size_t)offsets_[k + 1]; j++)
                if (std::fabs(mz_[j] - m) <= t) s += (double)it_[j];
            y[k] = s;
        }
        return true;
    }

    // LCDFile._binned + _fill: the plain sum per bin (file.cpp divides by the scan
    // count), bin m/z = intensity weighted mean, empty bins of the grid at zero
    bool sum_profiles(const std::vector<long>& idx, Spectrum& out, ms_progress_fn cb, void* user) override {
        out.profile = false;
        out.mz.clear();
        out.it.clear();
        if (idx.empty()) return true;
        const double binw = binw_ > 0 ? binw_ : 0.05;
        long long kmin = 0, kmax = 0;
        bool any = false;
        for (long i : idx) {
            if (i < 0 || (size_t)i + 1 >= offsets_.size()) throw std::invalid_argument("scan index out of range");
            for (size_t j = (size_t)offsets_[i]; j < (size_t)offsets_[i + 1]; j++) {
                const long long key = (long long)std::nearbyint((double)mz_[j] / binw);
                if (!any) { kmin = kmax = key; any = true; }
                else { kmin = std::min(kmin, key); kmax = std::max(kmax, key); }
            }
        }
        if (!any) return true;
        if (kmax - kmin + 1 > kMaxBins) {
            char buf[200];
            std::snprintf(buf, sizeof buf, "the bin width %g m/z is too small for m/z %g to %g (%lld bins, at most %lld)",
                          binw, kmin * binw, kmax * binw, kmax - kmin + 1, kMaxBins);
            throw Error(buf);
        }
        const size_t nb = (size_t)(kmax - kmin + 1);
        std::vector<double> s_i(nb, 0.0), s_m(nb, 0.0);
        std::vector<char> used(nb, 0);
        const long n = (long)idx.size();
        for (long k = 0; k < n; k++) {
            if (k % 200 == 0) tick(cb, user, k, n);
            const long i = idx[k];
            for (size_t j = (size_t)offsets_[i]; j < (size_t)offsets_[i + 1]; j++) {
                const double m = mz_[j], w = it_[j];
                const size_t b = (size_t)((long long)std::nearbyint(m / binw) - kmin);
                s_i[b] += w;
                s_m[b] += m * w;
                used[b] = 1;
            }
        }
        out.mz.resize(nb + 2);
        out.it.assign(nb + 2, 0.0);
        for (size_t g = 0; g < nb + 2; g++) out.mz[g] = (double)(kmin - 1 + (long long)g) * binw;
        for (size_t b = 0; b < nb; b++) {
            if (!used[b]) continue;
            out.mz[b + 1] = s_i[b] > 0 ? s_m[b] / s_i[b] : (double)(kmin + (long long)b) * binw;
            out.it[b + 1] = s_i[b];
        }
        return true;
    }

    double bin_width() const override { return binw_ > 0 ? binw_ : 0.05; }

    bool set_bin_width(double binw) override {
        if (!(binw >= 0) || binw > 10 || (binw > 0 && binw < 1e-6))
            throw std::invalid_argument("bin width must be 0 (the default 0.05 Da) or between 0.000001 and 10 Da");
        binw_ = binw;
        return true;
    }

private:
    std::string kind_, instrument_, empty_;
    long n_msms_ = 0;
    double binw_ = 0;
    std::vector<float> mz_, it_;
    std::vector<long long> offsets_;
    std::vector<ScanInfo> scans_, events_;
};

}  // namespace

std::unique_ptr<Reader> open_arrays(const std::string& kind, const std::string& instrument, long n, const double* rt,
                                    const int* event, const long long* offsets, const double* mz, const double* it,
                                    const double* tic, const double* bpc, int n_events, const int* ev_pol,
                                    const double* ev_lo, const double* ev_hi, long n_msms) {
    return std::unique_ptr<Reader>(new ArrayReader(kind, instrument, n, rt, event, offsets, mz, it, tic, bpc, n_events,
                                                   ev_pol, ev_lo, ev_hi, n_msms));
}

}  // namespace ms

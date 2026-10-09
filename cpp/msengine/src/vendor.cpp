// Data files read by the caller (ms_open_vendor): Agilent MassHunter .D, Waters
// MassLynx .raw and Thermo .raw are read in Python with the vendors' libraries or
// rainbow (hrms_vendor.py); this reader holds the scan list the caller gives and
// asks the caller for each spectrum, so that the averages, scans and mass
// chromatograms of these files run in file.cpp as for every other format.
// Mirrors hrms_vendor.VendorFile of the Python reference: the centroids of every
// scan are read when the file opens (the caller's centroids, or the profile
// centroided as hrms_data.centroid does), TIC and BPC from them where the caller
// gives none.
//
// Also the LZF decompression of Agilent's MSProfile.bin (ms_lzf_decompress;
// Python reference: _portable/lzf.py, decompress_py).
#include "file.h"
#include <algorithm>
#include <cmath>
#include <cstring>

namespace ms {
namespace {

class VendorReader : public Reader {
public:
    VendorReader(const char* kind, const char* instrument, long n_msms, const ms_vendor_scan* scans, long n,
                 const ms_vendor_scan* events, int n_events, bool common, ms_vendor_fn fn, void* user,
                 ms_progress_fn cb, void* cbuser)
        : kind_(kind && *kind ? kind : "vendor file"), instrument_(instrument ? instrument : ""),
          n_msms_(n_msms > 0 ? n_msms : 0), common_(common), fn_(fn), user_(user) {
        if (!fn) throw std::invalid_argument("ms_open_vendor: null spectrum callback");
        if (n <= 0 || !scans) throw std::invalid_argument("ms_open_vendor: no scans");
        if (n_events < 0 || (n_events > 0 && !events)) throw std::invalid_argument("ms_open_vendor: bad events");
        for (int e = 0; e < n_events; e++) {
            ScanInfo ev;
            ev.polarity = events[e].polarity < 0 ? -1 : 1;
            ev.mz_lo = std::isfinite(events[e].mz_lo) && events[e].mz_lo > 0 ? events[e].mz_lo : 0;
            ev.mz_hi = std::isfinite(events[e].mz_hi) && events[e].mz_hi > 0 ? events[e].mz_hi : 0;
            events_.push_back(ev);
        }
        scans_.resize((size_t)n);
        for (long i = 0; i < n; i++) {
            const ms_vendor_scan& s = scans[i];
            if (!std::isfinite(s.rt)) throw std::invalid_argument("ms_open_vendor: a scan time is not a finite number");
            if (i && s.rt < scans[i - 1].rt) throw std::invalid_argument("ms_open_vendor: the scans are not in time order");
            if (n_events > 0 && (s.event < 0 || s.event >= n_events))
                throw std::invalid_argument("ms_open_vendor: a scan names no event given");
            ScanInfo& o = scans_[i];
            o.rt = s.rt;
            o.polarity = n_events > 0 ? events_[s.event].polarity : (s.polarity < 0 ? -1 : 1);
            o.profile = s.profile != 0;
            o.event = n_events > 0 ? s.event : -1;
            o.mz_lo = std::isfinite(s.mz_lo) && s.mz_lo > 0 ? s.mz_lo : 0;
            o.mz_hi = std::isfinite(s.mz_hi) && s.mz_hi > 0 ? s.mz_hi : 0;
            o.tic = s.tic;
            o.bpc = s.bpc;
        }
        // centroids of every scan (mass chromatograms, centroid sums)
        coff_.assign(1, 0);
        std::vector<double> cm, ci;
        Spectrum sp;
        for (long i = 0; i < n; i++) {
            if (i % 200 == 0) tick(cb, cbuser, i, n);
            ScanInfo& o = scans_[i];
            fetch(i, 1, cm, ci);
            if (cm.empty()) {
                fetch(i, 0, sp.mz, sp.it);
                if (o.profile) centroid(sp.mz.data(), sp.it.data(), (long)sp.mz.size(), 0.002, cm, ci);
                else { cm.swap(sp.mz); ci.swap(sp.it); }
            }
            cmz_.insert(cmz_.end(), cm.begin(), cm.end());
            cit_.insert(cit_.end(), ci.begin(), ci.end());
            coff_.push_back(cmz_.size());
            if (!(o.tic >= 0) || !(o.bpc >= 0)) {   // negative or NaN: from the centroids
                double s = 0, b = 0;
                for (double v : ci) { s += v; b = std::max(b, v); }
                if (!(o.tic >= 0)) o.tic = s;
                if (!(o.bpc >= 0)) o.bpc = b;
            }
        }
        tick(cb, cbuser, n, n);
    }

    const char* kind() const override { return kind_.c_str(); }
    const std::string& instrument() const override { return instrument_; }
    const std::string& recalibration() const override { return empty_; }
    long n_msms() const override { return n_msms_; }
    const std::vector<ScanInfo>& scans() const override { return scans_; }
    const std::vector<ScanInfo>* events() const override { return events_.empty() ? nullptr : &events_; }
    bool common_axis() const override { return common_; }

    void spectrum(long i, Spectrum& out) override {
        check(i);
        fetch(i, 0, out.mz, out.it);
        out.profile = scans_[i].profile;
    }

    void centroids(long i, std::vector<double>& mz, std::vector<double>& it) override {
        check(i);
        mz.assign(cmz_.begin() + (long)coff_[i], cmz_.begin() + (long)coff_[i + 1]);
        it.assign(cit_.begin() + (long)coff_[i], cit_.begin() + (long)coff_[i + 1]);
    }

private:
    std::string kind_, instrument_, empty_;
    long n_msms_;
    bool common_;
    ms_vendor_fn fn_;
    void* user_;
    std::vector<ScanInfo> scans_, events_;
    std::vector<double> cmz_, cit_;
    std::vector<size_t> coff_;

    void check(long i) const {
        if (i < 0 || i >= (long)scans_.size()) throw std::invalid_argument("vendor file: scan index out of range");
    }

    // the caller's arrays of scan i, copied at once (they stay the caller's)
    void fetch(long i, int what, std::vector<double>& mz, std::vector<double>& it) {
        const double* pm = nullptr;
        const double* pi = nullptr;
        long n = 0;
        if (fn_(user_, i, what, &pm, &pi, &n) != 0)
            throw Error("scan " + std::to_string(i + 1) + " of the " + kind_ + " file could not be read");
        if (n < 0 || (n > 0 && (!pm || !pi))) throw Error("the reader of the " + kind_ + " file gave no spectrum");
        mz.assign(pm, pm + n);
        it.assign(pi, pi + n);
    }
};

}  // namespace

std::unique_ptr<Reader> open_vendor(const char* kind, const char* instrument, long n_msms, const ms_vendor_scan* scans,
                                    long n, const ms_vendor_scan* events, int n_events, bool common_axis,
                                    ms_vendor_fn fn, void* user, ms_progress_fn cb, void* cbuser) {
    return std::unique_ptr<Reader>(new VendorReader(kind, instrument, n_msms, scans, n, events, n_events, common_axis,
                                                    fn, user, cb, cbuser));
}

}  // namespace ms

extern "C" {

// liblzf's format: a control byte below 32 starts a run of ctrl + 1 literal bytes; otherwise the top
// three bits give the length (7: one more byte follows to add) and the rest with the next byte the
// distance back (+ 1) of a copy of length + 2 bytes from the output written so far.
MS_API int ms_lzf_decompress(const unsigned char* in, long n_in, unsigned char* out, long cap, long* n_out) {
    if (n_out) *n_out = 0;
    bool too_small = false;
    const int rc = ms::guarded([&] {
        if (!n_out || n_in < 0 || cap < 0 || (n_in > 0 && !in) || (cap > 0 && !out))
            throw std::invalid_argument("ms_lzf_decompress: bad arguments");
        long ip = 0, op = 0;
        while (ip < n_in) {
            const unsigned ctrl = in[ip++];
            if (ctrl < 32) {
                const long len = (long)ctrl + 1;
                if (ip + len > n_in) throw ms::Error("damaged LZF data (a literal run goes past the end)");
                if (op + len > cap) { too_small = true; return; }
                std::memcpy(out + op, in + ip, (size_t)len);
                op += len;
                ip += len;
            } else {
                long len = (long)(ctrl >> 5);
                long ref = op - (long)((ctrl & 0x1f) << 8) - 1;
                if (len == 7) {
                    if (ip >= n_in) throw ms::Error("damaged LZF data (truncated)");
                    len += in[ip++];
                }
                if (ip >= n_in) throw ms::Error("damaged LZF data (truncated)");
                ref -= in[ip++];
                len += 2;
                if (ref < 0) throw ms::Error("damaged LZF data (a copy starts before the output)");
                if (op + len > cap) { too_small = true; return; }
                for (long k = 0; k < len; k++) out[op + k] = out[ref + k];   // the copy may overlap itself
                op += len;
            }
        }
        *n_out = op;
    });
    if (rc == MS_OK && too_small) {
        ms::set_error("the LZF data need a larger output");
        return MS_NOT_FOUND;
    }
    return rc;
}

}  // extern "C"

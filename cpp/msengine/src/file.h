// Internal reader interface: one implementation per file format. file.cpp
// implements the ms_file API (events, chromatograms, averaging) on top of it.
//
// Threads: file.cpp calls the methods of one Reader one at a time (it holds the
// file's reader mutex around every call), so an implementation needs no locking
// of its own for the GUI's pattern (an average in a worker thread while the main
// thread asks for scans and chromatograms). The progress callback a method gets
// may be called with that mutex held; it is a wrapper of file.cpp that returns 0
// (cancel) once the file is being closed. Methods may throw (ms::Error,
// std::invalid_argument, ms::Cancelled through tick(), std::bad_alloc); file.cpp
// turns that into return codes. Readers must not throw from inside a callback that
// a vendor library calls (the exception cannot pass through the library's frames).
#pragma once
#include "common.h"
#include <memory>

namespace ms {

struct ScanInfo {
    double rt = 0;        // minutes
    int polarity = 1;     // +1 / -1
    bool profile = true;  // this scan has a profile spectrum
    double mz_lo = 0, mz_hi = 0;   // acquisition range, 0 = unknown
    double tic = 0, bpc = 0;
    int event = -1;       // index into Reader::events() when the reader defines events
    bool saturated = false;
};

class Reader {
public:
    virtual ~Reader() {}
    virtual const char* kind() const = 0;
    virtual const std::string& instrument() const = 0;          // "" if unknown
    virtual const std::string& recalibration() const = 0;       // "" if none
    virtual long n_msms() const = 0;
    virtual const std::vector<ScanInfo>& scans() const = 0;     // MS1 scans in time order (fixed after opening)
    // the profile (or, when the scan has none, the centroid) spectrum of scan i
    virtual void spectrum(long i, Spectrum& out) = 0;
    // the centroid (line) spectrum of scan i, for mass chromatograms
    virtual void centroids(long i, std::vector<double>& mz, std::vector<double>& it) = 0;
    // true when every profile scan shares the same m/z axis (sum point by point)
    virtual bool common_axis() const { return false; }
    // a reader may provide its own fast mass chromatogram (timsTOF SDK); default: from the centroids.
    // y has one value per entry of scans() (every MS1 scan, scans of the other polarity 0).
    // Shimadzu: PDA matrix (times x wavelengths, mAU) and sample information ("key\tvalue\n" lines)
    virtual bool pda(std::vector<double>&, std::vector<double>&, std::vector<double>&) { return false; }
    virtual std::string sample_info() { return ""; }
    // scan events as the file defines them (Shimadzu: the method's events, each with its polarity and
    // m/z range; the ScanInfo.event index refers to this list). Empty: events are the polarities.
    virtual const std::vector<ScanInfo>* events() const { return nullptr; }
    virtual bool xic(int /*polarity*/, double /*mz*/, double /*tol*/, std::vector<double>& /*y_per_scan*/) { return false; }
    // a reader may sum profile scans on its native axis (TSF/TDF: the digitizer sample
    // axis, converted to m/z with the calibration of the middle scan, as hrms_data does).
    // out = the plain sum (not divided by the count). false: sum spectrum(i) point by point.
    virtual bool sum_profiles(const std::vector<long>& /*idx*/, Spectrum& /*out*/, ms_progress_fn, void*) { return false; }
    virtual bool set_bin_width(double /*binw*/) { return false; }   // Shimadzu: bin width of the sums
    virtual double bin_width() const { return 0; }                // >0: sums are on bins of this width
};

// factories (mzml.cpp, bruker.cpp); throw ms::Error on failure
std::unique_ptr<Reader> open_mzml(const std::string& path, ms_progress_fn cb, void* user);
// Shimadzu .lcd (shimadzu.cpp); the reader also implements the optional
// pda() and sample_info() below
std::unique_ptr<Reader> open_shimadzu(const std::string& path, ms_progress_fn cb, void* user);
std::unique_ptr<Reader> open_bruker(const std::string& d_folder, const std::string& sdk_dir, ms_progress_fn cb, void* user);
bool is_bruker_d(const std::string& path, std::string* d_folder);   // accepts the folder, a file inside it, or a folder with one .d

}  // namespace ms

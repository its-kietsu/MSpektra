// Internal helpers of msengine.
#pragma once
#include "msengine.h"
#include <string>
#include <vector>
#include <new>
#include <stdexcept>
#include <cmath>

namespace ms {

// error of the last failed call in this thread (a fixed buffer: never allocates, never throws)
void set_error(const char* msg) noexcept;
void set_error(const std::string& msg) noexcept;
void clear_error() noexcept;
struct Error : std::runtime_error { using std::runtime_error::runtime_error; };
struct Cancelled : std::exception {   // the progress callback returned 0 (or the file is being closed)
    const char* msg;
    explicit Cancelled(const char* m = "cancelled") : msg(m) {}
    const char* what() const noexcept override { return msg; }
};

// progress helper: throws Cancelled when the callback returns 0
inline void tick(ms_progress_fn cb, void* user, long i, long n) {
    if (cb && !cb(user, i, n)) throw Cancelled();
}

int threads();  // the thread count in use (1 .. max_threads())
int max_threads();       // the largest count ms_set_threads accepts (4 x the processors, at least 8, at most 256)
double memory_total();   // physical memory in bytes (0 if unknown)
double memory_budget();  // what one deconvolution may use: half of it, at least 1 GB
std::string fmt_bytes(double b);  // "1.5 GB"

// the exception in flight -> return code + ms_last_error (call inside a catch block only)
int error_code_of_current_exception() noexcept;

// runs f() and turns exceptions into return codes + ms_last_error. The message is
// cleared first, so after the call ms_last_error() describes this call ("" on success).
template <class F> int guarded(F&& f) noexcept {
    clear_error();
    try { f(); return MS_OK; }
    catch (...) { return error_code_of_current_exception(); }
}

// for exported functions that return a value instead of a code: fallback when f throws
// (the message is set; ms_last_error is left as it is on success)
template <class T, class F> T guarded_value(T fallback, F&& f) noexcept {
    try { return f(); }
    catch (...) { error_code_of_current_exception(); return fallback; }
}

struct Spectrum {              // profile or centroid spectrum
    std::vector<double> mz, it;
    bool profile = true;
};

// ---- spectrum.cpp
void centroid(const double* mz, const double* it, long n, double rel, std::vector<double>& cmz, std::vector<double>& cit);
double apex(const double* mz, const double* it, long n, double x);
double resolution(const double* mz, const double* it, long n);

// Streaming sum of profile spectra whose m/z axes differ: bins as wide as the
// point spacing (relative), memory of one spectrum.
class LogBins {
public:
    // (a spacing below 1e-8 exists in no spectrum; it would make the bin keys overflow a 32 bit long)
    explicit LogBins(double r) : r_(r > 1e-8 ? r : 1e-8) {}
    void add(const double* mz, const double* it, long n);
    void result(std::vector<double>& mz, std::vector<double>& it) const;
private:
    double r_;
    long k0_ = 0; bool empty_ = true;
    std::vector<double> si_, sw_, sm_;
};

}  // namespace ms

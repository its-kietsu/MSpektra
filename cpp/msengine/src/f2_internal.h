// Internal helpers shared by unidec.cpp and isodec.cpp (agent F2): parameter
// strings, Python style number formatting, paths and process running.
#pragma once
#include "common.h"
#include "e2_internal.h"
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace ms {
namespace e2 {
// group_isotopes (peaks.cpp), and the mass of the first and last isotope peak of each species returned
std::vector<Peak> group_isotopes_ranges(const double* x, const double* y, size_t n, double threshold,
                                        std::vector<std::pair<double, double>>* ranges, double spacing = 1.00235);
}  // namespace e2
}  // namespace ms

namespace ms {
namespace unidec {

// "key=value" (or "key value") lines; has(): present and not empty
struct Params {
    std::map<std::string, std::string> kv;
    bool has(const std::string& k) const;
    std::string str(const std::string& k, const std::string& d = "") const;
    double num(const std::string& k, double d) const;
    long inum(const std::string& k, long d) const;
};
Params parse_params(const char* text);
std::string trim(const std::string& s);

std::string py_float(double v);               // repr() of a Python float
std::string md5_hex(const std::string& msg);
std::string abs_path(const std::string& p);
std::string join_path(const std::string& a, const std::string& b);
bool is_dir(const std::string& p);
void make_dirs(const std::string& path);
std::string temp_dir();
void unidec_place(std::string& folder, std::string& name);
// runs exe with args, stdout + stderr to log_path; polls cb(user, cb_i, cb_n) while it runs
// (0 kills the process and throws Cancelled); returns the exit code
long run_process(const std::string& exe, const std::vector<std::string>& args, const std::string& log_path,
                 ms_progress_fn cb, void* user, long cb_i, long cb_n);
std::vector<double> np_arange(double start, double stop, double step);
void gaussian_filter_reflect(std::vector<double>& v, double sigma);

}  // namespace unidec
}  // namespace ms

// UniDec: port of ms_deconv.unidec_deconvolute and of what UniDec's Python
// package does around the engine (engine.py process_data / run_unidec /
// pick_peaks, tools.dataprep, auto_peak_width, unidecstructure.config_export),
// with UniDec's engine itself built into the library (udengine/, a port of
// unidec.exe's C source; UniDec engine by Michael T. Marty, Marty et al.,
// Anal. Chem. 2015).
//
// Flow: the spectrum is restricted, the minimum intensity and baseline of
// run_method are applied, the data are binned as dataprep does (linear or
// constant resolution bins, normalised to 1), the configuration UniDec would
// write as <name>_conf.dat is read by the engine as unidec.exe reads it, and
// the engine runs in this process on the data as unidec.exe would have read
// them from <name>_input.dat. Its outputs are taken as the files _mass.txt,
// _fitdat.bin, _massgrid.bin and _error.txt and its peak lines would have
// given them (the text values rounded as printed). No file is written unless
// the parameter write_files is 1 (then the files UniDec writes are written to
// <work_dir>/<name>_unidecfiles). Peaks are picked as UniDec's peakdetect
// does, refined to the apex, scaled to the heights of the charge states in the
// spectrum, and grouped into species when the isotopes are resolved.
#include "e2_internal.h"
#include "f2_internal.h"
#include "udengine/udengine.h"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <exception>
#include <mutex>
#include <thread>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cstdint>
#include <cstdarg>
#include <cfloat>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <sys/stat.h>
#ifdef _WIN32
#include <windows.h>
#include <direct.h>
#else
#include <unistd.h>
#include <climits>
#endif

using ms::e2::Peak;
using ms::e2::Spec;

namespace ms {
namespace unidec {

// ============================================================ parameters ("key=value" lines)
bool Params::has(const std::string& k) const {
    auto it = kv.find(k);
    return it != kv.end() && !it->second.empty();
}
std::string Params::str(const std::string& k, const std::string& d) const {
    auto it = kv.find(k);
    return it == kv.end() || it->second.empty() ? d : it->second;
}
double Params::num(const std::string& k, double d) const {
    if (!has(k)) return d;
    const std::string& s = kv.at(k);
    char* e = nullptr;
    double v = strtod(s.c_str(), &e);
    if (e == s.c_str()) throw std::invalid_argument("parameter " + k + ": not a number: " + s);
    return v;
}
long Params::inum(const std::string& k, long d) const {
    if (!has(k)) return d;
    const double v = num(k, (double)d);
    if (!(std::fabs(v) < 2e9)) throw std::invalid_argument("parameter " + k + ": not a whole number in range: " + kv.at(k));
    return (long)v;
}

std::string trim(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && isspace((unsigned char)s[a])) a++;
    while (b > a && isspace((unsigned char)s[b - 1])) b--;
    return s.substr(a, b - a);
}

Params parse_params(const char* text) {
    Params p;
    if (!text) return p;
    std::string s = text;
    size_t pos = 0;
    while (pos <= s.size()) {
        size_t nl = s.find_first_of("\n;", pos);
        if (nl == std::string::npos) nl = s.size();
        std::string line = trim(s.substr(pos, nl - pos));
        pos = nl + 1;
        if (line.empty() || line[0] == '#') continue;
        size_t eq = line.find_first_of("=: \t");
        if (eq == std::string::npos) { p.kv[line] = ""; continue; }
        std::string k = trim(line.substr(0, eq)), v = trim(line.substr(eq + 1));
        if (!v.empty() && (v[0] == '=' || v[0] == ':')) v = trim(v.substr(1));
        if (v == "None" || (v == "none" && k != "binning")) v = "";
        p.kv[k] = v;
    }
    return p;
}

// ============================================================ small helpers shared with isodec.cpp
// Python repr() of a float (shortest round trip; ".0" on whole numbers; exponent below 1e-4 / from 1e16)
std::string py_float(double v) {
    if (v != v) return "nan";
    if (!std::isfinite(v)) return v > 0 ? "inf" : "-inf";
    if (v == 0) return std::signbit(v) ? "-0.0" : "0.0";
    char buf[64];
    int p = 1;
    for (; p <= 17; p++) {
        snprintf(buf, sizeof buf, "%.*e", p - 1, v);
        if (strtod(buf, nullptr) == v) break;
    }
    if (p > 17) { p = 17; snprintf(buf, sizeof buf, "%.16e", v); }
    // buf: [-]d.ddde[+-]XX
    std::string s = buf;
    bool neg = s[0] == '-';
    if (neg) s = s.substr(1);
    size_t e = s.find('e');
    std::string mant = s.substr(0, e);
    int ex = atoi(s.c_str() + e + 1);
    std::string digits;
    for (char c : mant) if (c != '.') digits += c;
    while (digits.size() > 1 && digits.back() == '0') digits.pop_back();
    std::string out;
    if (ex >= -4 && ex < 16) {
        if (ex >= 0) {
            std::string ip = digits.substr(0, std::min<size_t>(digits.size(), ex + 1));
            while ((int)ip.size() < ex + 1) ip += '0';
            std::string fp = digits.size() > (size_t)ex + 1 ? digits.substr(ex + 1) : "";
            out = ip + "." + (fp.empty() ? "0" : fp);
        } else {
            out = "0." + std::string(-ex - 1, '0') + digits;
        }
    } else {
        out = digits.substr(0, 1);
        if (digits.size() > 1) out += "." + digits.substr(1);
        char eb[16];
        snprintf(eb, sizeof eb, "e%c%02d", ex < 0 ? '-' : '+', std::abs(ex));
        out += eb;
    }
    return neg ? "-" + out : out;
}

// ---- MD5 (RFC 1321), for the shortened run names
std::string md5_hex(const std::string& msg) {
    static const uint32_t K[64] = {
        0xd76aa478, 0xe8c7b756, 0x242070db, 0xc1bdceee, 0xf57c0faf, 0x4787c62a, 0xa8304613, 0xfd469501,
        0x698098d8, 0x8b44f7af, 0xffff5bb1, 0x895cd7be, 0x6b901122, 0xfd987193, 0xa679438e, 0x49b40821,
        0xf61e2562, 0xc040b340, 0x265e5a51, 0xe9b6c7aa, 0xd62f105d, 0x02441453, 0xd8a1e681, 0xe7d3fbc8,
        0x21e1cde6, 0xc33707d6, 0xf4d50d87, 0x455a14ed, 0xa9e3e905, 0xfcefa3f8, 0x676f02d9, 0x8d2a4c8a,
        0xfffa3942, 0x8771f681, 0x6d9d6122, 0xfde5380c, 0xa4beea44, 0x4bdecfa9, 0xf6bb4b60, 0xbebfbc70,
        0x289b7ec6, 0xeaa127fa, 0xd4ef3085, 0x04881d05, 0xd9d4d039, 0xe6db99e5, 0x1fa27cf8, 0xc4ac5665,
        0xf4292244, 0x432aff97, 0xab9423a7, 0xfc93a039, 0x655b59c3, 0x8f0ccc92, 0xffeff47d, 0x85845dd1,
        0x6fa87e4f, 0xfe2ce6e0, 0xa3014314, 0x4e0811a1, 0xf7537e82, 0xbd3af235, 0x2ad7d2bb, 0xeb86d391};
    static const int R[64] = {7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22,
                              5, 9, 14, 20, 5, 9, 14, 20, 5, 9, 14, 20, 5, 9, 14, 20,
                              4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23,
                              6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21};
    uint32_t h0 = 0x67452301, h1 = 0xefcdab89, h2 = 0x98badcfe, h3 = 0x10325476;
    std::vector<uint8_t> m(msg.begin(), msg.end());
    uint64_t bits = (uint64_t)m.size() * 8;
    m.push_back(0x80);
    while (m.size() % 64 != 56) m.push_back(0);
    for (int i = 0; i < 8; i++) m.push_back((uint8_t)(bits >> (8 * i)));
    for (size_t off = 0; off < m.size(); off += 64) {
        uint32_t w[16];
        for (int i = 0; i < 16; i++)
            w[i] = (uint32_t)m[off + 4 * i] | ((uint32_t)m[off + 4 * i + 1] << 8) | ((uint32_t)m[off + 4 * i + 2] << 16) | ((uint32_t)m[off + 4 * i + 3] << 24);
        uint32_t a = h0, b = h1, c = h2, d = h3;
        for (int i = 0; i < 64; i++) {
            uint32_t f; int g;
            if (i < 16) { f = (b & c) | (~b & d); g = i; }
            else if (i < 32) { f = (d & b) | (~d & c); g = (5 * i + 1) % 16; }
            else if (i < 48) { f = b ^ c ^ d; g = (3 * i + 5) % 16; }
            else { f = c ^ (b | ~d); g = (7 * i) % 16; }
            uint32_t t = d; d = c; c = b;
            uint32_t x = a + f + K[i] + w[g];
            b = b + ((x << R[i]) | (x >> (32 - R[i])));
            a = t;
        }
        h0 += a; h1 += b; h2 += c; h3 += d;
    }
    char out[33];
    uint32_t hs[4] = {h0, h1, h2, h3};
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++) snprintf(out + 8 * i + 2 * j, 3, "%02x", (hs[i] >> (8 * j)) & 0xff);
    return std::string(out, 32);
}

// ---- file system
std::string abs_path(const std::string& p) {
#ifdef _WIN32
    char buf[4096];
    DWORD n = GetFullPathNameA(p.c_str(), sizeof buf, buf, nullptr);
    if (n == 0 || n >= sizeof buf) return p;
    return std::string(buf, n);
#else
    if (!p.empty() && p[0] == '/') return p;
    char cwd[PATH_MAX];
    if (!getcwd(cwd, sizeof cwd)) return p;
    return std::string(cwd) + "/" + p;
#endif
}

std::string join_path(const std::string& a, const std::string& b) {
    if (a.empty()) return b;
    char last = a.back();
#ifdef _WIN32
    if (last == '\\' || last == '/') return a + b;
    return a + "\\" + b;
#else
    if (last == '/') return a + b;
    return a + "/" + b;
#endif
}

bool is_dir(const std::string& p) {
    struct stat st;
    return stat(p.c_str(), &st) == 0 && (st.st_mode & S_IFDIR);
}

void make_dirs(const std::string& path) {
    if (path.empty() || is_dir(path)) return;
    size_t cut = path.find_last_of("/\\");
    if (cut != std::string::npos && cut > 0) {
        std::string parent = path.substr(0, cut);
        if (!(parent.size() == 2 && parent[1] == ':')) make_dirs(parent);
    }
#ifdef _WIN32
    _mkdir(path.c_str());
#else
    mkdir(path.c_str(), 0775);
#endif
    if (!is_dir(path)) throw std::runtime_error("cannot create the folder " + path);
}

std::string temp_dir() {
#ifdef _WIN32
    char buf[MAX_PATH + 1];
    DWORD n = GetTempPathA(sizeof buf, buf);
    std::string t = n ? std::string(buf, n) : std::string("C:\\Temp\\");
    while (!t.empty() && (t.back() == '\\' || t.back() == '/')) t.pop_back();
    return t;
#else
    const char* t = getenv("TMPDIR");
    return t && *t ? t : "/tmp";
#endif
}

// Windows path limit of UniDec (ms_deconv._unidec_place): the name as it is
// when it fits, else shortened with a short hash, else the temporary folder
void unidec_place(std::string& folder, std::string& name) {
    const int limit = 240;
    auto longest = [](const std::string& fo, const std::string& na) {
        return (int)join_path(join_path(abs_path(fo), na + "_unidecfiles"), na + "_manualfile.dat").size();
    };
    if (longest(folder, name) <= limit) return;
    std::string tag = md5_hex(name).substr(0, 6);
    for (int n = (int)name.size() - 1; n > 7; n--) {
        std::string s = name.substr(0, n);
        while (!s.empty() && strchr(" ._-", s.back())) s.pop_back();
        s += "_" + tag;
        if (longest(folder, s) <= limit) { name = s; return; }
    }
    folder = join_path(temp_dir(), "MS_Analysis_UniDec");
    name = "run_" + tag;
}

// ---- numpy helpers
// np.arange(start, stop, step): numpy fills start + i * delta with delta = (start + step) - start
std::vector<double> np_arange(double start, double stop, double step) {
    double len = std::ceil((stop - start) / step);
    if (len > 2e9) throw std::invalid_argument("too many points for this range and step");
    size_t n = len > 0 ? (size_t)len : 0;
    std::vector<double> v(n);
    double delta = (start + step) - start;
    for (size_t i = 0; i < n; i++) v[i] = start + (double)i * delta;
    return v;
}

// scipy.ndimage.gaussian_filter1d(v, sigma) with the defaults (truncate 4, mode reflect)
void gaussian_filter_reflect(std::vector<double>& v, double sigma) {
    long n = (long)v.size();
    if (n == 0 || !(sigma > 0)) return;
    if (sigma > 1e6) throw std::invalid_argument("smoothing width too large");
    int r = (int)(4.0 * sigma + 0.5);
    std::vector<double> w(2 * r + 1);
    double s2 = sigma * sigma, sum = 0;
    for (int i = -r; i <= r; i++) { w[i + r] = std::exp(-0.5 * (double)(i * i) / s2); sum += w[i + r]; }
    for (double& x : w) x /= sum;
    std::vector<double> out(n);
    for (long i = 0; i < n; i++) {
        double s = 0;
        for (int t = -r; t <= r; t++) {
            long j = i + t;
            while (j < 0 || j >= n) { if (j < 0) j = -j - 1; if (j >= n) j = 2 * n - 1 - j; }
            s += w[t + r] * v[j];
        }
        out[i] = s;
    }
    v.swap(out);
}

}  // namespace unidec
}  // namespace ms

namespace {

using namespace ms::unidec;

// ============================================================ UniDec configuration (unidecstructure.UniDecConfig)
struct UConfig {
    std::string version = "8.2.1";
    int imflag = 0, cdmsflag = 0;
    std::string infname, outfname;
    long numit = 100; int numz = 50, endz = 50, startz = 1;
    double zzsig = 1, psig = 1, beta = 0;
    int suppression_topn = 0; double suppression_topx = 0, suppression_percent = 0; int suppression_startit = 3, suppression_harmonic = 0, suppression_satellite = 0;
    double mzsig = 0.85; int psfun = 0, psfunz = 0, discreteplot = 0, smashflag = 0;
    double massub = 500000, masslb = 5000, msig = 0, molig = 0, massbins = 10, mtabsig = 0;
    int isomode = 0;
    double minmz = 0, maxmz = 0, subbuff = 0; int subtype = 2; double smooth = 0, mzbins = 0, peakwindow = 500, peakthresh = 0.1;
    int normthresh = 1; double peakplotthresh = 0.1, separation = 0.025, intthresh = 0, reductionpercent = 0;
    int aggressiveflag = 0, rawflag = 0; double adductmass = 1.007276467, nativezub = 1000, nativezlb = -1000;
    int poolflag = 2; double detectoreffva = 0, inflate = 1; int noiseflag = 0, linflag = 2;
    std::string cmap = "nipy_spectral", peakcmap = "rainbow", spectracmap = "rainbow";
    int publicationmode = 1, isotopemode = 0, peaknorm = 1, datanorm = 1, baselineflag = 1, orbimode = 0;
    int filterwidth = 20; double zerolog = -12;
    bool doubledec = false;
    double CDslope = 0.2074; int CDzbins = 1; double CDres = 0, CDScanStart = -1, CDScanEnd = -1, CDprethresh = 0, HTksmooth = 0;
    int CDScanCompress = 0, HTmaxscans = -1; double HToutputlb = -1, HToutputub = -1, HTtimepad = 0, HTanalysistime = 38.0;
    std::string HTxaxis = "Time"; int HTcycleindex = -1; double HTcycletime = 2.0, HTtimeshift = 7.0; std::string htbit = "5";
    double FTstart = 5, FTend = 1000; bool FTflatten = true; int FTapodize = 1; std::string demultiplexmode = "HT"; double FTsmooth = 0;
    int HTmaskn = 1000, HTwin = 5; double csig = 0, smoothdt = 0, subbufdt = 0; bool centroided = false, verbose = false;
    int zout = 0, pusher = 0; double ccsub = 25000, ccslb = 100, dtsig = 0.2, ccsbins = 100, nativeccsub = 20000, nativeccslb = -20000;
    int twaveflag = 0; double temp = 25, pressure = 2, volt = 50, gasmass = 4.002602, to = 0, driftlength = 0.18202;
    double tcal1 = 0.3293, tcal2 = 6.3597, tcal3 = 0, tcal4 = 0, edc = 1.57; int variablepw = 0; double minratio = 0;
    double error = 0;   // r squared after the run
};

std::string ival(long v) { return std::to_string(v); }
std::string bval(bool v) { return v ? "True" : "False"; }

// unidecstructure.config_export: "key value" lines in UniDec's order (minus the
// file keys it skips), then the IsoDec configuration it appends
std::string conf_text(const UConfig& c) {
    std::string s;
    auto add = [&](const char* k, const std::string& v) { if (!v.empty() && v != " ") { s += k; s += " "; s += v; s += "\n"; } };
    add("version", c.version); add("imflag", ival(c.imflag)); add("cdmsflag", ival(c.cdmsflag));
    add("input", c.infname); add("output", c.outfname); add("numit", ival(c.numit)); add("numz", ival(c.numz));
    add("endz", ival(c.endz)); add("startz", ival(c.startz)); add("zzsig", py_float(c.zzsig)); add("psig", py_float(c.psig));
    add("beta", py_float(c.beta)); add("suppression_topn", ival(c.suppression_topn)); add("suppression_topx", py_float(c.suppression_topx));
    add("suppression_percent", py_float(c.suppression_percent)); add("suppression_startit", ival(c.suppression_startit));
    add("suppression_harmonic", ival(c.suppression_harmonic)); add("suppression_satellite", ival(c.suppression_satellite));
    add("mzsig", py_float(c.mzsig)); add("psfun", ival(c.psfun)); add("zpsfn", ival(c.psfunz)); add("discreteplot", ival(c.discreteplot));
    add("smashflag", ival(c.smashflag)); add("massub", py_float(c.massub)); add("masslb", py_float(c.masslb)); add("msig", py_float(c.msig));
    add("molig", py_float(c.molig)); add("massbins", py_float(c.massbins)); add("mtabsig", py_float(c.mtabsig)); add("isomode", ival(c.isomode));
    add("minmz", py_float(c.minmz)); add("maxmz", py_float(c.maxmz)); add("subbuff", py_float(c.subbuff)); add("subtype", ival(c.subtype));
    add("smooth", py_float(c.smooth)); add("mzbins", py_float(c.mzbins)); add("peakwindow", py_float(c.peakwindow));
    add("peakthresh", py_float(c.peakthresh)); add("normthresh", ival(c.normthresh)); add("peakplotthresh", py_float(c.peakplotthresh));
    add("plotsep", py_float(c.separation)); add("intthresh", py_float(c.intthresh)); add("reductionpercent", py_float(c.reductionpercent));
    add("aggressive", ival(c.aggressiveflag)); add("rawflag", ival(c.rawflag)); add("adductmass", py_float(c.adductmass));
    add("nativezub", py_float(c.nativezub)); add("nativezlb", py_float(c.nativezlb)); add("poolflag", ival(c.poolflag));
    add("accvol", py_float(c.detectoreffva)); add("peakshapeinflate", py_float(c.inflate)); add("noiseflag", ival(c.noiseflag));
    add("linflag", ival(c.linflag)); add("cmap", c.cmap); add("peakcmap", c.peakcmap); add("spectracmap", c.spectracmap);
    add("publicationmode", ival(c.publicationmode)); add("isotopemode", ival(c.isotopemode)); add("peaknorm", ival(c.peaknorm));
    add("datanorm", ival(c.datanorm)); add("baselineflag", ival(c.baselineflag)); add("orbimode", ival(c.orbimode));
    add("integratelb", "-1"); add("integrateub", "-1"); add("filterwidth", ival(c.filterwidth)); add("zerolog", py_float(c.zerolog));
    add("doubledec", bval(c.doubledec)); add("CDslope", py_float(c.CDslope)); add("CDzbins", ival(c.CDzbins)); add("CDres", py_float(c.CDres));
    add("CDScanStart", py_float(c.CDScanStart)); add("CDScanEnd", py_float(c.CDScanEnd)); add("CDprethresh", py_float(c.CDprethresh));
    add("HTksmooth", py_float(c.HTksmooth)); add("CDScanCompress", ival(c.CDScanCompress)); add("HTmaxscans", ival(c.HTmaxscans));
    add("HToutlb", py_float(c.HToutputlb)); add("HToutub", py_float(c.HToutputub)); add("HTtimepad", py_float(c.HTtimepad));
    add("HTanalysistime", py_float(c.HTanalysistime)); add("HTxaxis", c.HTxaxis); add("HTcycleindex", ival(c.HTcycleindex));
    add("HTcylctime", py_float(c.HTcycletime)); add("HTtimeshift", py_float(c.HTtimeshift)); add("htbit", c.htbit);
    add("FTstart", py_float(c.FTstart)); add("FTend", py_float(c.FTend)); add("FTflatten", bval(c.FTflatten)); add("FTapodize", ival(c.FTapodize));
    add("demultiplexmode", c.demultiplexmode); add("FTsmooth", py_float(c.FTsmooth)); add("HTmaskn", ival(c.HTmaskn)); add("HTwin", ival(c.HTwin));
    add("csig", py_float(c.csig)); add("smoothdt", py_float(c.smoothdt)); add("subbufdt", py_float(c.subbufdt));
    add("centroided", ival(c.centroided ? 1 : 0)); add("verbose", ival(c.verbose ? 1 : 0)); add("zout", ival(c.zout)); add("pusher", ival(c.pusher));
    add("mindt", "-1"); add("maxdt", "-1"); add("ccsub", py_float(c.ccsub)); add("ccslb", py_float(c.ccslb)); add("dtsig", py_float(c.dtsig));
    add("ccsbins", py_float(c.ccsbins)); add("ubnativeccs", py_float(c.nativeccsub)); add("lbnativeccs", py_float(c.nativeccslb));
    add("twaveflag", ival(c.twaveflag)); add("temp", py_float(c.temp)); add("pressure", py_float(c.pressure)); add("volt", py_float(c.volt));
    add("gasmass", py_float(c.gasmass)); add("tnaught", py_float(c.to)); add("driftlength", py_float(c.driftlength));
    add("tcal1", py_float(c.tcal1)); add("tcal2", py_float(c.tcal2)); add("tcal3", py_float(c.tcal3)); add("tcal4", py_float(c.tcal4));
    add("edc", py_float(c.edc)); add("variablepw", ival(c.variablepw)); add("minratio", py_float(c.minratio));
    // IsoDecConfig.config_export (defaults; the engine ignores these)
    s += "iso_meanpeakspacing_thresh 0.01\niso_background_subtraction 0\niso_adduct_mass 1.007276467\niso_mass_diff_c 1.0033\n"
         "iso_peak_window 80\niso_phaseres 8\niso_matchtol 5\niso_minpeaks 3\niso_peak_thresh 0.0001\niso_css_thresh 0.7\n"
         "iso_maxshift 3\niso_mzwindowlb -1.05\niso_mzwindowub 4.05\niso_plusoneintwindowlb 0.1\niso_plusoneintwindowub 0.6\n"
         "iso_knockdown_rounds 5\niso_min_score_diff 0.1\niso_minareacovered 0.2\niso_minusoneaszero 1\niso_isotopethreshold 0.01\n"
         "iso_datathreshold 0.05\niso_zscore_threshold 0.95\niso_report_multiple_monoisos 1\niso_write_scans_without_precs 1\n"
         "iso_write_msalign 0\niso_write_tsv 1\n";
    return s;
}

void write_text(const std::string& path, const std::string& text) {
    FILE* f = fopen(path.c_str(), "w");
    if (!f) throw std::runtime_error("cannot write " + path);
    fputs(text.c_str(), f);
    fclose(f);
}

// the engine's configuration from the text of conf.dat, read as unidec.exe reads it (LoadConfig:
// "key value" lines, the keys matched as the engine matches them, then PostImport)
ms::udengine::Config engine_config(const std::string& conf) {
    ms::udengine::Config ec;
    size_t pos = 0;
    while (pos < conf.size()) {
        size_t nl = conf.find('\n', pos);
        if (nl == std::string::npos) nl = conf.size();
        const std::string line = conf.substr(pos, nl - pos);
        pos = nl + 1;
        const size_t sp = line.find(' ');
        if (sp == std::string::npos) continue;
        const std::string key = line.substr(0, sp);
        size_t v = sp;
        while (v < line.size() && line[v] == ' ') v++;
        if (key == "input" || key == "output") continue;   // file names: not used in the library
        ec.load_line(key, line.substr(v));
    }
    ec.post_import();
    return ec;
}

// m/z as unidec.exe read it: strtof of the "%.17g" text of the input file. (float)x is the same
// except where x lies within the 17 digit rounding of the midpoint between two floats.
float engine_mz(double x) {
    const float f = (float)x;
    if (!std::isfinite(f) || (double)f == x) return f;
    const float g = std::nextafter(f, x > (double)f ? std::numeric_limits<float>::infinity() : -std::numeric_limits<float>::infinity());
    const double mid = 0.5 * ((double)f + (double)g);
    if (std::fabs(x - mid) > 1e-15 * std::fabs(x)) return f;
    char b[40];
    snprintf(b, sizeof b, "%.17g", x);
    return strtof(b, nullptr);
}

// intensity as unidec.exe read it: strtof of the "%.6f" text of the input file
float engine_intensity(double y) {
    char b[64];
    snprintf(b, sizeof b, "%.6f", y);
    return strtof(b, nullptr);
}

// ============================================================ tools.dataprep
// tools.nearest: bisect_left, then the closer neighbour (with its edge rules)
long ud_nearest(const std::vector<double>& a, double t) {
    long n = (long)a.size();
    long i = std::lower_bound(a.begin(), a.end(), t) - a.begin();
    if (i <= 0) return 0;
    if (i >= n - 1) return n - 1;
    if (std::fabs(a[i] - t) > std::fabs(a[i - 1] - t)) i--;
    return i;
}

std::vector<double> nonlinear_axis(double start, double end, double res) {
    std::vector<double> axis;
    double i = start;
    axis.push_back(i);
    i += i / res;
    while (i < end) { axis.push_back(i); i += i / res; }
    return axis;
}

// tools.linearize with lintegrate (linflag 0 linear bins, 1 constant resolution)
Spec ud_linearize(const Spec& d, double binsize, int linflag) {
    size_t n = d.x.size();
    double firstpoint = std::ceil(d.x[0] / binsize) * binsize;
    double lastpoint = std::floor(d.x[n - 1] / binsize) * binsize;
    std::vector<double> intx = (linflag == 0 || linflag == 3) ? np_arange(firstpoint, lastpoint, binsize)
                                                               : nonlinear_axis(firstpoint, lastpoint, firstpoint / binsize);
    Spec out;
    out.x = intx;
    out.y.assign(intx.size(), 0.0);
    long l2 = (long)intx.size();
    if (l2 == 0) return out;
    for (size_t i = 0; i < n; i++) {
        double x = d.x[i], y = d.y[i];
        if (!(intx[0] < x && x < intx[l2 - 1])) continue;
        long index = ud_nearest(intx, x);
        if (intx[index] == x) out.y[index] += y;
        else if (intx[index] < x && index < l2 - 1) {
            long index2 = index + 1;
            double interpos = (x - intx[index]) / (intx[index2] - intx[index]);
            out.y[index] += (1 - interpos) * y;
            out.y[index2] += interpos * y;
        } else if (intx[index] > x && index > 0) {
            long index2 = index - 1;
            double interpos = (x - intx[index]) / (intx[index2] - intx[index]);
            out.y[index] += (1 - interpos) * y;
            out.y[index2] += interpos * y;
        }
    }
    return out;
}

// ms_deconv._engine_unique: points indistinguishable to the engine's float32
// reader merged (mean m/z, summed intensity)
Spec engine_unique(const Spec& d) {
    size_t n = d.x.size();
    if (n == 0) return d;
    std::vector<float> keys(n);
    bool increasing = true;
    for (size_t i = 0; i < n; i++) {
        keys[i] = (float)d.x[i];
        if (!std::isfinite(keys[i])) throw std::invalid_argument("m/z values exceed the UniDec engine's numeric range");
        if (i > 0 && !(keys[i] > keys[i - 1])) increasing = false;
    }
    if (increasing) return d;
    // runs of equal float32 keys (the data are sorted)
    Spec out;
    size_t i = 0;
    while (i < n) {
        size_t j = i;
        double sx = 0, sy = 0;
        while (j < n && keys[j] == keys[i]) { sx += d.x[j]; sy += d.y[j]; j++; }
        out.x.push_back(sx / (double)(j - i));
        out.y.push_back(sy);
        i = j;
    }
    return out;
}

// tools.dataprep for 1D MS data (crop, smooth, binning, background, normalise, threshold 0, duplicates)
Spec ud_dataprep(const Spec& raw, const UConfig& c) {
    Spec d;
    for (size_t i = 0; i < raw.x.size(); i++)
        if (raw.x[i] <= c.maxmz && raw.x[i] >= c.minmz) { d.x.push_back(raw.x[i]); d.y.push_back(raw.y[i]); }
    if (d.x.empty()) throw std::invalid_argument("No data left to deconvolute in the m/z range (check the m/z range and the minimum intensity)");
    if (c.smooth > 0) gaussian_filter_reflect(d.y, c.smooth);
    if (c.mzbins > 0 && c.linflag != 2) d = ud_linearize(d, c.mzbins, c.linflag);
    size_t n = d.x.size();
    double buff = std::fabs(c.subbuff);
    if (buff != 0 && n > 0) {
        if (c.subtype == 1) {           // datasimpsub: line between the averaged ends
            long b = (long)buff;
            double front = 0, back = 0;
            long nf = std::min<long>(b, (long)n);
            for (long i = 0; i < nf; i++) front += d.y[i];
            front = nf ? front / nf : 0;
            long b0 = std::max<long>(0, (long)n - b - 1), b1 = (long)n - 1;
            for (long i = b0; i < b1; i++) back += d.y[i];
            back = b1 > b0 ? back / (b1 - b0) : 0;
            for (size_t i = 0; i < n; i++) d.y[i] -= front + (back - front) / (double)n * (double)i;
        } else if (c.subtype == 2) {    // datacompsub: windowed minimum smoothed by a Gaussian of 2 x buff
            std::vector<double> mins(n);
            for (long i = 0; i < (long)n; i++) {
                long a = std::max(0L, i - (long)buff), b = std::min((long)n, i + (long)buff);
                double m = std::numeric_limits<double>::infinity();
                for (long j = a; j < b; j++) m = std::min(m, d.y[j]);
                mins[i] = b > a ? m : 0.0;
            }
            gaussian_filter_reflect(mins, buff * 2);
            for (size_t i = 0; i < n; i++) d.y[i] -= mins[i];
        } else if (c.subtype == 0) {
            double m = *std::min_element(d.y.begin(), d.y.end());
            for (double& v : d.y) v -= m;
        } else if (c.subtype == 3) {
            double m = *std::max_element(d.y.begin(), d.y.end());
            for (double& v : d.y) v = std::max(0.0, v - buff * m);
        } else {
            throw std::invalid_argument("UniDec background subtraction type " + std::to_string(c.subtype) + " is not supported");
        }
    }
    if (c.datanorm == 1 && n) {
        double mx = *std::max_element(d.y.begin(), d.y.end());
        if (mx != 0) for (double& v : d.y) v /= mx;
    }
    double thr = c.intthresh > 0 ? c.intthresh : 0.0;
    for (double& v : d.y) if (v < thr) v = 0.0;
    if (c.linflag == 2 && n > 2) {      // remove_middle_zeros
        Spec k;
        for (size_t i = 0; i < n; i++) {
            bool keep = i == 0 || i == n - 1 || d.y[i] != 0 || d.y[i - 1] != 0 || d.y[i + 1] != 0;
            if (keep) { k.x.push_back(d.x[i]); k.y.push_back(d.y[i]); }
        }
        d = k;
    }
    return engine_unique(d);
}

// ============================================================ tools.peakdetect
struct XY { double x, y; };

// local maxima within +- window points (a float, as peakwindow / massbins) above threshold x the top
std::vector<XY> ud_peakdetect(const double* x, const double* y, long n, double window, double threshold, bool norm = true) {
    std::vector<XY> peaks;
    if (n == 0) return peaks;
    double maxval = norm ? *std::max_element(y, y + n) : 1.0;
    for (long i = 0; i < n; i++) {
        if (!(y[i] > maxval * threshold)) continue;
        long start = (long)((double)i - window);      // int(): towards zero
        long end = (long)((double)i + window) + 1;
        if (start < 0) start = 0;
        if (end > n) end = n;
        double testmax = y[start];
        for (long j = start + 1; j < end; j++) testmax = std::max(testmax, y[j]);
        bool first = true;
        for (long j = start; j < i; j++) if (y[j] == y[i]) { first = false; break; }
        if (y[i] == testmax && first) peaks.push_back({x[i], y[i]});
    }
    return peaks;
}

// ============================================================ auto peak width (tools.auto_peak_width)
// radix 2 FFT (in place), for the autocorrelation
void fft(std::vector<std::complex<double>>& a, bool inverse) {
    size_t n = a.size();
    for (size_t i = 1, j = 0; i < n; i++) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    for (size_t len = 2; len <= n; len <<= 1) {
        double ang = 2 * M_PI / (double)len * (inverse ? 1 : -1);
        std::complex<double> wl(std::cos(ang), std::sin(ang));
        for (size_t i = 0; i < n; i += len) {
            std::complex<double> w(1);
            for (size_t j = 0; j < len / 2; j++) {
                std::complex<double> u = a[i + j], v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
                w *= wl;
            }
        }
    }
    if (inverse) for (auto& v : a) v /= (double)n;
}

// autocorrelation R(lag) for lag 0..n-1
std::vector<double> autocorrelation(const std::vector<double>& y) {
    size_t n = y.size(), m = 1;
    while (m < 2 * n) m <<= 1;
    std::vector<std::complex<double>> a(m);
    for (size_t i = 0; i < n; i++) a[i] = y[i];
    fft(a, false);
    for (auto& v : a) v = v * std::conj(v);
    fft(a, true);
    std::vector<double> r(n);
    for (size_t i = 0; i < n; i++) r[i] = a[i].real();
    return r;
}

// peak shape functions of fitting.py with area normalisation and a background (psfit)
double psfun_value(int psfun, double x, double s, double m, double a, double b) {
    if (b < 0 || a < 0 || s < 0) return 0.0;
    if (psfun == 0) {
        double sig = s / 2.35482;
        double aa = a / (sig * std::sqrt(2 * M_PI));
        return aa * std::exp(-(x - m) * (x - m) / (2.0 * sig * sig)) + b;
    }
    if (psfun == 1) {
        double aa = a * (1 / M_PI) * (s / 2.0);
        return aa / ((x - m) * (x - m) + (s / 2.0) * (s / 2.0)) + b;
    }
    double sig2 = s / (2 * std::sqrt(2 * std::log(2.0)));
    double a1 = a * ((1 / M_PI) / (s / 2.0)) / 0.83723895067;
    double a2 = a * 2.0 / (s * M_PI) / 0.83723895067;
    if (m < x) return a1 * ((s / 2.0) * (s / 2.0)) / ((x - m) * (x - m) + (s / 2.0) * (s / 2.0)) + b;
    return a2 * std::exp(-(x - m) * (x - m) / (2.0 * sig2 * sig2)) + b;
}

// ---- MINPACK lmdif (what scipy.optimize.curve_fit / leastsq runs), ported so that the
// peak width fit follows the same path as the Python: forward difference Jacobian,
// QR with column pivoting, the Levenberg Marquardt parameter of lmpar, trust region
// updates and the ftol / xtol / gtol / maxfev stopping rules.
namespace minpack {

double enorm(const double* v, int n) {
    double s = 0;
    for (int i = 0; i < n; i++) s += v[i] * v[i];
    return std::sqrt(s);
}

// a is m x n, column major (a[i + j * m])
void qrfac(int m, int n, std::vector<double>& a, std::vector<int>& ipvt, std::vector<double>& rdiag,
           std::vector<double>& acnorm, std::vector<double>& wa) {
    const double epsmch = DBL_EPSILON;
    for (int j = 0; j < n; j++) {
        acnorm[j] = enorm(&a[j * m], m);
        rdiag[j] = acnorm[j];
        wa[j] = rdiag[j];
        ipvt[j] = j;
    }
    int minmn = std::min(m, n);
    for (int j = 0; j < minmn; j++) {
        int kmax = j;
        for (int k = j + 1; k < n; k++) if (rdiag[k] > rdiag[kmax]) kmax = k;
        if (kmax != j) {
            for (int i = 0; i < m; i++) std::swap(a[i + j * m], a[i + kmax * m]);
            rdiag[kmax] = rdiag[j];
            wa[kmax] = wa[j];
            std::swap(ipvt[j], ipvt[kmax]);
        }
        double ajnorm = enorm(&a[j + j * m], m - j);
        if (ajnorm != 0) {
            if (a[j + j * m] < 0) ajnorm = -ajnorm;
            for (int i = j; i < m; i++) a[i + j * m] /= ajnorm;
            a[j + j * m] += 1;
            for (int k = j + 1; k < n; k++) {
                double sum = 0;
                for (int i = j; i < m; i++) sum += a[i + j * m] * a[i + k * m];
                double temp = sum / a[j + j * m];
                for (int i = j; i < m; i++) a[i + k * m] -= temp * a[i + j * m];
                if (rdiag[k] != 0) {
                    temp = a[j + k * m] / rdiag[k];
                    rdiag[k] *= std::sqrt(std::max(0.0, 1 - temp * temp));
                    double q = rdiag[k] / wa[k];
                    if (0.05 * q * q <= epsmch) {
                        rdiag[k] = enorm(&a[(j + 1) + k * m], m - j - 1);
                        wa[k] = rdiag[k];
                    }
                }
            }
        }
        rdiag[j] = -ajnorm;
    }
}

// r: n x n upper triangle in a column major m x n array (ldr = m)
void qrsolv(int n, std::vector<double>& r, int ldr, const std::vector<int>& ipvt, const std::vector<double>& diag,
            const std::vector<double>& qtb, std::vector<double>& x, std::vector<double>& sdiag, std::vector<double>& wa) {
    auto R = [&](int i, int j) -> double& { return r[i + j * ldr]; };
    for (int j = 0; j < n; j++) {
        for (int i = j; i < n; i++) R(i, j) = R(j, i);
        x[j] = R(j, j);
        wa[j] = qtb[j];
    }
    for (int j = 0; j < n; j++) {
        int l = ipvt[j];
        if (diag[l] != 0) {
            for (int k = j; k < n; k++) sdiag[k] = 0;
            sdiag[j] = diag[l];
            double qtbpj = 0;
            for (int k = j; k < n; k++) {
                if (sdiag[k] == 0) continue;
                double sn, cs;
                if (std::fabs(R(k, k)) < std::fabs(sdiag[k])) {
                    double cotan = R(k, k) / sdiag[k];
                    sn = 0.5 / std::sqrt(0.25 + 0.25 * cotan * cotan);
                    cs = sn * cotan;
                } else {
                    double tn = sdiag[k] / R(k, k);
                    cs = 0.5 / std::sqrt(0.25 + 0.25 * tn * tn);
                    sn = cs * tn;
                }
                R(k, k) = cs * R(k, k) + sn * sdiag[k];
                double temp = cs * wa[k] + sn * qtbpj;
                qtbpj = -sn * wa[k] + cs * qtbpj;
                wa[k] = temp;
                for (int i = k + 1; i < n; i++) {
                    temp = cs * R(i, k) + sn * sdiag[i];
                    sdiag[i] = -sn * R(i, k) + cs * sdiag[i];
                    R(i, k) = temp;
                }
            }
        }
        sdiag[j] = R(j, j);
        R(j, j) = x[j];
    }
    int nsing = n;
    for (int j = 0; j < n; j++) {
        if (sdiag[j] == 0 && nsing == n) nsing = j;
        if (nsing < n) wa[j] = 0;
    }
    for (int k = 0; k < nsing; k++) {
        int j = nsing - 1 - k;
        double sum = 0;
        for (int i = j + 1; i < nsing; i++) sum += R(i, j) * wa[i];
        wa[j] = (wa[j] - sum) / sdiag[j];
    }
    for (int j = 0; j < n; j++) x[ipvt[j]] = wa[j];
}

void lmpar(int n, std::vector<double>& r, int ldr, const std::vector<int>& ipvt, const std::vector<double>& diag,
           const std::vector<double>& qtb, double delta, double& par, std::vector<double>& x, std::vector<double>& sdiag) {
    const double dwarf = DBL_MIN;
    auto R = [&](int i, int j) -> double& { return r[i + j * ldr]; };
    std::vector<double> wa1(n), wa2(n);
    int nsing = n;
    for (int j = 0; j < n; j++) {
        wa1[j] = qtb[j];
        if (R(j, j) == 0 && nsing == n) nsing = j;
        if (nsing < n) wa1[j] = 0;
    }
    for (int k = 0; k < nsing; k++) {
        int j = nsing - 1 - k;
        wa1[j] /= R(j, j);
        double temp = wa1[j];
        for (int i = 0; i < j; i++) wa1[i] -= R(i, j) * temp;
    }
    for (int j = 0; j < n; j++) x[ipvt[j]] = wa1[j];
    int iter = 0;
    for (int j = 0; j < n; j++) wa2[j] = diag[j] * x[j];
    double dxnorm = enorm(wa2.data(), n);
    double fp = dxnorm - delta;
    if (fp <= 0.1 * delta) { par = 0; return; }
    double parl = 0;
    if (nsing >= n) {
        for (int j = 0; j < n; j++) { int l = ipvt[j]; wa1[j] = diag[l] * (wa2[l] / dxnorm); }
        for (int j = 0; j < n; j++) {
            double sum = 0;
            for (int i = 0; i < j; i++) sum += R(i, j) * wa1[i];
            wa1[j] = (wa1[j] - sum) / R(j, j);
        }
        double temp = enorm(wa1.data(), n);
        parl = ((fp / delta) / temp) / temp;
    }
    for (int j = 0; j < n; j++) {
        double sum = 0;
        for (int i = 0; i <= j; i++) sum += R(i, j) * qtb[i];
        wa1[j] = sum / diag[ipvt[j]];
    }
    double gnorm = enorm(wa1.data(), n);
    double paru = gnorm / delta;
    if (paru == 0) paru = dwarf / std::min(delta, 0.1);
    par = std::max(par, parl);
    par = std::min(par, paru);
    if (par == 0) par = gnorm / dxnorm;
    for (;;) {
        iter++;
        if (par == 0) par = std::max(dwarf, 0.001 * paru);
        double temp = std::sqrt(par);
        for (int j = 0; j < n; j++) wa1[j] = temp * diag[j];
        qrsolv(n, r, ldr, ipvt, wa1, qtb, x, sdiag, wa2);
        for (int j = 0; j < n; j++) wa2[j] = diag[j] * x[j];
        dxnorm = enorm(wa2.data(), n);
        temp = fp;
        fp = dxnorm - delta;
        if (std::fabs(fp) <= 0.1 * delta || (parl == 0 && fp <= temp && temp < 0) || iter == 10) break;
        for (int j = 0; j < n; j++) { int l = ipvt[j]; wa1[j] = diag[l] * (wa2[l] / dxnorm); }
        for (int j = 0; j < n; j++) {
            wa1[j] /= sdiag[j];
            temp = wa1[j];
            for (int i = j + 1; i < n; i++) wa1[i] -= R(i, j) * temp;
        }
        temp = enorm(wa1.data(), n);
        double parc = ((fp / delta) / temp) / temp;
        if (fp > 0) parl = std::max(parl, par);
        if (fp < 0) paru = std::min(paru, par);
        par = std::max(parl, par + parc);
    }
    if (iter == 0) par = 0;
}

// returns info (1..4 converged, 5 maxfev, 6..8 tolerances too small), fcn(x, fvec)
template <class F> int lmdif(F&& fcn, int m, int n, std::vector<double>& x, double ftol, double xtol, double gtol, int maxfev,
                             double factor = 100.0) {
    const double epsmch = DBL_EPSILON;
    std::vector<double> fvec(m), wa4(m), fjac((size_t)m * n), diag(n, 0.0), qtf(n), wa1(n), wa2(n), wa3(n), rdiag(n), acnorm(n);
    std::vector<int> ipvt(n);
    int nfev = 0, info = 0;
    fcn(x.data(), fvec.data()); nfev++;
    double fnorm = enorm(fvec.data(), m);
    double par = 0, xnorm = 0, delta = 0;
    int iter = 1;
    const double eps = std::sqrt(epsmch);
    for (;;) {
        // fdjac2
        for (int j = 0; j < n; j++) {
            double temp = x[j];
            double h = eps * std::fabs(temp);
            if (h == 0) h = eps;
            x[j] = temp + h;
            fcn(x.data(), wa4.data()); nfev++;
            x[j] = temp;
            for (int i = 0; i < m; i++) fjac[i + j * m] = (wa4[i] - fvec[i]) / h;
        }
        qrfac(m, n, fjac, ipvt, rdiag, acnorm, wa1);
        if (iter == 1) {
            for (int j = 0; j < n; j++) { diag[j] = acnorm[j]; if (acnorm[j] == 0) diag[j] = 1; }
            for (int j = 0; j < n; j++) wa3[j] = diag[j] * x[j];
            xnorm = enorm(wa3.data(), n);
            delta = factor * xnorm;
            if (delta == 0) delta = factor;
        }
        for (int i = 0; i < m; i++) wa4[i] = fvec[i];
        for (int j = 0; j < n; j++) {
            if (fjac[j + j * m] != 0) {
                double sum = 0;
                for (int i = j; i < m; i++) sum += fjac[i + j * m] * wa4[i];
                double temp = -sum / fjac[j + j * m];
                for (int i = j; i < m; i++) wa4[i] += fjac[i + j * m] * temp;
            }
            fjac[j + j * m] = rdiag[j];
            qtf[j] = wa4[j];
        }
        double gnorm = 0;
        if (fnorm != 0) {
            for (int j = 0; j < n; j++) {
                int l = ipvt[j];
                if (acnorm[l] == 0) continue;
                double sum = 0;
                for (int i = 0; i <= j; i++) sum += fjac[i + j * m] * (qtf[i] / fnorm);
                gnorm = std::max(gnorm, std::fabs(sum / acnorm[l]));
            }
        }
        if (gnorm <= gtol) return 4;
        for (int j = 0; j < n; j++) diag[j] = std::max(diag[j], acnorm[j]);
        for (;;) {
            lmpar(n, fjac, m, ipvt, diag, qtf, delta, par, wa1, wa2);
            for (int j = 0; j < n; j++) { wa1[j] = -wa1[j]; wa2[j] = x[j] + wa1[j]; wa3[j] = diag[j] * wa1[j]; }
            double pnorm = enorm(wa3.data(), n);
            if (iter == 1) delta = std::min(delta, pnorm);
            fcn(wa2.data(), wa4.data()); nfev++;
            double fnorm1 = enorm(wa4.data(), m);
            double actred = -1;
            if (0.1 * fnorm1 < fnorm) actred = 1 - (fnorm1 / fnorm) * (fnorm1 / fnorm);
            for (int j = 0; j < n; j++) {
                wa3[j] = 0;
                double temp = wa1[ipvt[j]];
                for (int i = 0; i <= j; i++) wa3[i] += fjac[i + j * m] * temp;
            }
            double temp1 = enorm(wa3.data(), n) / fnorm;
            double temp2 = (std::sqrt(par) * pnorm) / fnorm;
            double prered = temp1 * temp1 + temp2 * temp2 / 0.5;
            double dirder = -(temp1 * temp1 + temp2 * temp2);
            double ratio = prered != 0 ? actred / prered : 0;
            if (ratio <= 0.25) {
                double temp = actred >= 0 ? 0.5 : 0.5 * dirder / (dirder + 0.5 * actred);
                if (0.1 * fnorm1 >= fnorm || temp < 0.1) temp = 0.1;
                delta = temp * std::min(delta, pnorm / 0.1);
                par /= temp;
            } else if (par == 0 || ratio >= 0.75) {
                delta = pnorm / 0.5;
                par = 0.5 * par;
            }
            if (ratio >= 1e-4) {
                for (int j = 0; j < n; j++) { x[j] = wa2[j]; wa2[j] = diag[j] * x[j]; }
                for (int i = 0; i < m; i++) fvec[i] = wa4[i];
                xnorm = enorm(wa2.data(), n);
                fnorm = fnorm1;
                iter++;
            }
            if (std::fabs(actred) <= ftol && prered <= ftol && 0.5 * ratio <= 1) info = 1;
            if (delta <= xtol * xnorm) info = 2;
            if (std::fabs(actred) <= ftol && prered <= ftol && 0.5 * ratio <= 1 && info == 2) info = 3;
            if (info != 0) return info;
            if (nfev >= maxfev) info = 5;
            if (std::fabs(actred) <= epsmch && prered <= epsmch && 0.5 * ratio <= 1) info = 6;
            if (delta <= epsmch * xnorm) info = 7;
            if (gnorm <= epsmch) info = 8;
            if (info != 0) return info;
            if (ratio >= 1e-4) break;
        }
    }
}

}  // namespace minpack

// scipy curve_fit(psfun_fit, x, y, p0) with the defaults: p = [fwhm, mid, amplitude, background];
// false when the fit did not converge within maxfev (curve_fit raises RuntimeError)
bool lm_fit(int psfun, const std::vector<double>& x, const std::vector<double>& y, std::vector<double>& p) {
    int m = (int)x.size();
    auto fcn = [&](const double* q, double* f) {
        for (int i = 0; i < m; i++) f[i] = psfun_value(psfun, x[i], q[0], q[1], q[2], q[3]) - y[i];
    };
    int info = minpack::lmdif(fcn, m, 4, p, 1.49012e-8, 1.49012e-8, 0.0, 200 * (4 + 1));
    return info >= 1 && info <= 4;
}

struct PeakFit { std::vector<double> p; double error; bool ok; };

// fitting.isolated_peak_fit
PeakFit isolated_peak_fit(const std::vector<double>& x, const std::vector<double>& y, int psfun) {
    size_t n = x.size();
    size_t imax = std::max_element(y.begin(), y.end()) - y.begin();
    double midguess = x[imax];
    double bguess = *std::min_element(y.begin(), y.end());
    // weighted_std_2(x, y - b)
    double sw = 0, swx = 0;
    for (size_t i = 0; i < n; i++) { double w = y[i] - bguess; sw += w; swx += w * x[i]; }
    double avg = swx / sw, var = 0;
    for (size_t i = 0; i < n; i++) { double w = y[i] - bguess; var += w * (x[i] - avg) * (x[i] - avg); }
    double sigguess = std::sqrt(var / sw);
    auto maxfit = [&](double s, double m, double a, double b) {
        double mx = -std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < n; i++) mx = std::max(mx, psfun_value(psfun, x[i], s, m, a, b));
        return mx;
    };
    double ymax = y[imax];
    double aguess = ymax / maxfit(sigguess, midguess, 1, bguess);
    aguess = aguess * ymax / maxfit(sigguess, midguess, aguess, bguess);
    PeakFit f;
    f.p = {sigguess, midguess, aguess, bguess};
    f.ok = std::isfinite(sigguess) && std::isfinite(aguess) && lm_fit(psfun, x, y, f.p);
    f.error = 0;
    for (size_t i = 0; i < n; i++) {
        double d = y[i] - psfun_value(psfun, x[i], f.p[0], f.p[1], f.p[2], f.p[3]);
        f.error += d * d;
    }
    return f;
}

// tools.auto_peak_width: FWHM of the tallest peak, the width of its window from the
// first autocorrelation peak; the peak shape (0 Gaussian, 1 Lorentzian, 2 split) with
// the smallest error. Returns false when the fit failed (Python: exception).
bool auto_peak_width(const Spec& d, double& fwhm, int& psfun_out) {
    using ms::e2::np_sum;
    long n = (long)d.x.size();
    if (n < 3) return false;
    long maxpos1 = std::max_element(d.y.begin(), d.y.end()) - d.y.begin();
    double maxval = d.x[maxpos1];
    if (maxval == 0) { fwhm = 1; psfun_out = 0; return true; }
    std::vector<double> R = autocorrelation(d.y);
    double rmax = 0;
    for (double v : R) rmax = std::max(rmax, v);
    if (rmax == 0) { fwhm = 0; return true; }     // no autocorrelation peaks -> width 0
    for (double& v : R) v /= rmax;
    // x spacing: mean of the spacing around the tallest point
    long start = (long)std::max((double)maxpos1 - (double)n / 10.0, 0.0);
    long end = (long)std::min((double)(n - 1), (double)maxpos1 + (double)n / 10.0);
    long cs = start, ce = end;
    if (ce - cs < 20) { cs = 0; ce = n; }
    std::vector<double> diffs;
    for (long i = cs + 1; i < ce; i++) diffs.push_back(d.x[i] - d.x[i - 1]);
    if (diffs.empty()) return false;
    double xdiff = np_sum(diffs.data(), diffs.size()) / (double)diffs.size();
    // 'same' output of the full autocorrelation: index j <-> lag j - maxpos
    long maxpos = n - 1 - (n - 1) / 2;
    double c0 = (double)maxpos * xdiff;
    std::vector<double> ax, ay;
    // lags from 2 points on, by index: UniDec's test cx > xdiff let lag 1 through when the
    // rounding made cx a hair larger than xdiff
    for (long j = maxpos + 2; j < n; j++) {
        double cx = (double)j * xdiff - c0;
        ax.push_back(cx); ay.push_back(R[j - maxpos]);
    }
    // a small threshold: with 0, the round off of the FFT (1e-16 of the top) between the
    // peaks of a clean spectrum was taken for autocorrelation peaks
    const double autocorr_threshold = 1e-6;
    std::vector<XY> cpeaks = ud_peakdetect(ax.data(), ay.data(), (long)ax.size(), 10, autocorr_threshold, true);
    if (cpeaks.empty()) { fwhm = 0; return true; }
    auto select = [&](double sig, std::vector<double>& ix, std::vector<double>& iy) {
        ix.clear(); iy.clear();
        for (long i = 0; i < n; i++) if (d.x[i] < maxval + sig && d.x[i] > maxval - sig) { ix.push_back(d.x[i]); iy.push_back(d.y[i]); }
    };
    std::vector<double> ix, iy;
    double sig = cpeaks[0].x / 2.0;
    select(sig, ix, iy);
    if (ix.size() < 6) {
        sig = cpeaks[0].x;
        select(sig, ix, iy);
        if (ix.size() < 6 && cpeaks.size() > 1) { sig = cpeaks[1].x; select(sig, ix, iy); }
    }
    if (ix.size() < 3) return false;
    PeakFit fits[3];
    for (int k = 0; k < 3; k++) {
        fits[k] = isolated_peak_fit(ix, iy, k);
        if (!fits[k].ok || !std::isfinite(fits[k].error)) return false;
    }
    int best = 0;
    for (int k = 1; k < 3; k++) if (fits[k].error < fits[best].error) best = k;
    if (getenv("MS_DEBUG")) {
        fprintf(stderr, "apw: xdiff %.17g cpeak0 %.17g sig %.17g npts %zu\n", xdiff, cpeaks[0].x, sig, ix.size());
        for (int k = 0; k < 3; k++) fprintf(stderr, "apw fit %d: %.17g %.17g %.17g %.17g err %.17g\n", k, fits[k].p[0], fits[k].p[1], fits[k].p[2], fits[k].p[3], fits[k].error);
    }
    fwhm = std::nearbyint(fits[best].p[0] * 1e5) / 1e5;
    psfun_out = best;
    return true;
}

// ============================================================ what the engine can and cannot do
// Everything here follows the source of the bundled engine (UniDec 8.2.1, unidec/src: udmain.c,
// udcore.c, udstruct.c, udio.c) for the configuration this driver writes (conf.dat of write_conf).

// The engine's neighbourhood of a mass x charge point (SetUpBlur, for zzsig >= 0 and msig >= 0 as the
// window allows): charges j - zh .. j + zh with zh = int(zzsig), masses (unused offsets, molig is 0).
struct EngineBlur { int zlength = 3, mlength = 1; int closeness() const { return zlength * mlength; } };
EngineBlur engine_blur(double zzsig, double msig) {
    EngineBlur b;
    float zs = (float)zzsig, ms = (float)msig;
    if (zs >= 0 && ms >= 0) {
        b.zlength = 1 + 2 * (int)zs;
        b.mlength = 1 + 2 * (int)ms;
    } else {
        b.zlength = zs != 0 ? 1 + 2 * (int)(3 * std::fabs(zs) + 0.5f) : 1;
        b.mlength = ms != 0 ? 1 + 2 * (int)(3 * std::fabs(ms) + 0.5f) : 1;
    }
    return b;
}

// Peak memory (bytes) of the engine for L input points, Z charges, M masses of the output axis, C
// neighbours (closeness) and W points in the widest peak shape window (0 = linearised data, 1-d
// peak shape). From the allocations of the engine: input and mass table (mtab float, barr char per
// cell), the blur (closeind int and closearray float per neighbour and cell), blur, newblur and
// oldblur (float per cell), a second barr (char per cell), the peak shape (float per point and
// window point), small per point arrays, and either the temporary copy of point smoothing or
// suppression (float per cell, each iteration) or the outputs (mass grid, float per mass and charge).
double engine_memory(double L, double Z, double M, double C, double W, bool transient) {
    const double cells = L * Z;
    const double base = cells * (4.0 + 1.0 + 1.0 + 12.0 + 8.0 * C) + 4.0 * L * std::max(W, 1.0) + 40.0 * L;
    const double outputs = 4.0 * M * Z + 16.0 * M;
    return base + std::max(outputs, transient ? 4.0 * cells : 0.0) + 80e6;   // + the program itself (measured under Wine: 70 MB)
}

// The engine's 32 bit index arithmetic (int i * numz + j and so on) limits the arrays it indexes.
// Returns an empty string when the sizes fit, else what does not fit.
std::string engine_index_problem(double L, double Z, double M, double C, double W, bool speedy) {
    const double lim = 2147483647.0;
    if (C * L * Z > lim) return "the engine's neighbour table (" + std::to_string((long long)(C * L * Z)) + " entries)";
    if (L * Z > lim) return "the engine's m/z x charge grid";
    if (M * Z > lim) return "the engine's mass x charge grid";
    if (!speedy && L * std::max(W, 1.0) > lim) return "the engine's peak shape table";
    return "";
}

// The engine's nearest point search (udtools.c nearfast)
long engine_nearfast(const std::vector<float>& d, float point) {
    long start = 0, length = (long)d.size() - 1;
    long diff = length - start;
    while (diff > 1) {
        long mid = start + (length - start) / 2;
        if (point < d[mid]) length = mid;
        else if (point == d[mid]) return mid;
        else if (point > d[mid]) start = mid;
        diff = length - start;
    }
    return std::fabs(point - d[start]) >= std::fabs(point - d[length]) ? length : start;
}

// Widest peak shape window of the engine for data that are not linearised (SetStartsEnds): from the
// point nearest to x - 6 sigma to the one nearest to x + 6 sigma (sigma = FWHM / 2.35482 for the
// Gaussian, the FWHM for the other shapes), reflected at the low end as the engine does.
long engine_window_points(const std::vector<double>& x, double mzsig, int psfun) {
    float sig = (float)mzsig;
    if (psfun == 0) sig = sig / 2.35482f;
    const size_t n = x.size();
    if (n < 2) return 1;
    std::vector<float> d(n);
    for (size_t i = 0; i < n; i++) d[i] = (float)x[i];
    const float win = 6.0f * std::fabs(sig);
    long best = 1;
    for (size_t i = 0; i < n; i++) {
        const float lo = d[i] - win, hi = d[i] + win;
        const long start = lo < d[0] ? -engine_nearfast(d, 2.0f * d[0] - lo) : engine_nearfast(d, lo);
        const long end = hi > d[n - 1] ? (long)n - 1 + engine_nearfast(d, 2.0f * d[0] - hi) : engine_nearfast(d, hi);
        best = std::max(best, end - start);
    }
    return best;
}

// Which points of the mass x charge grid the engine keeps. SetLimits allows a point when its mass is
// inside the mass range; MakeSparseBlur then removes every point that has fewer than two allowed
// neighbours within 2 sigma of the predicted m/z (itself and the same mass at charge z +- 1, ...);
// KillB removes the points of zero intensity. When nothing is left the engine reports "Setup is
// bad. No points are allowed", writes its outputs from an empty fit array and Windows ends it with
// 0xC0000409 (fwrite of a null buffer in WriteDecon). The removal runs in parallel and changes the
// allowed set while it runs, so the outcome is decided here for sure only in two cases: a point of
// non zero intensity whose neighbour at z +- 1 has it as its own neighbour back (both keep two
// neighbours whatever the order: sure), or no point of non zero intensity with any allowed
// neighbour (nothing can be left).
struct EnginePoints {
    long in_range = 0;     // points of non zero intensity with a mass in the range at some charge
    bool sure = false;     // a point is certainly kept
    bool maybe = false;    // a point may be kept (depends on the order of the engine's threads)
};
EnginePoints engine_points(const std::vector<double>& x, const std::vector<double>& y, const UConfig& c) {
    EnginePoints r;
    const size_t L = x.size();
    if (L < 2) return r;
    // the values as the engine reads them: m/z "%.17g" and intensity "%.6f" into floats
    std::vector<float> mz(L), it(L);
    for (size_t i = 0; i < L; i++) {
        mz[i] = (float)x[i];
        char b[64];
        snprintf(b, sizeof b, "%.6f", y[i]);
        it[i] = strtof(b, nullptr);
    }
    const float adduct = strtof(py_float(c.adductmass).c_str(), nullptr);
    const float lb = strtof(py_float(c.masslb).c_str(), nullptr), ub = strtof(py_float(c.massub).c_str(), nullptr);
    float sig = strtof(py_float(c.mzsig).c_str(), nullptr);
    if (c.psfun == 0) sig = sig / 2.35482f;
    const float thr = sig * 2.0f;
    const EngineBlur bl = engine_blur(c.zzsig, c.msig);
    if ((float)c.zzsig < 0 || (float)c.msig < 0) { r.sure = r.maybe = true; return r; }   // not predicted
    const int numz = c.numz, zh = (bl.zlength - 1) / 2;
    auto zval = [&](int j) { return c.startz + j; };
    auto mass = [&](size_t i, int j) { return (float)zval(j) * (mz[i] - adduct); };
    auto allowed = [&](size_t i, int j) { float m = mass(i, j); return m < ub && m > lb; };
    // neighbour of (i, j) at charge index j2 for the same mass: its index, or -1
    auto neighbour = [&](size_t i, int j, int j2) -> long {
        if (j2 < 0 || j2 >= numz || zval(j2) == 0) return -1;
        const float point = (mass(i, j) + (float)zval(j2) * adduct) / (float)zval(j2);
        if (point < mz[0] - thr || point > mz[L - 1] + thr) return -1;
        const long ind = engine_nearfast(mz, point);
        if (!allowed((size_t)ind, j2) || !(std::fabs(point - mz[ind]) < thr)) return -1;
        return ind;
    };
    for (size_t i = 0; i < L; i++) {
        if (!(it[i] > 0)) continue;
        for (int j = 0; j < numz; j++) {
            if (!allowed(i, j)) continue;
            r.in_range++;
            if (bl.mlength >= 3) { r.sure = r.maybe = true; return r; }   // counts itself mlength times
            const bool self = neighbour(i, j, j) >= 0;
            int num = self ? 1 : 0;
            for (int dz = -zh; dz <= zh; dz++) {
                if (dz == 0) continue;
                long k = neighbour(i, j, j + dz);
                if (k < 0) continue;
                num++;
                if (self && neighbour((size_t)k, j + dz, j) == (long)i && neighbour((size_t)k, j + dz, j + dz) >= 0) {
                    r.sure = r.maybe = true;
                    return r;
                }
            }
            if (num >= 2) r.maybe = true;
        }
    }
    return r;
}

}  // namespace

// ============================================================ result object
struct ms_unidec_result {
    std::vector<double> mass_x, mass_y, fit_x, fit_y, z, zdist;
    std::vector<double> pk_mass, pk_height, pk_area, pk_score, pk_apex;
    std::vector<int> pk_niso;
    double r2 = 0, uniscore = 0;
    std::string notes, folder;
    double q_score = 0; int q_level = 2; std::string q_text;   // reliability (ms_unidec_quality)
};

namespace {

struct UPeak { double mass, height, area, score, apex; int n_iso; };
struct EScore { double mass, inten, dscore; };   // a peak line of the engine

// ============================================================ reliability of the result (ms_deconv.unidec_quality)
// From UniDec's own scores (Kostelic and Marty, Methods Mol. Biol. 2022, 2500, 159-180; UD_score.c): the
// level is decided on the intensity weighted DScore of the listed masses and on R squared separately
// (artifacts can fit well with low DScores; real masses with unsuited settings keep good DScores with a
// low R squared); the score is the UniScore. Thresholds from the test spectra (see ms_deconv).
const double QUALITY_D_POOR = 0.2, QUALITY_D_GOOD = 0.35, QUALITY_R2_POOR = 0.2, QUALITY_R2_GOOD = 0.6, QUALITY_D_LOW = 0.2;

struct Quality { double score; int level; std::string text; };   // level 0 good, 1 fair, 2 poor

std::string fmt(const char* f, ...) {
    char b[1200];
    va_list ap;
    va_start(ap, f);
    vsnprintf(b, sizeof b, f, ap);
    va_end(ap);
    return b;
}

Quality unidec_quality(const std::vector<UPeak>& peaks, double uniscore, double r2, double mzsig, double step, int z_main,
                       double mass_lo, double mass_hi, int z_lo, int z_hi, bool resolved, double width_set, double width_data) {
    Quality q;
    q.score = std::min(std::max(uniscore, 0.0), 1.0);
    if (peaks.empty()) {
        q.level = 2;
        q.text = "Unreliable result: UniDec found no mass above the peak threshold. Check the m/z, mass and charge ranges and the minimum intensity.";
        return q;
    }
    double sh = 0, shd = 0;
    size_t itop = 0;
    for (size_t k = 0; k < peaks.size(); k++) {
        double h2 = peaks[k].height * peaks[k].height;
        sh += h2; shd += h2 * peaks[k].score;
        if (peaks[k].height > peaks[itop].height) itop = k;
    }
    const double wd = shd / std::max(sh, 1e-300);
    const int n = (int)peaks.size();
    int n_low = 0;
    for (size_t k = 0; k < peaks.size(); k++) if (peaks[k].score < QUALITY_D_LOW && k != itop) n_low++;
    const double pct = 100.0 * std::max(r2, 0.0);
    q.level = (wd < QUALITY_D_POOR || r2 < QUALITY_R2_POOR) ? 2 : (wd >= QUALITY_D_GOOD && r2 >= QUALITY_R2_GOOD) ? 0 : 1;
    const bool iso_resolved = z_main > 0 && mzsig > 0 && mzsig < 0.8 * 1.00235 / z_main;
    if (q.level == 2 && wd < QUALITY_D_POOR) {
        if (resolved && r2 >= QUALITY_R2_POOR && n <= 10)
            q.text = fmt("Unreliable result: UniDec gives the masses low scores (mean DScore %.2f) because their isotope peaks overlap in the mass spectrum. Use a mass step of 0.5 to 1 Da for average masses, or a smaller peak width (now %.3g m/z) if the isotopes are resolved in the data.", wd, mzsig);
        else if (n >= 10)
            q.text = fmt("Unreliable result: UniDec gives its %d masses low scores (mean DScore %.2f), the pattern of a few peaks read at many charge states (harmonics), not of charge state series. The m/z range probably holds no protein of %g to %g Da at charges %d to %d: check where its charge states are and set the m/z, mass and charge ranges to match.", n, wd, mass_lo, mass_hi, z_lo, z_hi);
        else
            q.text = fmt("Unreliable result: UniDec gives the masses low scores (mean DScore %.2f): they are not supported by consistent charge state series. Check the m/z, mass and charge ranges, the minimum intensity and the peak width.", wd);
    } else if (q.level == 2) {
        std::string hint;
        if (width_set > 0 && width_data > 0 && width_set > 3.0 * width_data)
            hint = fmt("The peak width of %.3g m/z is much wider than the peaks of the data (about %.3g m/z): set it to 0 (measured).", width_set, width_data);
        else if (width_set > 0 && width_data > 0 && width_set < width_data / 3.0)
            hint = fmt("The peak width of %.3g m/z is much narrower than the peaks of the data (about %.3g m/z): set it to 0 (measured).", width_set, width_data);
        else
            hint = "Most of the signal is not explained by masses in the range: check the mass, charge and m/z ranges.";
        q.text = fmt("Unreliable result: UniDec explains only %.0f %% of the spectrum (R squared %.2f). ", pct, r2) + hint;
    } else if (q.level == 1 && r2 < QUALITY_R2_GOOD) {
        std::string hint = iso_resolved && step >= 0.1
            ? fmt("The isotope peaks are resolved in the data but a mass step of %g Da cannot follow them: for isotope resolved masses use 0.01 to 0.05 Da over a narrow mass range.", step)
            : std::string("Part of the signal is not explained: widen the mass or charge range, or check the baseline and other compounds in the m/z range.");
        q.text = fmt("The masses are supported by charge state series (mean DScore %.2f) but UniDec explains only %.0f %% of the spectrum (R squared %.2f). ", wd, pct, r2) + hint;
    } else if (q.level == 1) {
        q.text = fmt("UniDec's scores are moderate (mean DScore %.2f, R squared %.2f). Masses with a score below %.1f may be artifacts; narrower mass and charge ranges usually help.", wd, r2, QUALITY_D_LOW);
    } else {
        q.text = fmt("UniDec's scores are good (UniScore %.2f, mean DScore %.2f, R squared %.2f).", q.score, wd, r2);
        if (n_low)
            q.text += fmt(" %d minor mass%s ha%s a low score (below %.1f) and may be an artifact%s.", n_low, n_low > 1 ? "es" : "",
                          n_low > 1 ? "ve" : "s", QUALITY_D_LOW, n_low > 1 ? "s" : "");
    }
    return q;
}

// ============================================================ the engine in this process
// The engine runs on a thread of its own while this thread reports its iterations to the progress
// callback (10 to 90 of 100) and passes a cancellation (the callback returns 0) on to it within about
// 50 ms. (Its own thread also keeps the OpenMP thread count of the engine out of the caller's.)
void run_engine(const ms::udengine::Config& ec, const std::vector<float>& mz, const std::vector<float>& it, int threads,
                ms::udengine::Output& eo, ms_progress_fn cb, void* user) {
    std::atomic<int> cancel(0);
    std::atomic<long> iter(0);
    ms::udengine::Control ctl;
    ctl.cancel = &cancel;
    ctl.iteration = &iter;
    ctl.threads = threads;
    std::mutex mu;
    std::condition_variable cv;
    bool done = false;
    std::exception_ptr err;
    std::thread th([&] {
        try {
            ms::udengine::run(ec, mz.data(), it.data(), (int)mz.size(), eo, ctl);
        } catch (...) {
            err = std::current_exception();
        }
        std::lock_guard<std::mutex> lk(mu);
        done = true;
        cv.notify_one();
    });
    bool cancelled = false;
    try {
        const long numit = std::max(1, std::abs(ec.numit));
        std::unique_lock<std::mutex> lk(mu);
        while (!done) {
            cv.wait_for(lk, std::chrono::milliseconds(50));
            if (done || cancelled || !cb) continue;
            lk.unlock();
            const long i = 10 + 80 * std::min(iter.load(std::memory_order_relaxed), numit) / numit;
            if (!cb(user, i, 100)) {
                cancel.store(1);
                cancelled = true;
            }
            lk.lock();
        }
    } catch (...) {
        cancel.store(1);
        th.join();
        throw;
    }
    th.join();
    if (cancelled) throw ms::Cancelled();
    if (err) std::rethrow_exception(err);
}

// write_files 1: the outputs as unidec.exe writes them (udio.c WriteDecon, WriteGlobalOutputs; the
// peak lines of its log)
void write_engine_outputs(const std::string& prefix, const ms::udengine::Output& eo, const std::vector<double>& mx,
                          const std::vector<double>& my, const std::vector<EScore>& dscores, const ms::udengine::Config& ec) {
    auto bin = [&](const std::string& suffix, const std::vector<float>& v) {
        FILE* f = fopen((prefix + suffix).c_str(), "wb");
        if (!f) throw std::runtime_error("cannot write " + prefix + suffix);
        if (!v.empty()) fwrite(v.data(), sizeof(float), v.size(), f);
        fclose(f);
    };
    bin("_fitdat.bin", eo.fitdat);
    if (ec.baselineflag == 1) bin("_baseline.bin", eo.baseline);
    if (ec.rawflag == 0 || ec.rawflag == 1) bin("_massgrid.bin", eo.massgrid);
    std::string t;
    char b[160];
    for (size_t i = 0; i < mx.size(); i++) {
        snprintf(b, sizeof b, "%.6f %.6f\n", mx[i], my[i]);
        t += b;
    }
    write_text(prefix + "_mass.txt", t);
    snprintf(b, sizeof b, "error = %.6f\ntime = %.6f\niterations = %d\nuniscore = %.6f\n", ms::udengine::printed6(eo.error, false),
             eo.seconds, eo.iterations, ms::udengine::printed6(eo.uniscore, false));
    t = b;
    snprintf(b, sizeof b, "mzsig = %.6f\nzzsig = %.6f\nbeta = %.6f\npsig = %.6f\n", ms::udengine::printed6(ec.mzsig, false),
             ms::udengine::printed6(ec.zsig, false), ms::udengine::printed6(ec.beta, false), ms::udengine::printed6(ec.psig, false));
    write_text(prefix + "_error.txt", t + b);
    t = "UniDec engine in the msengine library\n";
    for (const EScore& d : dscores) {
        snprintf(b, sizeof b, "Peak: Mass: %.6f Int: %.6f DScore: %.6f \n", d.mass, d.inten, d.dscore);
        t += b;
    }
    snprintf(b, sizeof b, "Average Peaks Score (UniScore): %.6f\n", ms::udengine::printed6(eo.uniscore, false));
    write_text(prefix + "_log.txt", t + b);
}

void run_unidec(const double* mz, const double* it, long n0, std::string folder,
                std::string name, const Params& p, ms_unidec_result& out, ms_progress_fn cb, void* user) {
    using namespace ms::e2;
    // ---- check_grid (unidec)
    double lo = p.num("mass_lo", 0), hi = p.num("mass_hi", 0), step = p.num("mass_step", 1.0);
    double z0 = p.num("z_lo", 1), z1 = p.num("z_hi", 50);
    if (!std::isfinite(lo) || !std::isfinite(hi) || !std::isfinite(step) || !std::isfinite(z0) || !std::isfinite(z1))
        throw std::invalid_argument("Mass, charge and step values must be finite numbers");
    if (lo <= 0 || z0 < 1 || z1 < z0 || z0 != std::floor(z0) || z1 != std::floor(z1))
        throw std::invalid_argument("Use positive masses and whole-number charges in increasing order");
    if (z1 > 10000) throw std::invalid_argument("Use charges from 1 to 10000");
    if (step <= 0) throw std::invalid_argument("The mass step must be larger than 0");
    if (step < 0.001) {
        char b[256];
        snprintf(b, sizeof b, "A mass step of %g Da is finer than any mass spectrometer resolves (it would only make the calculation extremely slow); use 0.005 Da or more", step);
        throw std::invalid_argument(b);
    }
    if (hi <= lo) throw std::invalid_argument("The upper mass must be larger than the lower mass");
    double nm = (hi - lo) / step;
    int nz = (int)z1 - (int)z0 + 1;
    if (nz < 2)
        throw std::invalid_argument("UniDec needs at least two charge states (its engine stops with a single one, e.g. 1 to 1). Use a charge range such as 1 to 2, or Maximum entropy for one charge state.");
    {
        // the engine stores masses as 32 bit floats: below two float steps at the upper mass its mass
        // axis can no longer hold the step (equal or uneven neighbours)
        const double ulp = std::ldexp(1.0, (int)std::floor(std::log2(hi)) - 23);
        if (step < 2.0 * ulp) {
            const double nice[] = {0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0};
            double use = 1.0;
            for (double v : nice) if (v >= 2.0 * ulp) { use = v; break; }
            char b[400];
            snprintf(b, sizeof b, "A mass step of %g Da is finer than the UniDec engine can represent at %g Da (it stores masses as 32 bit numbers, %.2g Da apart there). Use a mass step of %g Da or more.",
                     step, hi, ulp, use);
            throw std::invalid_argument(b);
        }
    }
    {
        // charge and mass smoothing below 1 leave every point of the engine without the neighbour it
        // requires (MakeSparseBlur keeps a point only with two allowed neighbours), the engine then
        // stops with Windows error 0xC0000409 whatever the data
        const EngineBlur bl = engine_blur(p.num("zzsig", 1.0), p.num("msig", 0.0));
        if (bl.closeness() < 2)
            throw std::invalid_argument("UniDec's engine needs a charge smoothing of 1 or more (or a mass smoothing of 1 Da or more): with both below 1 it removes every point and stops (Windows code 0xC0000409). Use Charge smoothing 1, the UniDec default.");
    }
    {
        // the mass x charge grid alone (before the data are known): its float per cell, with the
        // engine's int index; the whole need is checked again when the engine input is ready
        const double cells = nm * nz;
        const double need = 4.0 * cells + 16.0 * nm;
        const double budget = ms::memory_budget();
        if (cells > 2147483647.0 || need > budget) {
            char b[600];
            snprintf(b, sizeof b, "Mass range %g to %g Da every %g Da with %d charges gives %.0f x %d grid points and would need more than %s of memory%s. Use a larger mass step (isotope peaks need 0.01 to 0.02 Da, envelopes 0.5 to 1 Da) or narrower mass and charge ranges.",
                     lo, hi, step, nz, nm, nz, ms::fmt_bytes(need).c_str(),
                     cells > 2147483647.0 ? " (more than the UniDec engine can index)" : ("; this computer can give it " + ms::fmt_bytes(budget)).c_str());
            throw std::invalid_argument(b);
        }
    }
    // ---- run_method: restriction, minimum intensity, baseline
    double mz_lo = p.num("mz_lo", 0), mz_hi = p.num("mz_hi", 0);
    Spec spec = restrict_spec(mz, it, n0, mz_lo, mz_hi);
    size_t n = spec.x.size();
    double base = n ? *std::max_element(spec.y.begin(), spec.y.end()) : 0.0;
    double min_int = p.num("min_intensity", 0);
    bool pct = p.inum("min_intensity_pct", 0) != 0 || p.str("min_intensity_unit") == "pct";
    double thr = 0;
    if (min_int > 0) thr = pct ? min_int / 100.0 * base : min_int;
    std::vector<char> keep;
    long npk = 0;
    if (thr > 0 && n) {
        auto pm = peak_mask(spec, thr);
        keep = pm.first; npk = pm.second;
    } else {
        keep.assign(n, 0);
        for (size_t i = 0; i < n; i++) keep[i] = spec.y[i] > 0;
        if (n > 2) for (size_t i = 1; i + 1 < n; i++) if (spec.y[i] > spec.y[i - 1] && spec.y[i] >= spec.y[i + 1] && spec.y[i] > 0) npk++;
    }
    long points_used = 0;
    for (char k : keep) points_used += k;
    std::string summary;
    if (n) {
        char b[512];
        snprintf(b, sizeof b, "m/z %.2f to %.2f", spec.x.front(), spec.x.back());
        summary = b;
        if (thr > 0) {
            double pc = base > 0 ? 100.0 * thr / base : 0.0;
            snprintf(b, sizeof b, ", peaks from %s up (%.3g %% of the base peak %s): %ld peaks kept whole, %ld of %ld points",
                     fmt_int(thr).c_str(), pc, fmt_int(base).c_str(), npk, points_used, (long)n);
            summary += b;
        } else {
            snprintf(b, sizeof b, ", no minimum intensity (base peak %s)", fmt_int(base).c_str());
            summary += b;
        }
    }
    if (p.inum("baseline", 0)) subtract_baseline(spec, p.num("baseline_width", 15.0) > 0 ? p.num("baseline_width", 15.0) : 15.0);
    double peak_width = p.num("peak_width", 0);
    int psfun = (int)p.inum("psfun", 0);
    bool have_psfun = p.has("psfun");
    std::string binning = p.str("binning", "auto");
    double mzbins_p = p.num("mzbins", 0);
    bool native_bins = false;
    if (thr > 0 && n > 10) {
        // UniDec with gaps of zeros: the peak width is measured before the weak peaks
        // are removed, and the data are put on bins of their own spacing
        if (!(peak_width > 0)) {
            double fw = 0; int ps = psfun;
            bool okw = auto_peak_width(spec, fw, ps);
            if (getenv("MS_DEBUG")) fprintf(stderr, "auto_peak_width (run_method) ok %d fwhm %.17g psfun %d\n", (int)okw, fw, ps);
            if (okw && std::isfinite(fw) && fw > 0 && fw < 0.2 * (spec.x.back() - spec.x.front())) {
                peak_width = fw;
                psfun = ps; have_psfun = true;
            }
        }
        if ((binning == "auto" || binning == "none") && n <= 60000) {
            binning = "linear";
            mzbins_p = median_diff(spec.x.data(), n);
            native_bins = true;
        }
    }
    if (thr > 0) {
        bool anypos = false;
        for (size_t i = 0; i < n; i++) { if (!keep[i]) spec.y[i] = 0.0; else if (spec.y[i] > 0) anypos = true; }
        if (npk < 1 || points_used < 5 || !anypos)
            throw std::invalid_argument("No data above the minimum intensity (" + fmt_int(thr) + "; the base peak in the m/z range is " +
                                        fmt_int(base) + "). Lower the minimum intensity.");
    }

    // ---- unidec_deconvolute
    if (n < 10) throw std::invalid_argument("Too few data points in the m/z range");
    if (name.empty()) name = "spectrum";
    // the files of a UniDec run (conf.dat, the engine input and its outputs) only when asked for:
    // MS Analysis reads none of them, the engine runs on the data in memory
    const bool write_files = p.inum("write_files", 0) != 0;
    std::string udir;
    ms::tick(cb, user, 2, 100);
    UConfig c;
    if (write_files) {
        unidec_place(folder, name);
        make_dirs(folder);
        folder = abs_path(folder);
        udir = join_path(folder, name + "_unidecfiles");
        make_dirs(udir);
        c.outfname = join_path(udir, name);
    } else {
        c.outfname = name;
    }
    c.infname = c.outfname + "_input.dat";
    std::string confname = c.outfname + "_conf.dat";
    c.minmz = spec.x.front(); c.maxmz = spec.x.back();
    c.startz = (int)z0; c.endz = (int)z1; c.numz = c.endz - c.startz + 1;
    c.masslb = lo; c.massub = hi; c.massbins = step;
    c.peakwindow = std::max(p.num("peak_window", 10.0), 2.0 * step);
    c.peakthresh = p.num("peak_threshold", p.num("peak_thresh", 0.1));
    {   // carrier_mass
        double a = p.num("adduct_mass", 0);
        double sign = p.num("sign", 1);
        if (a == 0 || !std::isfinite(a)) a = PROTON * (sign < 0 ? -1 : 1);
        c.adductmass = a;
    }
    if (p.has("numit")) c.numit = p.inum("numit", 100);
    if (p.has("zzsig")) c.zzsig = p.num("zzsig", 1);
    if (p.has("psig")) c.psig = p.num("psig", 1);
    if (p.has("beta")) c.beta = p.num("beta", 0);
    if (have_psfun) c.psfun = psfun;
    if (p.has("msig")) c.msig = p.num("msig", 0);
    if (p.has("isotopemode")) c.isotopemode = (int)p.inum("isotopemode", 0);
    if (p.has("poolflag")) c.poolflag = (int)p.inum("poolflag", 2);
    if (p.has("smooth")) c.smooth = p.num("smooth", 0);
    if (p.has("subbuff")) c.subbuff = p.num("subbuff", 0);
    if (p.has("subtype")) c.subtype = (int)p.inum("subtype", 2);
    if (p.has("peaknorm")) c.peaknorm = (int)p.inum("peaknorm", 1);
    if (binning == "none") { c.linflag = 2; c.mzbins = 0; }
    else if (binning == "linear") { c.linflag = 0; c.mzbins = mzbins_p > 0 ? mzbins_p : 0.01; }
    else if (binning == "resolution") {
        c.linflag = 1;
        c.mzbins = mzbins_p > 0 ? mzbins_p : spec.x.front() / (3.0 * std::max(estimate_resolution(spec.x.data(), spec.y.data(), n), 1000.0));
    } else if (n > 60000) {   // automatic
        if (step >= 0.1) {
            c.linflag = 0;
            c.mzbins = std::max(median_diff(spec.x.data(), n), (spec.x.back() - spec.x.front()) / 60000.0);
        } else {
            c.linflag = 1;
            c.mzbins = spec.x.front() / (3.0 * std::max(estimate_resolution(spec.x.data(), spec.y.data(), n), 1000.0));
        }
    }
    double width = peak_width > 0 ? peak_width : 0.0;
    if (width > 0) c.mzsig = width;
    c.intthresh = 0.0;
    int psfun0 = c.psfun;
    ms::tick(cb, user, 5, 100);
    // ---- process_data and the engine input
    Spec data2 = ud_dataprep(spec, c);
    {
        bool anypos = false;
        for (double v : data2.y) if (v > 0) { anypos = true; break; }
        if (data2.x.size() < 3 || !anypos)
            throw std::invalid_argument("No data left to deconvolute in the m/z range (check the m/z range and the minimum intensity)");
    }
    if (width <= 0) {
        double fw = 0; int ps = c.psfun;
        // measured on evenly spaced data: without bins the zeros between the peaks are removed
        // from data2, and the autocorrelation (lags in points) then mixes gaps with steps
        // (about 40 times too wide on clean spectra)
        Spec even;
        if (c.linflag == 2) {
            UConfig ce = c;
            ce.linflag = 0;
            ce.mzbins = 0;
            even = ud_dataprep(spec, ce);
        }
        bool ok = auto_peak_width(c.linflag == 2 ? even : data2, fw, ps);
        if (ok) { c.mzsig = fw; c.psfun = ps; }
        else { c.psfun = psfun0; c.mzsig = 0; }   // Python: the failed width leaves the default 0.85; measured below instead
        if (getenv("MS_DEBUG")) fprintf(stderr, "auto_peak_width ok %d fwhm %.17g psfun %d\n", (int)ok, c.mzsig, c.psfun);
        // the automatic width can fail on weak or noisy spectra
        if (!std::isfinite(c.mzsig) || c.mzsig <= 0 || c.mzsig > 0.2 * (spec.x.back() - spec.x.front()))
            c.mzsig = estimate_peak_width(spec.x.data(), spec.y.data(), n);
    }
    if (c.linflag == 0 && c.mzbins > 0 && c.mzsig < 2.5 * c.mzbins && !native_bins) c.mzsig = 2.5 * c.mzbins;  // a peak narrower than the bins crashes the engine
    double width_x = 0;   // m/z where the width shown in the notes applies (constant resolution bins)
    if (c.linflag == 1 && c.mzbins > 0 && data2.x.size() > 1) {
        // constant resolution bins: the engine's peak shape is a kernel in points set at the first bin
        // (MakePeakShape1D), so the width it applies grows with m/z; the width measured (or set) at the
        // tallest peak is converted to the first point, at least two bins (2.7 forced two of the widest
        // bins: on 295 to 3005 m/z ten times too wide at the low end, R squared -2.5)
        size_t it = std::max_element(data2.y.begin(), data2.y.end()) - data2.y.begin();
        width_x = data2.x[it];
        if (width_x > 0) c.mzsig *= data2.x.front() / width_x;
        c.mzsig = std::max(c.mzsig, 2.0 * c.mzbins);
    }
    // ---- what the engine will do with this input: its memory, its index limits, and whether any
    // point survives its neighbour rule (else it stops with 0xC0000409)
    std::string plan_note;
    auto check_engine = [&]() {
        const bool speedy = c.linflag != 2;
        const EngineBlur bl = engine_blur(c.zzsig, c.msig);
        const double L = (double)data2.x.size(), Z = c.numz, C = bl.closeness();
        // mass axis: the range plus the peak shape margin at the highest charge (SetupOutputs)
        float sig = (float)c.mzsig;
        if (c.psfun == 0) sig /= 2.35482f;
        const double M = (c.massub - c.masslb) / c.massbins + 2.0 * (6.0 * std::fabs((double)sig) * c.endz) / c.massbins + 2.0;
        const double W = speedy ? 1.0 : (double)engine_window_points(data2.x, c.mzsig, c.psfun);
        const std::string idx = engine_index_problem(L, Z, M, C, W, speedy);
        const double need = engine_memory(L, Z, M, C, W, c.psig >= 1 || c.beta > 0) + 40.0 * (double)n0;
        const double budget = ms::memory_budget();
        if (!idx.empty() || need > budget) {
            char b[900];
            snprintf(b, sizeof b, "UniDec would work on %.0f data points x %d charges (m/z %.2f to %.2f, %s) and %.0f masses: its engine would need about %s of memory%s. Use Data reduction: Automatic (bins), a narrower m/z or charge range, or a larger mass step.",
                     L, (int)Z, data2.x.front(), data2.x.back(), c.linflag == 2 ? "no data reduction" : "binned", M, ms::fmt_bytes(need).c_str(),
                     !idx.empty() ? (" and " + idx + " would exceed its 32 bit index").c_str() : ("; this computer can give it " + ms::fmt_bytes(budget)).c_str());
            throw std::invalid_argument(b);
        }
        return engine_points(data2.x, data2.y, c);
    };
    EnginePoints ep = check_engine();
    if (!ep.sure && c.linflag == 2 && ep.in_range > 0) {
        // data with gaps (zero runs removed, or isolated peaks): the masses lack the neighbouring
        // charge state the engine requires, and nothing may be left. On bins of the data's own
        // spacing every m/z has a point again (zero intensity ones count as neighbours), which is
        // what the engine expects; the peak width measured on the data is kept.
        const double bw = median_diff(spec.x.data(), n);
        c.linflag = 0;
        c.mzbins = bw;
        native_bins = true;
        data2 = ud_dataprep(spec, c);
        bool anypos = false;
        for (double v : data2.y) if (v > 0) { anypos = true; break; }
        if (data2.x.size() < 3 || !anypos)
            throw std::invalid_argument("No data left to deconvolute in the m/z range (check the m/z range and the minimum intensity)");
        char b[400];
        snprintf(b, sizeof b, "data put on linear bins of %.4g m/z (their own spacing): without data reduction their gaps left no mass with the neighbouring charge state the UniDec engine needs, and it would have stopped", bw);
        plan_note = b;
        ep = check_engine();
    }
    if (!ep.sure) {
        char b[1200];
        if (ep.in_range == 0)
            snprintf(b, sizeof b, "No data point in m/z %.2f to %.2f gives a mass between %g and %g Da at charges %d to %d. Check the m/z, mass and charge ranges.",
                     data2.x.front(), data2.x.back(), c.masslb, c.massub, c.startz, c.endz);
        else
            snprintf(b, sizeof b, "UniDec cannot deconvolute this spectrum with these settings: no peak in m/z %.2f to %.2f has the same mass at the next higher or lower charge (masses %g to %g Da, charges %d to %d), so the UniDec engine would remove every point and stop (Windows code 0xC0000409). The m/z range probably holds no charge state series, only isolated peaks (for example singly charged ions, or a minimum intensity above the protein signal). Check the ranges and the minimum intensity.",
                     data2.x.front(), data2.x.back(), c.masslb, c.massub, c.startz, c.endz);
        throw std::invalid_argument(b);
    }
    const std::string conf = conf_text(c);
    const size_t nd = data2.x.size();
    if (write_files) {
        write_text(confname, conf);
        FILE* f = fopen(c.infname.c_str(), "w");
        if (!f) throw std::runtime_error("cannot write " + c.infname);
        for (size_t i = 0; i < nd; i++) fprintf(f, "%.17g %.6f\n", data2.x[i], data2.y[i]);
        fclose(f);
    }
    ms::tick(cb, user, 10, 100);
    // ---- the engine, in this process: the configuration and the data as unidec.exe read them from
    // conf.dat and the input file
    const ms::udengine::Config ec = engine_config(conf);
    std::vector<float> emz(nd), eint(nd);
    for (size_t i = 0; i < nd; i++) {
        emz[i] = engine_mz(data2.x[i]);
        eint[i] = engine_intensity(data2.y[i]);
    }
    long threads = p.inum("threads", 0);   // unidec.exe used all logical processors unless told otherwise
    if (threads <= 0) threads = (long)std::max(1u, std::thread::hardware_concurrency());
    threads = std::min<long>(threads, ms::max_threads());
    ms::udengine::Output eo;
    run_engine(ec, emz, eint, (int)threads, eo, cb, user);
    std::vector<float>().swap(emz);
    std::vector<float>().swap(eint);
    ms::tick(cb, user, 90, 100);
    // ---- unidec_imports: what _mass.txt, _error.txt and the peak lines of unidec.exe gave ("%f" of
    // its C runtime, read back), the fit and the mass grid as the floats of _fitdat.bin and _massgrid.bin
    std::vector<double> mx(eo.mlen), my(eo.mlen);
    for (int i = 0; i < eo.mlen; i++) {
        mx[i] = ms::udengine::printed6(eo.massaxis[i], false);
        my[i] = ms::udengine::printed6(eo.massaxisval[i], false);
    }
    if (mx.empty()) throw std::runtime_error("The UniDec engine gave no mass spectrum.");
    std::vector<float>& massgrid = eo.massgrid;
    std::vector<float>& fitdat = eo.fitdat;
    {
        const double sse = ms::udengine::printed6(eo.error, false);
        double mean = np_sum(data2.y.data(), nd) / (double)nd;
        std::vector<double> sq(nd);
        for (size_t i = 0; i < nd; i++) sq[i] = (data2.y[i] - mean) * (data2.y[i] - mean);
        double ss = np_sum(sq.data(), nd);
        c.error = ss > 0 ? 1.0 - sse / ss : 0.0;
    }
    out.uniscore = ms::udengine::printed6(eo.uniscore, false);
    // DScores of the engine's peaks ("Peak: Mass: m Int: i DScore: d" of UD_score.c, printed for DScores above 0)
    std::vector<EScore> dscores;
    for (size_t k = 0; k < eo.peakx.size(); k++)
        if (eo.dscores[k] > 0)
            dscores.push_back({ms::udengine::printed6(eo.peakx[k], false), ms::udengine::printed6(eo.peaky[k], false),
                               ms::udengine::printed6(eo.dscores[k], false)});
    if (write_files) write_engine_outputs(c.outfname, eo, mx, my, dscores, ec);
    // ---- pick_peaks (tools.peakdetect, then the normalisation of setup_peaks)
    size_t nm_pts = mx.size();
    std::vector<XY> raw_peaks = ud_peakdetect(mx.data(), my.data(), (long)nm_pts, c.peakwindow / c.massbins, c.peakthresh, c.normthresh != 0);
    double massdatnormtop = *std::max_element(my.begin(), my.end());
    if (!raw_peaks.empty()) {
        double norm = 1;
        double pmax = 0, psum = 0;
        for (const XY& q : raw_peaks) { pmax = std::max(pmax, q.y); psum += q.y; }
        if (c.peaknorm == 1) norm = pmax / 100.0;
        else if (c.peaknorm == 2) norm = psum / 100.0;
        else norm = pmax / (massdatnormtop != 0 ? massdatnormtop : 1.0);
        if (norm == 0) norm = 1;
        for (XY& q : raw_peaks) q.y /= norm;
        for (double& v : my) v /= norm;
    }
    // ---- intensities in data units: the tallest mass gets the summed heights of its charge states
    std::vector<int> zs;
    for (int z = c.startz; z <= c.endz; z++) zs.push_back(z);
    double factor = 1.0;
    size_t itop = std::max_element(my.begin(), my.end()) - my.begin();
    double mymax = my[itop];
    bool have_zdist = false;
    std::vector<double> zdist(zs.size(), 0.0);
    if (nm_pts && mymax > 0) {
        double top = mx[itop];
        double a = c.adductmass;
        if (a == 1.0 || a == -1.0) a = PROTON * a;
        double w = std::max(c.mzsig, 1e-4);
        std::vector<double> heights(zs.size());
        for (size_t k = 0; k < zs.size(); k++) {
            double z = zs[k];
            double m = (top + z * a) / z;
            size_t lo_i = std::lower_bound(spec.x.begin(), spec.x.end(), m - w) - spec.x.begin();
            size_t hi_i = std::lower_bound(spec.x.begin(), spec.x.end(), m + w) - spec.x.begin();
            double h = 0;
            if (hi_i > lo_i) { h = spec.y[lo_i]; for (size_t i = lo_i + 1; i < hi_i; i++) h = std::max(h, spec.y[i]); }
            heights[k] = h;
        }
        std::vector<float> zrel;
        float zrel_max = 0;
        if (massgrid.size() == nm_pts * zs.size()) {
            zrel.assign(zs.size(), 0.0f);
            for (size_t i = 0; i < nm_pts; i++)
                for (size_t k = 0; k < zs.size(); k++) zrel[k] += massgrid[i * zs.size() + k];
            for (float v : zrel) zrel_max = std::max(zrel_max, v);
        }
        std::vector<char> use(zs.size());
        for (size_t k = 0; k < zs.size(); k++) {
            use[k] = heights[k] > 0;
            if (!zrel.empty() && zrel_max > 0) use[k] = use[k] && zrel[k] >= 0.05f * zrel_max;
        }
        std::vector<double> hu;
        double hmax_use = 0, hmax = 0;
        bool any_use = false;
        for (size_t k = 0; k < zs.size(); k++) {
            hmax = std::max(hmax, heights[k]);
            if (use[k]) { hu.push_back(heights[k]); hmax_use = std::max(hmax_use, heights[k]); any_use = true; }
        }
        double s_h = np_sum(hu.data(), hu.size());
        if (s_h > 0) factor = s_h / mymax;
        have_zdist = true;
        if (!zrel.empty() && zrel_max > 0 && hmax > 0) {
            double scale = any_use ? hmax_use : hmax;
            for (size_t k = 0; k < zs.size(); k++) zdist[k] = (double)(zrel[k] / zrel_max) * scale;
        } else {
            for (size_t k = 0; k < zs.size(); k++) zdist[k] = use[k] ? heights[k] : 0.0;
        }
    }
    for (double& v : my) v *= factor;
    // ---- fit on the data scale
    if (fitdat.size() == nd && nd) {
        double smax = *std::max_element(spec.y.begin(), spec.y.end());
        double dmax = *std::max_element(data2.y.begin(), data2.y.end());
        double ymax = smax / std::max(dmax, 1e-30);
        out.fit_x = data2.x;
        out.fit_y.resize(nd);
        for (size_t i = 0; i < nd; i++) out.fit_y[i] = (double)(float)(fitdat[i] * (float)ymax);
    }
    // ---- peaks: apex between the grid points, data units, one entry per apex
    std::vector<UPeak> peaks;
    for (const XY& q : raw_peaks) {
        UPeak u;
        u.mass = refine_peak_mass(mx.data(), my.data(), nm_pts, q.x);
        u.height = q.y * factor;
        u.area = 0.0;
        u.score = 0.0;
        u.apex = q.x;
        u.n_iso = 0;
        for (const auto& ds : dscores) if (std::fabs(ds.mass - q.x) <= 0.5 * c.massbins) { u.score = ds.dscore; break; }
        peaks.push_back(u);
    }
    std::stable_sort(peaks.begin(), peaks.end(), [](const UPeak& a, const UPeak& b) { return a.height > b.height; });
    std::vector<UPeak> kept;
    for (const UPeak& q : peaks) {
        bool ok = true;
        for (const UPeak& r : kept) if (!(std::fabs(q.mass - r.mass) > 0.3 * c.massbins)) { ok = false; break; }
        if (ok) kept.push_back(q);
    }
    std::sort(kept.begin(), kept.end(), [](const UPeak& a, const UPeak& b) { return a.mass < b.mass; });
    peaks = kept;
    bool grouped_used = false;
    if (c.massbins <= 0.25 && peaks.size() >= 5) {
        long close = 0;
        for (size_t i = 1; i < peaks.size(); i++) if (std::fabs(peaks[i].mass - peaks[i - 1].mass - 1.00235) < 0.06) close++;
        if (close >= 4) {
            std::vector<std::pair<double, double>> rng;
            std::vector<Peak> g = group_isotopes_ranges(mx.data(), my.data(), nm_pts, p.num("peak_threshold", p.num("peak_thresh", 0.1)), &rng);
            if (!g.empty()) {
                peaks.clear();
                for (size_t k = 0; k < g.size(); k++) {
                    // the species: the DScores of its isotope peaks, weighted by their intensities
                    double sw = 0, sd = 0;
                    for (const auto& ds : dscores)
                        if (ds.mass >= rng[k].first - 0.5 * c.massbins && ds.mass <= rng[k].second + 0.5 * c.massbins) { sw += ds.inten; sd += ds.inten * ds.dscore; }
                    peaks.push_back({g[k].mass, g[k].height, g[k].area, sw > 0 ? sd / sw : 0.0, g[k].apex, g[k].n_iso});
                }
                grouped_used = true;
            }
        }
    }
    // ---- notes
    {
        // the engine of UniDec 8.2.1 has no isotope mode (read and ignored: identical outputs for 0, 1 and 2)
        const char* iso_txt = c.isotopemode != 0 ? ", isotope mode ignored (the UniDec 8.2.1 engine has none)" : "";
        char b[1024];
        const double shown = width_x > 0 ? c.mzsig * width_x / data2.x.front() : c.mzsig;
        std::string wtxt = width_x > 0 ? fmt(" at m/z %.0f (constant resolution)", width_x) : std::string();
        snprintf(b, sizeof b, "UniDec: charge %d to %d, mass %g to %g Da every %g Da, peak width %.3g m/z%s%s", c.startz, c.endz, c.masslb, c.massub,
                 c.massbins, shown, wtxt.c_str(), iso_txt);
        out.notes = b;
        if (!plan_note.empty()) out.notes += "; " + plan_note;
        if (!summary.empty()) out.notes += "; " + summary;
        out.notes += fmt("; UniDec engine run in the library (%ld threads)", threads);
        int z_main = 0;
        if (have_zdist) {
            size_t kz = 0;
            for (size_t k = 1; k < zdist.size(); k++) if (zdist[k] > zdist[kz]) kz = k;
            if (!zdist.empty() && zdist[kz] > 0) z_main = zs[kz];
        }
        // a width set by hand (not the one measured before the minimum intensity was applied)
        const double width_set = std::max(p.num("peak_width", 0.0), 0.0);
        const double width_data = width_set > 0 ? estimate_peak_width(spec.x.data(), spec.y.data(), n) : 0.0;
        Quality q = unidec_quality(peaks, out.uniscore, c.error, shown, c.massbins, z_main, c.masslb, c.massub, c.startz, c.endz,
                                   grouped_used, width_set, width_data);
        out.q_score = q.score; out.q_level = q.level; out.q_text = q.text;
        if (q.level == 2) out.notes = q.text + "; " + out.notes;
    }
    out.r2 = c.error;
    out.folder = udir;
    out.mass_x = mx; out.mass_y = my;
    for (size_t k = 0; k < zs.size(); k++) { out.z.push_back(zs[k]); out.zdist.push_back(have_zdist ? zdist[k] : 0.0); }
    for (const UPeak& q : peaks) {
        out.pk_mass.push_back(q.mass); out.pk_height.push_back(q.height); out.pk_area.push_back(q.area);
        out.pk_score.push_back(q.score); out.pk_apex.push_back(q.apex); out.pk_niso.push_back(q.n_iso);
    }
}

}  // namespace

extern "C" {

MS_API int ms_unidec(const double* mz, const double* it, long n, const char* engine_exe, const char* work_dir,
                     const char* name, const char* params, ms_unidec_result** out, ms_progress_fn cb, void* user) {
    if (out) *out = nullptr;
    return ms::guarded([&] {
        (void)engine_exe;   // since 3.1 the engine runs in the library (kept for source compatibility; may be NULL)
        if (!mz || !it || n < 0 || !out) throw std::invalid_argument("ms_unidec: bad arguments");
        std::unique_ptr<ms_unidec_result> r(new ms_unidec_result);
        Params p = parse_params(params);
        run_unidec(mz, it, n, work_dir ? work_dir : ".", name ? name : "spectrum", p, *r, cb, user);
        *out = r.release();
    });
}

MS_API void ms_unidec_free(ms_unidec_result* r) { delete r; }

MS_API long ms_unidec_mass(const ms_unidec_result* r, const double** mass, const double** it) {
    if (!r) return 0;
    if (mass) *mass = r->mass_x.data();
    if (it) *it = r->mass_y.data();
    return (long)r->mass_x.size();
}
MS_API long ms_unidec_fit(const ms_unidec_result* r, const double** mz, const double** it) {
    if (!r) return 0;
    if (mz) *mz = r->fit_x.data();
    if (it) *it = r->fit_y.data();
    return (long)r->fit_x.size();
}
MS_API long ms_unidec_zdist(const ms_unidec_result* r, const double** z, const double** it) {
    if (!r) return 0;
    if (z) *z = r->z.data();
    if (it) *it = r->zdist.data();
    return (long)r->z.size();
}
MS_API long ms_unidec_peaks(const ms_unidec_result* r, const double** mass, const double** height, const double** area,
                            const double** score, const double** apex, const int** n_iso) {
    if (!r) return 0;
    if (mass) *mass = r->pk_mass.data();
    if (height) *height = r->pk_height.data();
    if (area) *area = r->pk_area.data();
    if (score) *score = r->pk_score.data();
    if (apex) *apex = r->pk_apex.data();
    if (n_iso) *n_iso = r->pk_niso.data();
    return (long)r->pk_mass.size();
}
MS_API double ms_unidec_r2(const ms_unidec_result* r) { return r ? r->r2 : 0.0; }
MS_API double ms_unidec_quality(const ms_unidec_result* r, int* level, const char** text) {
    if (level) *level = r ? r->q_level : 2;
    if (text) *text = r ? r->q_text.c_str() : "";
    return r ? r->q_score : 0.0;
}
MS_API double ms_unidec_uniscore(const ms_unidec_result* r) { return r ? r->uniscore : 0.0; }
MS_API const char* ms_unidec_notes(const ms_unidec_result* r) { return r ? r->notes.c_str() : ""; }
MS_API const char* ms_unidec_folder(const ms_unidec_result* r) { return r ? r->folder.c_str() : ""; }

MS_API const char* ms_unidec_engine(void) {
    return "UniDec 8.2.1 engine built into msengine (UniDec engine by Michael T. Marty; Marty et al., Anal. Chem. 2015, "
           "87, 4370-4376, DOI 10.1021/acs.analchem.5b00140)";
}

}

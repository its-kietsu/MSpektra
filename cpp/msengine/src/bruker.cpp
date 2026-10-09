// Bruker .d reader (agent E3): analysis.baf through Bruker's Baf2Sql library
// (baf2sql_c), analysis.tsf and analysis.tdf through the TDF SDK (timsdata).
// Mirrors hrms_data.BrukerD / BrukerTSF / BrukerTDF of the Python reference:
//   - MS1 scans = the rows with the smallest MS level (MsLevel / MsMsType),
//     ordered by time then id; rt = Time / 60; tic / bpc from the table (from
//     the line spectra when the table has none)
//   - one event per polarity with the acquisition m/z range
//   - TSF / TDF: profiles on the digitizer sample axis converted with
//     tsf_index_to_mz / tims_index_to_mz; sum_profiles() adds frames on the
//     sample axis and converts with the middle frame (as Python _sum does)
//   - TDF: ion mobility summed (tims_extract_profile_for_frame over all
//     scans), mass chromatograms through tims_extract_chromatograms
//   - the recalibration note of DataAnalysis (calibration.sqlite)
// The vendor libraries are loaded at run time (LoadLibrary / dlopen) from
// sdk_dir, or next to this library, each with its timsdata/ and baf2sql/
// sub folders.
#include "file.h"
#include "bruker_sys.h"
#include <algorithm>
#include <mutex>
#include <memory>
#include <unordered_map>
#include <cstdio>

namespace ms {

using sys::join;

// =========================================================== SDK loading
static std::string lib_name(const char* base) {
#ifdef _WIN32
    return std::string(base) + ".dll";
#else
    return std::string("lib") + base + ".so";
#endif
}

// candidate folders: sdk_dir (or the library's folder) and its sub folder
static std::string find_lib(const std::string& sdk_dir, const char* base, const char* sub) {
    std::vector<std::string> dirs;
    std::string root = sdk_dir.empty() ? sys::module_dir() : sdk_dir;
    if (!root.empty()) {
        dirs.push_back(root);
        dirs.push_back(join(root, sub));
    }
    if (sdk_dir.empty()) {    // the module folder's parent as well (python/ next to the library)
        std::string parent = sys::dirname(root);
        if (!parent.empty()) { dirs.push_back(join(parent, sub)); }
    }
    std::string name = lib_name(base);
    for (const std::string& d : dirs) {
        std::string p = join(d, name);
        if (sys::is_file(p)) return p;
    }
    std::string tried;
    for (const std::string& d : dirs) tried += (tried.empty() ? "" : ", ") + d;
    throw Error("Bruker's " + name + " was not found (looked in " + tried + ")");
}

extern "C" {
typedef void (*tims_profile_cb)(int64_t id, uint32_t n, const int32_t* values, void* user);
struct TimsChromJob {
    int64_t id;
    double time_begin, time_end, mz_min, mz_max, ook0_min, ook0_max;
};
typedef uint32_t (*tims_job_gen)(TimsChromJob* job, void* user);
typedef uint32_t (*tims_chrom_sink)(int64_t id, uint32_t n, const int64_t* frame_ids, const uint64_t* values, void* user);
}

struct TimsSdk {
    sys::DynLib lib;
    uint64_t (*tsf_open)(const char*, uint32_t) = nullptr;
    void (*tsf_close)(uint64_t) = nullptr;
    uint32_t (*tsf_get_last_error_string)(char*, uint32_t) = nullptr;
    uint32_t (*tsf_has_recalibrated_state)(uint64_t) = nullptr;
    int32_t (*tsf_read_line_spectrum_v2)(uint64_t, int64_t, double*, float*, int32_t) = nullptr;
    int32_t (*tsf_read_profile_spectrum_v2)(uint64_t, int64_t, uint32_t*, int32_t) = nullptr;
    uint32_t (*tsf_index_to_mz)(uint64_t, int64_t, const double*, double*, uint32_t) = nullptr;
    uint64_t (*tims_open)(const char*, uint32_t) = nullptr;
    void (*tims_close)(uint64_t) = nullptr;
    uint32_t (*tims_get_last_error_string)(char*, uint32_t) = nullptr;
    uint32_t (*tims_has_recalibrated_state)(uint64_t) = nullptr;
    uint32_t (*tims_index_to_mz)(uint64_t, int64_t, const double*, double*, uint32_t) = nullptr;
    uint32_t (*tims_extract_profile_for_frame)(uint64_t, int64_t, uint32_t, uint32_t, tims_profile_cb, void*) = nullptr;
    uint32_t (*tims_extract_chromatograms)(uint64_t, tims_job_gen, tims_chrom_sink, void*) = nullptr;

    explicit TimsSdk(const std::string& path) : lib(path) {
        lib.get(tsf_open, "tsf_open");
        lib.get(tsf_close, "tsf_close");
        lib.get(tsf_get_last_error_string, "tsf_get_last_error_string");
        lib.get(tsf_has_recalibrated_state, "tsf_has_recalibrated_state");
        lib.get(tsf_read_line_spectrum_v2, "tsf_read_line_spectrum_v2");
        lib.get(tsf_read_profile_spectrum_v2, "tsf_read_profile_spectrum_v2");
        lib.get(tsf_index_to_mz, "tsf_index_to_mz");
        lib.get(tims_open, "tims_open");
        lib.get(tims_close, "tims_close");
        lib.get(tims_get_last_error_string, "tims_get_last_error_string");
        lib.get(tims_has_recalibrated_state, "tims_has_recalibrated_state");
        lib.get(tims_index_to_mz, "tims_index_to_mz");
        lib.get(tims_extract_profile_for_frame, "tims_extract_profile_for_frame");
        lib.get(tims_extract_chromatograms, "tims_extract_chromatograms");
    }
    std::string tsf_error() const {
        uint32_t n = tsf_get_last_error_string(nullptr, 0);
        std::string buf(n > 0 ? n : 1, '\0');
        tsf_get_last_error_string(&buf[0], n);
        buf.resize(strlen(buf.c_str()));
        return buf.empty() ? "unknown TDF SDK error" : buf;
    }
    std::string tims_error() const {
        uint32_t n = tims_get_last_error_string(nullptr, 0);
        std::string buf(n > 0 ? n : 1, '\0');
        tims_get_last_error_string(&buf[0], n);
        buf.resize(strlen(buf.c_str()));
        return buf.empty() ? "unknown TDF SDK error" : buf;
    }
};

struct BafSdk {
    sys::DynLib lib;
    uint32_t (*get_sqlite_cache_filename_v2)(char*, uint32_t, const char*, int) = nullptr;
    uint64_t (*array_open_storage)(int, const char*) = nullptr;
    void (*array_close_storage)(uint64_t) = nullptr;
    int (*array_get_num_elements)(uint64_t, uint64_t, uint64_t*) = nullptr;
    int (*array_read_double)(uint64_t, uint64_t, double*) = nullptr;
    uint32_t (*get_last_error_string)(char*, uint32_t) = nullptr;

    explicit BafSdk(const std::string& path) : lib(path) {
        lib.get(get_sqlite_cache_filename_v2, "baf2sql_get_sqlite_cache_filename_v2");
        lib.get(array_open_storage, "baf2sql_array_open_storage");
        lib.get(array_close_storage, "baf2sql_array_close_storage");
        lib.get(array_get_num_elements, "baf2sql_array_get_num_elements");
        lib.get(array_read_double, "baf2sql_array_read_double");
        lib.get(get_last_error_string, "baf2sql_get_last_error_string");
    }
    std::string error() const {
        uint32_t n = get_last_error_string(nullptr, 0);
        std::string buf(n > 0 ? n : 1, '\0');
        get_last_error_string(&buf[0], n);
        buf.resize(strlen(buf.c_str()));
        return buf.empty() ? "unknown Baf2Sql error" : buf;
    }
};

// loaded once per process (the first sdk_dir wins, as the Python module does)
static std::mutex g_sdk_mutex;
static std::shared_ptr<TimsSdk> g_tims;
static std::shared_ptr<BafSdk> g_baf;

static std::shared_ptr<TimsSdk> tims_sdk(const std::string& sdk_dir) {
    std::lock_guard<std::mutex> lock(g_sdk_mutex);
    if (!g_tims) g_tims = std::make_shared<TimsSdk>(find_lib(sdk_dir, "timsdata", "timsdata"));
    return g_tims;
}
static std::shared_ptr<BafSdk> baf_sdk(const std::string& sdk_dir) {
    std::lock_guard<std::mutex> lock(g_sdk_mutex);
    if (!g_baf) g_baf = std::make_shared<BafSdk>(find_lib(sdk_dir, "baf2sql_c", "baf2sql"));
    return g_baf;
}

// ============================================================ .d folders
static bool is_tsf(const std::string& d) {
    return sys::is_dir(d) && sys::is_file(join(d, "analysis.tsf")) && sys::is_file(join(d, "analysis.tsf_bin"));
}
static bool is_tdf(const std::string& d) {
    return sys::is_dir(d) && sys::is_file(join(d, "analysis.tdf")) && sys::is_file(join(d, "analysis.tdf_bin"));
}
static bool is_d(const std::string& d) {
    return sys::is_dir(d) && (sys::is_file(join(d, "analysis.baf")) || is_tsf(d) || is_tdf(d));
}

// hrms_data.find_d_folder: the folder, a file inside it, or a folder with one .d
static std::string find_d_folder(const std::string& path0) {
    std::string path = sys::abspath(path0);
    if (sys::is_file(path)) {
        std::string base = sys::lower(sys::basename(path));
        if (base.compare(0, 9, "analysis.") == 0 || is_d(sys::dirname(path))) path = sys::dirname(path);
    }
    if (is_d(path)) return path;
    if (sys::is_dir(path)) {
        std::vector<std::string> subs;
        for (const std::string& n : sys::listdir(path))
            if (sys::ends_with(sys::lower(n), ".d") && is_d(join(path, n))) subs.push_back(join(path, n));
        if (subs.size() == 1) return subs[0];
    }
    return "";
}

bool is_bruker_d(const std::string& path, std::string* d_folder) {
    std::string d = find_d_folder(path);
    if (d.empty()) return false;
    if (d_folder) *d_folder = d;
    return true;
}

// hrms_data._recal_note: who recalibrated the file and when. Of several
// calibration states the latest (highest Id) is the one the Bruker libraries
// use ("ORDER BY Id DESC LIMIT 1" in timsdata), so its date is shown.
static std::string recal_note(const std::string& d) {
    try {
        sys::Db db(join(d, "calibration.sqlite"));
        std::map<std::string, std::string> info;
        db.query("SELECT KeyName, Value FROM CalibrationInfo WHERE KeyName IN ('CalibrationDateTime', 'CalibrationSoftware') "
                 "ORDER BY CalibrationState DESC",
                 [&](sqlite3_stmt* st) { info.emplace(sys::Db::text(st, 0), sys::Db::text(st, 1)); });
        std::string date = info.count("CalibrationDateTime") ? info["CalibrationDateTime"].substr(0, 10) : "";
        if (date.size() == 10 && date[4] == '-') date = date.substr(8, 2) + "." + date.substr(5, 2) + "." + date.substr(0, 4);
        std::string sw = info.count("CalibrationSoftware") ? info["CalibrationSoftware"] : "DataAnalysis";
        std::string out;
        for (const std::string& x : {sw, date}) if (!x.empty()) out += (out.empty() ? "" : ", ") + x;
        return out;
    } catch (...) {
        return "DataAnalysis";
    }
}

static double to_double(const std::string& s, bool* ok = nullptr) {
    char* end = nullptr;
    double v = strtod(s.c_str(), &end);
    bool good = !s.empty() && end && *end == '\0';
    if (ok) *ok = good;
    return good ? v : 0.0;
}

// =========================================================== base reader
class BrukerBase : public Reader {
public:
    const std::string& instrument() const override { return instrument_; }
    const std::string& recalibration() const override { return recal_; }
    long n_msms() const override { return n_msms_; }
    const std::vector<ScanInfo>& scans() const override { return scans_; }

    // centroid (line) spectrum of scan i, cached after the first read
    void centroids(long i, std::vector<double>& mz, std::vector<double>& it) override {
        check(i);
        std::lock_guard<std::mutex> lock(mutex_);
        if (!cached_[i]) {
            read_lines(i, line_mz_[i], line_it_[i]);
            cached_[i] = 1;
        }
        mz = line_mz_[i];
        it = line_it_[i];
    }

protected:
    std::string path_, instrument_, recal_;
    long n_msms_ = 0;
    std::vector<ScanInfo> scans_;
    std::map<std::string, std::string> props_;
    std::mutex mutex_;                         // the SDK handle and the line cache
    std::vector<char> cached_;
    std::vector<std::vector<double>> line_mz_, line_it_;

    void check(long i) const {
        if (i < 0 || i >= (long)scans_.size()) throw std::invalid_argument("scan index out of range");
    }
    void init_cache() {
        cached_.assign(scans_.size(), 0);
        line_mz_.resize(scans_.size());
        line_it_.resize(scans_.size());
    }
    // reads the line spectrum of scan i (no lock held by the caller's cache logic when called from tic fallback)
    virtual void read_lines(long i, std::vector<double>& mz, std::vector<double>& it) = 0;

    // events: one per polarity, sorted (negative first, as Python sorts 0/1 with 1 = '-')
    // (the file layer groups scans by ScanInfo.polarity; nothing to do here)

    // hrms_data: when the table has no TIC, the line spectra give tic and bpc
    void tic_from_lines(ms_progress_fn cb, void* user) {
        bool any = false;
        for (const ScanInfo& s : scans_) if (s.tic != 0) { any = true; break; }
        if (any) return;
        long n = (long)scans_.size();
        for (long i = 0; i < n; i++) {
            std::vector<double> mz, it;
            centroids(i, mz, it);
            double sum = 0, top = 0;
            for (double v : it) { sum += v; if (v > top) top = v; }
            scans_[i].tic = sum;
            scans_[i].bpc = it.empty() ? 0.0 : top;
            if (i % 200 == 0) tick(cb, user, i, n);
        }
    }
};

// =================================================================== BAF
class BrukerBaf : public BrukerBase {
public:
    BrukerBaf(const std::string& d, const std::string& sdk_dir, ms_progress_fn cb, void* user) : sdk_(baf_sdk(sdk_dir)) {
        path_ = d;
        work_ = d;
        std::string cache;
        try {
            cache = cache_file(d);
        } catch (const Error&) {
            // write protected folder: work on a temporary copy named after the
            // folder and its analysis.baf (path, size, time), as the Python version
            std::string name = sys::basename(d);
            size_t dot = name.rfind('.');
            if (dot != std::string::npos && dot > 0) name = name.substr(0, dot);
            sys::Stat st = sys::stat_path(join(d, "analysis.baf"));
            std::string tag = "copy";
            if (st.exists) {
                char buf[64];
                snprintf(buf, sizeof buf, "|%lld|%lld", (long long)st.size, (long long)st.mtime);
                tag = sys::md5_hex(sys::abspath(d) + buf).substr(0, 8);
            }
            std::string tmp = join(join(sys::temp_dir(), "MSAnalysis_d"), name + "_" + tag + ".d");
            if (!sys::is_dir(tmp) && (!sys::mkdir_p(sys::dirname(tmp)) || !sys::copy_tree_complete(d, tmp)))
                throw Error("could not copy " + d + " to " + tmp);
            work_ = tmp;
            cache = cache_file(tmp);
        }
        handle_ = sdk_->array_open_storage(0, work_.c_str());
        if (!handle_) throw Error(sdk_->error());
        // the destructor does not run when the constructor throws (cancelled, no spectra, SQLite or
        // SDK error): the storage is closed here then, or the files stay locked until the program ends
        struct CloseOnThrow {
            BrukerBaf* r;
            ~CloseOnThrow() { if (r && r->handle_) { r->sdk_->array_close_storage(r->handle_); r->handle_ = 0; } }
        } close_on_throw{this};
        struct Row {
            int64_t id; double rt, tic, bpc; int64_t pmz, pit, lmz, lit; double lo, hi; bool has_range; int64_t pol, level;
        };
        std::vector<Row> rows;
        {
            sys::Db db(cache);
            props_ = db.key_values("SELECT Key, Value FROM Properties");
            db.query("SELECT s.Id, s.Rt, s.SumIntensity, s.MaxIntensity, s.ProfileMzId, s.ProfileIntensityId, s.LineMzId, "
                     "s.LineIntensityId, s.MzAcqRangeLower, s.MzAcqRangeUpper, ak.Polarity, ak.MsLevel, ak.ScanMode "
                     "FROM Spectra s LEFT JOIN AcquisitionKeys ak ON s.AcquisitionKey = ak.Id ORDER BY s.Rt, s.Id",
                     [&](sqlite3_stmt* st) {
                         Row r;
                         r.id = sys::Db::integer(st, 0);
                         r.rt = sys::Db::num(st, 1);
                         r.tic = sys::Db::num(st, 2);
                         r.bpc = sys::Db::num(st, 3);
                         r.pmz = sys::Db::integer(st, 4);
                         r.pit = sys::Db::integer(st, 5);
                         r.lmz = sys::Db::integer(st, 6);
                         r.lit = sys::Db::integer(st, 7);
                         r.has_range = !sys::Db::is_null(st, 8) && !sys::Db::is_null(st, 9);
                         r.lo = sys::Db::num(st, 8);
                         r.hi = sys::Db::num(st, 9);
                         r.pol = sys::Db::integer(st, 10, -1);
                         r.level = sys::Db::integer(st, 11, 0);
                         rows.push_back(r);
                     });
        }
        if (rows.empty()) throw Error("No spectra in " + d);
        int64_t lmin = rows[0].level;
        for (const Row& r : rows) lmin = std::min(lmin, r.level);
        for (const Row& r : rows) {
            if (r.level != lmin) { n_msms_++; continue; }
            ScanInfo s;
            s.rt = r.rt / 60.0;
            s.tic = r.tic;
            s.bpc = r.bpc;
            s.polarity = r.pol == 1 ? -1 : 1;
            s.profile = r.pmz && r.pit;
            if (r.has_range) { s.mz_lo = r.lo; s.mz_hi = r.hi; }
            scans_.push_back(s);
            ids_.push_back({r.pmz, r.pit, r.lmz, r.lit});
        }
        auto p = props_.find("InstrumentName");
        if (p != props_.end()) instrument_ = p->second;
        init_cache();
        tic_from_lines(cb, user);
        recal_ = baf_recalibration();
        tick(cb, user, (long)scans_.size(), (long)scans_.size());
        close_on_throw.r = nullptr;   // opened: the destructor closes it
    }
    ~BrukerBaf() override {
        if (handle_) sdk_->array_close_storage(handle_);
    }
    const char* kind() const override { return "Bruker .d"; }

    void spectrum(long i, Spectrum& out) override {
        check(i);
        const Ids& r = ids_[i];
        if (r.pmz && r.pit) {
            std::lock_guard<std::mutex> lock(mutex_);
            read(r.pmz, out.mz);
            read(r.pit, out.it);
            out.profile = true;
            return;
        }
        centroids(i, out.mz, out.it);
        out.profile = false;
    }

protected:
    void read_lines(long i, std::vector<double>& mz, std::vector<double>& it) override {
        const Ids& r = ids_[i];
        if (r.lmz && r.lit) {
            read(r.lmz, mz);
            read(r.lit, it);
        } else if (r.pmz && r.pit) {
            std::vector<double> pm, pi;
            read(r.pmz, pm);
            read(r.pit, pi);
            centroid(pm.data(), pi.data(), (long)pm.size(), 0.002, mz, it);
        } else {
            mz.clear();
            it.clear();
        }
    }

private:
    struct Ids { int64_t pmz, pit, lmz, lit; };
    std::shared_ptr<BafSdk> sdk_;
    std::string work_;
    uint64_t handle_ = 0;
    std::vector<Ids> ids_;

    // A recalibration saved by DataAnalysis (calibration.sqlite or Calibrator.ami
    // in the .d folder) is applied by Baf2Sql when it reads the spectra (the
    // storage is opened with it); it is recognised by comparing the m/z of the
    // first scan with those of the calibration stored in analysis.baf
    // (hrms_data.BrukerD._baf_recal).
    std::string baf_recalibration() {
        long first = -1;
        for (size_t i = 0; i < ids_.size(); i++)
            if ((ids_[i].lmz && ids_[i].lit) || (ids_[i].pmz && ids_[i].pit)) { first = (long)i; break; }
        if (first < 0) return "";
        const int64_t ident = (ids_[first].lmz && ids_[first].lit) ? ids_[first].lmz : ids_[first].pmz;
        const uint64_t raw = sdk_->array_open_storage(1, work_.c_str());
        if (!raw) return "";
        std::vector<double> a, b;
        uint64_t n = 0;
        bool ok = sdk_->array_get_num_elements(raw, (uint64_t)ident, &n) && n > 0;
        if (ok) {
            a.resize(n);
            ok = sdk_->array_read_double(raw, (uint64_t)ident, a.data()) != 0;
        }
        sdk_->array_close_storage(raw);
        if (!ok) return "";
        try { read(ident, b); } catch (const Error&) { return ""; }
        if (a.size() != b.size()) return "";
        bool differs = false;
        for (size_t i = 0; i < a.size() && !differs; i++)
            if (std::fabs(a[i] - b[i]) > 1e-9 * std::max(std::fabs(b[i]), 1.0)) differs = true;
        if (!differs) return "";
        if (sys::is_file(join(path_, "calibration.sqlite"))) return recal_note(path_);
        return "DataAnalysis";
    }

    std::string cache_file(const std::string& d) {
        uint32_t n = sdk_->get_sqlite_cache_filename_v2(nullptr, 0, d.c_str(), 0);
        if (n == 0) throw Error(sdk_->error());
        std::string buf(n, '\0');
        sdk_->get_sqlite_cache_filename_v2(&buf[0], n, d.c_str(), 0);
        buf.resize(strlen(buf.c_str()));
        return buf;
    }
    void read(int64_t ident, std::vector<double>& out) {
        uint64_t n = 0;
        if (!sdk_->array_get_num_elements(handle_, (uint64_t)ident, &n)) throw Error(sdk_->error());
        out.resize(n);
        if (n && !sdk_->array_read_double(handle_, (uint64_t)ident, out.data())) throw Error(sdk_->error());
    }
};

// =================================================================== TSF
class BrukerTsf : public BrukerBase {
public:
    BrukerTsf(const std::string& d, const std::string& sdk_dir, ms_progress_fn cb, void* user, bool tdf = false)
        : sdk_(tims_sdk(sdk_dir)), tdf_(tdf) {
        path_ = d;
        // the recalibration of DataAnalysis is used when there is one (as for analysis.baf)
        handle_ = tdf_ ? sdk_->tims_open(d.c_str(), 1) : sdk_->tsf_open(d.c_str(), 1);
        if (!handle_) throw Error("The file could not be opened: " + sdk_error());
        // the destructor does not run when the constructor throws (cancelled, no spectra, SQLite or
        // SDK error): the handle is closed here then, or the files stay locked until the program ends
        struct CloseOnThrow {
            BrukerTsf* r;
            ~CloseOnThrow() { if (r) r->close_handle(); }
        } close_on_throw{this};
        uint32_t recal = tdf_ ? sdk_->tims_has_recalibrated_state(handle_) : sdk_->tsf_has_recalibrated_state(handle_);
        recal_ = recal ? recal_note(d) : "";
        struct Row { int64_t id; double t; std::string pol; int64_t msms; double tic, bpc; int64_t nscans; };
        std::vector<Row> rows;
        {
            sys::Db db(join(d, tdf_ ? "analysis.tdf" : "analysis.tsf"));
            props_ = db.key_values("SELECT Key, Value FROM GlobalMetadata");
            db.query(std::string("SELECT Id, Time, Polarity, MsMsType, SummedIntensities, MaxIntensity") +
                         (tdf_ ? ", NumScans" : "") + " FROM Frames ORDER BY Time, Id",
                     [&](sqlite3_stmt* st) {
                         Row r;
                         r.id = sys::Db::integer(st, 0);
                         r.t = sys::Db::num(st, 1);
                         r.pol = sys::Db::text(st, 2);
                         r.msms = sys::Db::integer(st, 3, 0);
                         r.tic = sys::Db::num(st, 4);
                         r.bpc = sys::Db::num(st, 5);
                         r.nscans = tdf_ ? sys::Db::integer(st, 6, 0) : 0;
                         rows.push_back(r);
                     });
        }
        if (rows.empty()) throw Error("No spectra in " + d);
        int64_t lmin = rows[0].msms;
        for (const Row& r : rows) lmin = std::min(lmin, r.msms);
        bool ok_lo, ok_hi;
        double lo = to_double(prop("MzAcqRangeLower"), &ok_lo), hi = to_double(prop("MzAcqRangeUpper"), &ok_hi);
        if (!ok_lo || !ok_hi) lo = hi = 0;
        for (const Row& r : rows) {
            if (r.msms != lmin) { n_msms_++; continue; }
            ScanInfo s;
            s.rt = r.t / 60.0;
            s.tic = r.tic;
            s.bpc = r.bpc;
            s.polarity = r.pol == "-" ? -1 : 1;
            s.mz_lo = lo;
            s.mz_hi = hi;
            scans_.push_back(s);
            ids_.push_back(r.id);
            nscans_.push_back(r.nscans);
            time_s_.push_back(r.t);
        }
        for (ScanInfo& s : scans_) s.profile = true;
        if (tdf_) {
            has_profile_ = true;
            has_line_ = false;
        } else {
            has_profile_ = prop("HasProfileSpectra", "1") == "1";
            has_line_ = prop("HasLineSpectra", "0") == "1";
            for (ScanInfo& s : scans_) s.profile = has_profile_;
        }
        const double nsamp = to_double(prop("DigitizerNumSamples", "0"));
        n_samples_ = (nsamp > 0 && nsamp < 1e9) ? (long)nsamp : 0;   // ((long)NaN is undefined; the SDK says what it needs anyway)
        prof_buf_ = std::max(n_samples_, 1024L);
        instrument_ = prop("InstrumentName");
        init_cache();
        if (!tdf_) tic_from_lines(cb, user);
        tick(cb, user, (long)scans_.size(), (long)scans_.size());
        close_on_throw.r = nullptr;   // opened: the destructor closes it
    }
    ~BrukerTsf() override { close_handle(); }
    const char* kind() const override { return tdf_ ? "Bruker .d (timsTOF)" : "Bruker .d (TSF)"; }

    // the profile of frame i on its own calibration (hrms_data._profile_frame)
    void spectrum(long i, Spectrum& out) override {
        check(i);
        if (!has_profile_) {
            centroids(i, out.mz, out.it);
            out.profile = false;
            return;
        }
        std::lock_guard<std::mutex> lock(mutex_);
        profile_raw(i, raw_);
        out.it.assign(raw_.begin(), raw_.end());
        index_axis(out.mz, out.it.size());
        to_mz(i, out.mz, out.mz);
        out.profile = true;
    }

    // hrms_data.BrukerTSF._sum: frames added on the digitizer sample axis,
    // converted with the calibration of the middle frame
    bool sum_profiles(const std::vector<long>& idx, Spectrum& out, ms_progress_fn cb, void* user) override {
        if (!has_profile_) return false;
        out.mz.clear();
        out.it.clear();
        out.profile = true;
        if (idx.empty()) return true;
        std::lock_guard<std::mutex> lock(mutex_);
        std::vector<double> acc;
        long n = (long)idx.size();
        for (long k = 0; k < n; k++) {
            if (k % 10 == 0) tick(cb, user, k, n);
            check(idx[k]);
            profile_raw(idx[k], raw_);
            if (raw_.empty()) continue;
            if (raw_.size() > acc.size()) acc.resize(raw_.size(), 0.0);
            for (size_t j = 0; j < raw_.size(); j++) acc[j] += raw_[j];
        }
        if (acc.empty()) return true;
        long mid = idx[idx.size() / 2];
        out.it = std::move(acc);
        index_axis(out.mz, out.it.size());
        to_mz(mid, out.mz, out.mz);
        return true;
    }

    // timsTOF: the SDK's chromatogram extraction (hrms_data.BrukerTDF._xic_raw)
    bool xic(int polarity, double mz, double tol, std::vector<double>& y) override {
        if (!tdf_) return false;
        double lo = mz - tol, hi = mz + tol;
        struct Ctx {
            TimsChromJob job;
            bool given = false;
            std::vector<int64_t> frames;
            std::vector<double> values;
            bool failed = false;
        } ctx;
        ctx.job = {1, -1.0, time_s_.back() + 1e6, std::min(lo, hi), std::max(lo, hi), 0.0, 1e3};
        tims_job_gen gen = [](TimsChromJob* job, void* user) -> uint32_t {
            Ctx* c = (Ctx*)user;
            if (c->given) return 2;
            *job = c->job;
            c->given = true;
            return 1;
        };
        tims_chrom_sink sink = [](int64_t, uint32_t n, const int64_t* frames, const uint64_t* values, void* user) -> uint32_t {
            Ctx* c = (Ctx*)user;
            try {   // no exception may pass through the frames of the SDK
                if (n && frames && values) {
                    c->frames.assign(frames, frames + n);
                    c->values.resize(n);
                    for (uint32_t i = 0; i < n; i++) c->values[i] = (double)values[i];
                }
            } catch (...) { c->failed = true; return 0; }
            return 1;
        };
        {
            std::lock_guard<std::mutex> lock(mutex_);
            if (!sdk_->tims_extract_chromatograms(handle_, gen, sink, &ctx)) throw Error(ctx.failed ? "out of memory" : sdk_error());
            if (ctx.failed) throw std::bad_alloc();
        }
        y.assign(scans_.size(), 0.0);
        if (frame_index_.empty())
            for (size_t i = 0; i < ids_.size(); i++) frame_index_.emplace(ids_[i], (long)i);
        for (size_t k = 0; k < ctx.frames.size(); k++) {
            auto f = frame_index_.find(ctx.frames[k]);
            if (f == frame_index_.end()) continue;
            if (scans_[f->second].polarity != polarity) continue;
            y[f->second] += ctx.values[k];
        }
        return true;
    }

protected:
    void read_lines(long i, std::vector<double>& mz, std::vector<double>& it) override {
        mz.clear();
        it.clear();
        if (has_line_) line_raw(i, mz, it);
        if (mz.empty() && has_profile_) {
            profile_raw(i, raw_);
            std::vector<double> pm(raw_.size()), pi(raw_.begin(), raw_.end());
            index_axis(pm, pm.size());
            to_mz(i, pm, pm);
            centroid(pm.data(), pi.data(), (long)pm.size(), 0.002, mz, it);
        }
    }

private:
    std::shared_ptr<TimsSdk> sdk_;
    bool tdf_;
    uint64_t handle_ = 0;
    void close_handle() {
        if (!handle_) return;
        if (tdf_) sdk_->tims_close(handle_); else sdk_->tsf_close(handle_);
        handle_ = 0;
    }
    std::vector<int64_t> ids_;
    std::vector<int64_t> nscans_;
    std::vector<double> time_s_;
    std::unordered_map<int64_t, long> frame_index_;
    bool has_profile_ = true, has_line_ = false;
    long n_samples_ = 0, prof_buf_ = 1024, line_buf_ = 65536;
    std::vector<uint32_t> prof32_;
    std::vector<double> raw_;             // profile intensities of the last frame read (under mutex_)
    std::vector<double> line_idx_;
    std::vector<float> line_val_;

    std::string prop(const char* key, const char* dflt = "") const {
        auto p = props_.find(key);
        return p == props_.end() ? std::string(dflt) : p->second;
    }
    std::string sdk_error() const { return tdf_ ? sdk_->tims_error() : sdk_->tsf_error(); }

    static void index_axis(std::vector<double>& idx, size_t n) {
        idx.resize(n);
        for (size_t j = 0; j < n; j++) idx[j] = (double)j;
    }
    // digitizer indices -> m/z with the calibration of frame i (in place allowed)
    void to_mz(long i, const std::vector<double>& index, std::vector<double>& out) {
        if (&index != &out) out.resize(index.size());
        if (index.empty()) return;
        uint32_t ok = tdf_ ? sdk_->tims_index_to_mz(handle_, ids_[i], index.data(), out.data(), (uint32_t)index.size())
                           : sdk_->tsf_index_to_mz(handle_, ids_[i], index.data(), out.data(), (uint32_t)index.size());
        if (!ok) throw Error(sdk_error());
    }
    // the profile intensities of frame i on the digitizer sample axis
    void profile_raw(long i, std::vector<double>& out) {
        if (tdf_) {
            struct Ctx { std::vector<double>* out; bool got = false; bool failed = false; } ctx{&out, false, false};
            out.clear();
            tims_profile_cb cb = [](int64_t, uint32_t n, const int32_t* v, void* user) {
                Ctx* c = (Ctx*)user;
                try {   // no exception may pass through the frames of the SDK
                    if (v) c->out->assign(v, v + n); else c->out->clear();
                    c->got = true;
                } catch (...) { c->failed = true; }
            };
            uint32_t end = (uint32_t)std::min<int64_t>(std::max<int64_t>(nscans_[i], 1), 0x7FFFFFFF);
            if (!sdk_->tims_extract_profile_for_frame(handle_, ids_[i], 0, end, cb, &ctx)) throw Error(sdk_error());
            if (ctx.failed) throw std::bad_alloc();
            return;
        }
        for (;;) {
            prof32_.resize(prof_buf_);
            int32_t n = sdk_->tsf_read_profile_spectrum_v2(handle_, ids_[i], prof32_.data(), (int32_t)prof_buf_);
            if (n < 0) throw Error(sdk_error());
            if (n > prof_buf_) { prof_buf_ = n; continue; }
            out.assign(prof32_.begin(), prof32_.begin() + n);
            return;
        }
    }
    void line_raw(long i, std::vector<double>& mz, std::vector<double>& it) {
        for (;;) {
            line_idx_.resize(line_buf_);
            line_val_.resize(line_buf_);
            int32_t n = sdk_->tsf_read_line_spectrum_v2(handle_, ids_[i], line_idx_.data(), line_val_.data(), (int32_t)line_buf_);
            if (n < 0) throw Error(sdk_error());
            if (n > line_buf_) { line_buf_ = n; continue; }
            line_idx_.resize(n);
            to_mz(i, line_idx_, mz);
            it.assign(line_val_.begin(), line_val_.begin() + n);
            return;
        }
    }
};

// ============================================================== factory
std::unique_ptr<Reader> open_bruker(const std::string& path, const std::string& sdk_dir, ms_progress_fn cb, void* user) {
    std::string d = find_d_folder(path);
    if (d.empty()) throw Error(path + " is not a Bruker .d folder (analysis.baf, analysis.tsf or analysis.tdf)");
    if (!sys::is_file(join(d, "analysis.baf"))) {
        if (is_tsf(d)) return std::unique_ptr<Reader>(new BrukerTsf(d, sdk_dir, cb, user, false));
        if (is_tdf(d)) return std::unique_ptr<Reader>(new BrukerTsf(d, sdk_dir, cb, user, true));
    }
    return std::unique_ptr<Reader>(new BrukerBaf(d, sdk_dir, cb, user));
}

}  // namespace ms

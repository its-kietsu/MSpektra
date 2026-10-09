// The ms_file API (events, chromatograms, averaging) on top of ms::Reader.
// Mirrors hrms_data.BrukerD.chromatogram / _sum / _sum_lines / average /
// scan_spectrum of the Python reference.
//
// Threads and handles. The GUI averages in a worker thread while its main thread
// asks for scans and chromatograms and may close the file (or open another one)
// meanwhile; the bridge may also close a file from a garbage collector thread.
//   - Every open file is in a registry. A call with NULL, with a closed file or
//     with any other value that is not an open file fails (MS_BAD_ARG, "" or 0)
//     without touching memory; ms_close of such a value does nothing.
//   - A call holds a reference to the file while it runs. ms_close takes the file
//     out of the registry at once; calls still running stop at their next progress
//     tick (MS_CANCELLED) and the last of them frees the file.
//   - The reader is used under the file's mutex, one reader call at a time. The
//     loops of this file take it per scan, so a scan or chromatogram asked for
//     during an average waits for one scan, not for the whole average (readers that
//     sum on their own, TSF, TDF and Shimadzu, hold it for their sum).
//   - Arrays and strings handed out belong to the calling thread: they stay valid
//     until the same thread makes the same call on the same file again (for
//     ms_event_scans: for the same event), or makes a file call after the file was
//     closed, or ends. A worker thread's average therefore stays readable even if
//     the main thread closes the file right after the call returned.
//   - A progress callback must not call functions on the file it reports for
//     (MS_BAD_ARG: the reader is busy in the same thread); it may call ms_close.
#include "file.h"
#include "fileio.h"
#include <algorithm>
#include <atomic>
#include <cctype>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <mutex>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#ifdef _OPENMP
#include <omp.h>
#endif

struct ms_file {
    std::unique_ptr<ms::Reader> reader;
    std::string summary;
    struct Event {
        int polarity = 1;
        std::vector<long> idx;                 // scan indices (reader order)
        std::vector<double> rt, tic, bpc;
        std::vector<unsigned char> flags;      // bit 0: saturated readings left out (since 2.8)
        double lo = 0, hi = 0;                 // 0 = unknown
    };
    std::vector<Event> events;                 // fixed after ms_open
    bool has_profile = false;
    // concurrency (see the top of this file)
    uint64_t serial = 0;                       // unique per opened file, never reused
    std::mutex reader_mu;                      // one reader call at a time
    std::atomic<bool> closing{false};          // ms_close was called
    int busy = 0;                              // calls running on this file (under the registry mutex)
};

namespace {

// ================================================================== registry
struct Registry {
    std::mutex mu;
    std::unordered_map<const ms_file*, uint64_t> live;    // open files
    std::unordered_set<uint64_t> serials;                 // their serial numbers
    uint64_t next = 1;
};
Registry& registry() {
    static Registry* r = new Registry;   // never destroyed: other threads may still call while the process ends
    return *r;
}

// the file a progress callback of this thread is reporting for (re-entry guard)
thread_local const ms_file* t_callback_file = nullptr;

// A reference to an open file for the duration of a call. ms_close while the
// call runs: the file is freed when the last reference goes.
class FileRef {
public:
    // uses_reader: the call needs the reader (not allowed from the file's own progress callback,
    // where the reader is busy in this very thread)
    explicit FileRef(const ms_file* f, bool uses_reader = true) {
        if (!f) { why_ = "null file"; return; }
        if (uses_reader && f == t_callback_file) { why_ = "called from a progress callback of the same file"; return; }
        Registry& r = registry();
        std::lock_guard<std::mutex> lk(r.mu);
        if (!r.live.count(f)) { why_ = "not an open file (closed already?)"; return; }
        f_ = const_cast<ms_file*>(f);
        f_->busy++;
    }
    ~FileRef() {
        if (!f_) return;
        bool last = false;
        {
            std::lock_guard<std::mutex> lk(registry().mu);
            last = --f_->busy == 0 && f_->closing.load();
        }
        if (last) delete f_;
    }
    FileRef(const FileRef&) = delete;
    FileRef& operator=(const FileRef&) = delete;
    explicit operator bool() const { return f_ != nullptr; }
    ms_file* operator->() const { return f_; }
    ms_file* get() const { return f_; }
    // throws std::invalid_argument("<fn>: <why>") when there is no open file
    void require(const char* fn) const {
        if (!f_) throw std::invalid_argument(std::string(fn) + ": " + why_);
    }
private:
    ms_file* f_ = nullptr;
    const char* why_ = "";
};

// ============================================================ thread outputs
struct Outputs {
    uint64_t serial = 0;
    std::vector<double> xic_y, avg_mz, avg_it, scan_mz, scan_it, pda_wl, pda_t, pda_a;
    std::vector<std::vector<double>> ev_rt, ev_tic, ev_bpc;
    std::string summary, recal, info;
};
thread_local std::vector<std::unique_ptr<Outputs>> t_outputs;

Outputs& outputs_of(const ms_file* f) {
    // arrays of files closed since are dropped first (their validity ended with the close)
    std::vector<std::unique_ptr<Outputs>> dead;
    {
        Registry& r = registry();
        std::lock_guard<std::mutex> lk(r.mu);
        for (size_t i = 0; i < t_outputs.size();) {
            if (t_outputs[i]->serial != f->serial && !r.serials.count(t_outputs[i]->serial)) {
                dead.push_back(std::move(t_outputs[i]));
                t_outputs.erase(t_outputs.begin() + (long)i);
            } else {
                i++;
            }
        }
    }
    dead.clear();
    for (auto& o : t_outputs)
        if (o->serial == f->serial) return *o;
    t_outputs.emplace_back(new Outputs);
    t_outputs.back()->serial = f->serial;
    return *t_outputs.back();
}

// ======================================================= progress callbacks
// What the readers and the loops below get as their callback: stops (returns 0)
// once the file is being closed, and marks the thread while the caller's
// callback runs (re-entry guard).
struct Progress {
    ms_progress_fn cb = nullptr;
    void* user = nullptr;
    const ms_file* f = nullptr;
};
int progress_wrap(void* u, long i, long n) {
    const Progress* p = static_cast<const Progress*>(u);
    if (p->f && p->f->closing.load(std::memory_order_relaxed)) return 0;
    if (!p->cb) return 1;
    struct Mark {
        const ms_file* prev;
        explicit Mark(const ms_file* f) : prev(t_callback_file) { t_callback_file = f; }
        ~Mark() { t_callback_file = prev; }
    } mark(p->f);
    const int r = p->cb(p->user, i, n);
    if (p->f && p->f->closing.load(std::memory_order_relaxed)) return 0;
    return r;
}

// ================================================================== helpers
bool ends_with_ci(const std::string& s, const char* suffix) {
    const size_t n = std::strlen(suffix);
    if (s.size() < n) return false;
    for (size_t i = 0; i < n; i++)
        if (std::tolower((unsigned char)s[s.size() - n + i]) != std::tolower((unsigned char)suffix[i])) return false;
    return true;
}

// a data file given as a path must be a readable regular file (not a folder named .mzML)
void require_regular_file(const std::string& p, const char* what) {
    switch (ms::path_type(p)) {
        case 1: return;
        case 0: throw ms::Error("cannot open " + p + " (no such file)");
        case 2: throw ms::Error(p + " is a folder, not " + what);
        default: throw ms::Error(p + " is not a regular file");
    }
}

const ms_file::Event& event_of(const ms_file* f, int event) {
    if (event < 0 || event >= (int)f->events.size()) throw std::invalid_argument("no such scan event");
    return f->events[event];
}

// scans of the event with rt in [t0, t1]
std::vector<long> scans_in(const ms_file::Event& ev, double t0, double t1) {
    const double lo = std::min(t0, t1), hi = std::max(t0, t1);
    std::vector<long> out;
    for (size_t k = 0; k < ev.idx.size(); k++)
        if (ev.rt[k] >= lo && ev.rt[k] <= hi) out.push_back(ev.idx[k]);
    return out;
}

// hrms_data._clean_axis: points with a finite m/z above 0, sorted by m/z (stable)
void clean_axis(std::vector<double>& m, std::vector<double>& it) {
    size_t w = 0;
    bool sorted = true;
    for (size_t i = 0; i < m.size(); i++) {
        if (!(m[i] > 0) || !std::isfinite(m[i])) continue;
        if (w && m[i] < m[w - 1]) sorted = false;
        m[w] = m[i];
        it[w] = it[i];
        w++;
    }
    m.resize(w);
    it.resize(w);
    if (sorted) return;
    std::vector<size_t> o(w);
    std::iota(o.begin(), o.end(), (size_t)0);
    std::stable_sort(o.begin(), o.end(), [&](size_t a, size_t b) { return m[a] < m[b]; });
    std::vector<double> m2(w), i2(w);
    for (size_t k = 0; k < w; k++) { m2[k] = m[o[k]]; i2[k] = it[o[k]]; }
    m.swap(m2);
    it.swap(i2);
}

// hrms_data._interp_gap at one x: linear interpolation of the profile (xp, fp) that reads a gap of a
// zero trimmed profile (two points more than 1.5 x the spacing next to them apart) as zero
// (j: the last index with xp[j] <= x, or -1)
inline double interp_gap_at(double x, long j, const double* xp, const double* fp, long n) {
    if (j < 0 || !(x <= xp[n - 1])) return 0.0;
    if (xp[j] == x) return fp[j];
    if (j >= n - 1) return 0.0;
    const double h = xp[j + 1] - xp[j];
    double ref = std::numeric_limits<double>::infinity();
    if (j >= 1) ref = xp[j] - xp[j - 1];
    if (j + 2 <= n - 1) ref = std::min(ref, xp[j + 2] - xp[j + 1]);
    if (h > 1.5 * ref) return 0.0;
    return fp[j] + (fp[j + 1] - fp[j]) * (x - xp[j]) / h;
}

double interp_gap(double x, const double* xp, const double* fp, long n) {
    if (n <= 0) return 0.0;
    return interp_gap_at(x, (long)(std::upper_bound(xp, xp + n, x) - xp) - 1, xp, fp, n);   // searchsorted right - 1
}

// hrms_data._resample_add: adds the profile (m, it) to the sum acc on the m/z axis grid (both clean),
// read on the grid by linear interpolation; its points farther than 0.75 x its point spacing from every
// grid point are added to the grid first (with the sum read there from the grid)
void resample_add(std::vector<double>& grid, std::vector<double>& acc, std::vector<double>& m, std::vector<double>& it) {
    clean_axis(m, it);
    const long n = (long)m.size(), g = (long)grid.size();
    if (n == 0) return;
    if (n == g && std::equal(grid.begin(), grid.end(), m.begin())) {
        for (long k = 0; k < g; k++) acc[k] += it[k];
        return;
    }
    const double inf = std::numeric_limits<double>::infinity();
    std::vector<long> pos;   // (grid index before which it goes) of the new points, in order
    std::vector<long> which;
    long j = 0;
    for (long i = 0; i < n; i++) {
        double sp = 0.0;
        if (n > 1) sp = std::min(i > 0 ? m[i] - m[i - 1] : inf, i < n - 1 ? m[i + 1] - m[i] : inf);
        while (j < g && grid[j] < m[i]) j++;   // searchsorted left (m sorted)
        const double dl = j > 0 ? m[i] - grid[j - 1] : inf, dr = j < g ? grid[j] - m[i] : inf;
        if (std::min(dl, dr) > 0.75 * sp) { pos.push_back(j); which.push_back(i); }
    }
    if (!pos.empty()) {
        const long nn = (long)pos.size();
        std::vector<double> ng, na;
        ng.reserve(g + nn);
        na.reserve(g + nn);
        long q = 0;
        for (long k = 0; k <= g; k++) {
            for (; q < nn && pos[q] == k; q++) {   // np.insert: before grid[k]
                ng.push_back(m[which[q]]);
                na.push_back(interp_gap(m[which[q]], grid.data(), acc.data(), g));
            }
            if (k < g) { ng.push_back(grid[k]); na.push_back(acc[k]); }
        }
        grid.swap(ng);
        acc.swap(na);
    }
    const long a = (long)(std::lower_bound(grid.begin(), grid.end(), m[0]) - grid.begin());
    const long b = (long)(std::upper_bound(grid.begin(), grid.end(), m[n - 1]) - grid.begin());
    long q = 0;   // the last index with m[q] <= grid[k] (grid[k] >= m[0] here)
    for (long k = a; k < b; k++) {
        while (q + 1 < n && m[q + 1] <= grid[k]) q++;
        acc[k] += interp_gap_at(grid[k], q, m.data(), it.data(), n);
    }
}

// hrms_data._merge_sticks: summed centroids (sorted) less than ppm apart from the one before are one
// ion (its centroid scatters from scan to scan): one stick at the intensity weighted m/z, intensities summed
const double MERGE_STICKS_PPM = 3.0;
void merge_sticks(std::vector<double>& mz, std::vector<double>& it, double ppm) {
    const size_t n = mz.size();
    if (n < 2) return;
    std::vector<double> om, oi;
    om.reserve(n);
    oi.reserve(n);
    size_t i = 0;
    while (i < n) {
        double si = 0.0, sw = 0.0, sm = 0.0;
        size_t j = i;
        do {
            const double w = std::max(it[j], 1e-30);
            si += it[j];
            sw += w;
            sm += mz[j] * w;
            j++;
        } while (j < n && !(mz[j] - mz[j - 1] > mz[j - 1] * (ppm * 1e-6)));
        if (j == i + 1) { om.push_back(mz[i]); oi.push_back(it[i]); }   // one stick: as it is
        else { om.push_back(sm / sw); oi.push_back(si); }
        i = j;
    }
    mz.swap(om);
    it.swap(oi);
}

// LogBins (1 ppm) of the centroids of the scans, the sticks of one ion merged when several scans are
// summed (= BrukerD._sum_lines)
void sum_lines(ms_file* f, const std::vector<long>& scans, ms_progress_fn cb, void* user,
               std::vector<double>& mz, std::vector<double>& it) {
    ms::LogBins bins(1e-6);
    std::vector<double> cm, ci;
    const long n = (long)scans.size();
    long used = 0;
    for (long k = 0; k < n; k++) {
        if (k % 50 == 0) ms::tick(cb, user, k, n);
        {
            std::lock_guard<std::mutex> lk(f->reader_mu);
            f->reader->centroids(scans[k], cm, ci);
        }
        if (!cm.empty()) { bins.add(cm.data(), ci.data(), (long)cm.size()); used++; }
    }
    bins.result(mz, it);
    if (used > 1 && mz.size() > 1) merge_sticks(mz, it, MERGE_STICKS_PPM);
}

// summed spectrum of the scans (= BrukerD._sum)
void sum_scans(ms_file* f, const std::vector<long>& scans, ms_progress_fn cb, void* user,
               std::vector<double>& mz, std::vector<double>& it) {
    mz.clear();
    it.clear();
    std::vector<double> grid, acc;
    bool have_grid = false, clean = false;
    ms::Spectrum sp;
    bool own = false, common = false;
    {
        // a reader may sum on its own (TSF/TDF: digitizer samples, Shimadzu: bins)
        std::lock_guard<std::mutex> lk(f->reader_mu);
        own = f->reader->sum_profiles(scans, sp, cb, user);
        common = f->reader->common_axis();
    }
    if (own) {
        mz.swap(sp.mz);
        it.swap(sp.it);
        return;
    }
    const long n = (long)scans.size();
    const int nt = ms::threads();
    for (long k = 0; k < n; k++) {
        if (k % 10 == 0) ms::tick(cb, user, k, n);
        {
            std::lock_guard<std::mutex> lk(f->reader_mu);
            f->reader->spectrum(scans[k], sp);
        }
        if (sp.mz.empty()) continue;
        if (sp.it.size() != sp.mz.size()) throw ms::Error("a spectrum of the file has m/z and intensity arrays of different length");
        if (!sp.profile) { sum_lines(f, scans, cb, user, mz, it); return; }
        if (!have_grid) {
            grid.swap(sp.mz);
            acc.swap(sp.it);
            have_grid = true;
        } else if (sp.mz.size() == grid.size() && (common || std::equal(grid.begin(), grid.end(), sp.mz.begin()))) {
            const long m = (long)grid.size();
            const double* src = sp.it.data();
            double* dst = acc.data();
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt) if (m > 100000)
#endif
            for (long j = 0; j < m; j++) dst[j] += src[j];
        } else {
            // axes differ (zero trimmed profiles): each spectrum read on the axis of the sum by linear
            // interpolation, the axis taking in its points where it has none
            if (!clean) { clean_axis(grid, acc); clean = true; }
            resample_add(grid, acc, sp.mz, sp.it);
        }
    }
    (void)nt;
    if (!have_grid) return;
    mz.swap(grid);
    it.swap(acc);
}

// np.interp(x, xp, fp, left=0, right=0) for increasing xp
double interp0(double x, const std::vector<double>& xp, const std::vector<double>& fp) {
    const size_t n = xp.size();
    if (n == 0 || x < xp[0] || x > xp[n - 1]) return 0.0;
    if (n == 1) return fp[0];
    size_t j = (size_t)(std::upper_bound(xp.begin(), xp.end(), x) - xp.begin());
    if (j == 0) return fp[0];
    j--;
    if (j >= n - 1) return fp[n - 1];
    const double slope = (fp[j + 1] - fp[j]) / (xp[j + 1] - xp[j]);
    return slope * (x - xp[j]) + fp[j];
}

// the message of a cancellation caused by ms_close
void closed_cancel(const FileRef& ref) {
    if (ref->closing.load()) throw ms::Cancelled("cancelled: the file was closed");
}

// ms_open, ms_open_arrays and ms_open_vendor: the reader given (arrays, vendor) or the one of the
// path's format, then the events, the summary and the registry (throws on failure)
ms_file* open_file(const std::string& p, std::unique_ptr<ms::Reader> given, const char* sdk_dir, ms_progress_fn cb,
                   void* user) {
        ms_file* out = nullptr;
        std::unique_ptr<ms_file> f(new ms_file());
        Progress pr;
        pr.cb = cb;
        pr.user = user;
        std::string d;
        if (given) {
            f->reader = std::move(given);
        } else if (ends_with_ci(p, ".mzml")) {
            require_regular_file(p, "an mzML file");
            f->reader = ms::open_mzml(p, progress_wrap, &pr);
        } else if (ends_with_ci(p, ".lcd")) {
            require_regular_file(p, "a Shimadzu .lcd file");
            f->reader = ms::open_shimadzu(p, progress_wrap, &pr);
        } else if (ms::is_bruker_d(p, &d)) {
            f->reader = ms::open_bruker(d, sdk_dir ? sdk_dir : "", progress_wrap, &pr);
        } else {
            throw ms::Error("Not a supported data file (Bruker .d, Shimadzu .lcd or mzML): " + p);
        }
        const std::vector<ms::ScanInfo>& scans = f->reader->scans();
        if (scans.empty()) throw ms::Error("No MS scans in " + p);
        // events: those of the reader (Shimadzu: the method's scan events), else the polarities present, positive first
        const std::vector<ms::ScanInfo>* defs = f->reader->events();
        if (defs && !defs->empty()) {
            for (size_t e = 0; e < defs->size(); e++) {
                ms_file::Event ev;
                ev.polarity = (*defs)[e].polarity < 0 ? -1 : 1;
                ev.lo = (*defs)[e].mz_lo;
                ev.hi = (*defs)[e].mz_hi;
                for (size_t i = 0; i < scans.size(); i++) {
                    const ms::ScanInfo& s = scans[i];
                    if (s.event != (int)e) continue;
                    ev.idx.push_back((long)i);
                    ev.rt.push_back(s.rt);
                    ev.tic.push_back(s.tic);
                    ev.bpc.push_back(s.bpc);
                    ev.flags.push_back(s.saturated ? 1 : 0);
                    if (s.profile) f->has_profile = true;
                }
                f->events.push_back(std::move(ev));
            }
        }
        std::vector<int> pols;
        if (f->events.empty())
        for (const ms::ScanInfo& s : scans) {
            const int pp = s.polarity < 0 ? -1 : 1;
            if (std::find(pols.begin(), pols.end(), pp) == pols.end()) pols.push_back(pp);
            if (s.profile) f->has_profile = true;
        }
        std::sort(pols.begin(), pols.end(), [](int a, int b) { return a > b; });
        for (int pp : pols) {  // (empty when the reader defined the events)
            ms_file::Event ev;
            ev.polarity = pp;
            bool have_lo = false, have_hi = false;
            for (size_t i = 0; i < scans.size(); i++) {
                const ms::ScanInfo& s = scans[i];
                if ((s.polarity < 0 ? -1 : 1) != pp) continue;
                ev.idx.push_back((long)i);
                ev.rt.push_back(s.rt);
                ev.tic.push_back(s.tic);
                ev.bpc.push_back(s.bpc);
                ev.flags.push_back(s.saturated ? 1 : 0);
                if (s.mz_lo > 0) { ev.lo = have_lo ? std::min(ev.lo, s.mz_lo) : s.mz_lo; have_lo = true; }
                if (s.mz_hi > 0) { ev.hi = have_hi ? std::max(ev.hi, s.mz_hi) : s.mz_hi; have_hi = true; }
            }
            f->events.push_back(std::move(ev));
        }
        // summary line, as the window shows it
        std::string ev;
        for (size_t k = 0; k < f->events.size(); k++) {
            if (k) ev += ", ";
            ev += f->events[k].polarity > 0 ? "Positive" : "Negative";
        }
        // wording per format, as hrms_data's summary() methods
        const std::string kind = f->reader->kind();
        char buf[256];
        if (kind == "Shimadzu .lcd") {  // lcms_data.LCDFile.summary()
            std::string names;
            for (size_t k = 0; k < f->events.size(); k++) {
                if (k) names += ", ";
                names += f->events[k].polarity > 0 ? "Positive" : "Negative";
                if (f->events.size() > 1) names += " (event " + std::to_string(k + 1) + ")";
            }
            std::snprintf(buf, sizeof buf, "%ld scans, %.1f min, %d scan event(s): ", (long)scans.size(),
                          scans.back().rt, (int)f->events.size());
            f->summary = buf + names;
        } else {
            const bool tsf = kind == "Bruker .d (TSF)", tdf = kind == "Bruker .d (timsTOF)";
            std::snprintf(buf, sizeof buf, "%ld MS %s, %.1f min, ", (long)scans.size(), (tsf || tdf) ? "spectra" : "scans",
                          scans.back().rt);
            const std::string& inst = f->reader->instrument();
            f->summary = (inst.empty() ? "" : inst + ": ") + buf + ev;
            if (f->reader->n_msms()) {
                std::snprintf(buf, sizeof buf, ", %ld MS/MS %s (not shown)", f->reader->n_msms(),
                              tdf ? "frames" : tsf ? "spectra" : "scans");
                f->summary += buf;
            }
            if (tsf) f->summary += std::string(", TSF format, ") + (f->has_profile ? "profile" : "centroid only");
            else if (tdf) f->summary += ", timsTOF format (ion mobility summed), profile";
            else f->summary += f->has_profile ? ", profile" : ", centroid only";
            const std::string& recal = f->reader->recalibration();
            if (!recal.empty()) f->summary += ", recalibrated (" + recal + ")";
        }
        // the file is open: into the registry
        Registry& r = registry();
        std::lock_guard<std::mutex> lk(r.mu);
        f->serial = r.next++;
        r.serials.insert(f->serial);
        r.live.emplace(f.get(), f->serial);
        out = f.release();
        return out;
}

}  // namespace

extern "C" {

MS_API ms_file* ms_open(const char* path, const char* sdk_dir, ms_progress_fn cb, void* user) {
    ms_file* out = nullptr;
    ms::guarded([&] {
        if (!path) throw std::invalid_argument("ms_open: null path");
        const std::string p(path);
        if (p.empty()) throw std::invalid_argument("ms_open: empty path");
        out = open_file(p, nullptr, sdk_dir, cb, user);
    });
    return out;
}

MS_API ms_file* ms_open_arrays(const char* kind, const char* instrument, long n_scans, const double* rt,
                               const int* event, const long long* offsets, const double* mz, const double* it,
                               const double* tic, const double* bpc, int n_events, const int* ev_polarity,
                               const double* ev_lo, const double* ev_hi, long n_msms) {
    ms_file* out = nullptr;
    ms::guarded([&] {
        const std::string k(kind && *kind ? kind : "arrays");
        out = open_file(k, ms::open_arrays(k, instrument ? instrument : "", n_scans, rt, event, offsets, mz, it, tic,
                                           bpc, n_events, ev_polarity, ev_lo, ev_hi, n_msms),
                        nullptr, nullptr, nullptr);
    });
    return out;
}

MS_API ms_file* ms_open_vendor(const char* kind, const char* instrument, long n_msms, const ms_vendor_scan* scans,
                               long n, const ms_vendor_scan* events, int n_events, int common_axis,
                               ms_vendor_fn fn, void* user, ms_progress_fn cb, void* cbuser) {
    ms_file* out = nullptr;
    ms::guarded([&] {
        Progress pr;
        pr.cb = cb;
        pr.user = cbuser;
        const std::string k(kind && *kind ? kind : "vendor file");
        out = open_file(k, ms::open_vendor(kind, instrument, n_msms, scans, n, events, n_events, common_axis != 0, fn,
                                           user, progress_wrap, &pr),
                        nullptr, nullptr, nullptr);
    });
    return out;
}

MS_API void ms_close(ms_file* f) {
    if (!f) return;
    try {
        bool now = false;
        {
            Registry& r = registry();
            std::lock_guard<std::mutex> lk(r.mu);
            auto it = r.live.find(f);
            if (it == r.live.end()) return;     // closed already, or not a file: nothing to do
            r.serials.erase(it->second);
            r.live.erase(it);
            f->closing = true;                  // running calls stop at their next progress tick
            now = f->busy == 0;                 // else the last running call frees it
        }
        if (now) delete f;
    } catch (...) {
    }
}

MS_API const char* ms_file_kind(const ms_file* f) {
    return ms::guarded_value<const char*>("", [&]() -> const char* {
        FileRef ref(f, false);
        return ref ? ref->reader->kind() : "";   // a string literal of the reader
    });
}
MS_API const char* ms_file_summary(const ms_file* f) {
    return ms::guarded_value<const char*>("", [&]() -> const char* {
        FileRef ref(f, false);
        if (!ref) return "";
        Outputs& o = outputs_of(ref.get());
        o.summary = ref->summary;
        return o.summary.c_str();
    });
}
MS_API int ms_file_has_profile(const ms_file* f) {
    return ms::guarded_value<int>(0, [&] {
        FileRef ref(f, false);
        return ref && ref->has_profile ? 1 : 0;
    });
}
MS_API long ms_file_n_scans(const ms_file* f) {
    return ms::guarded_value<long>(0, [&] {
        FileRef ref(f, false);
        return ref ? (long)ref->reader->scans().size() : 0L;
    });
}
MS_API long ms_file_n_msms(const ms_file* f) {
    return ms::guarded_value<long>(0, [&] {
        FileRef ref(f, false);
        return ref ? ref->reader->n_msms() : 0L;
    });
}
MS_API int ms_file_n_events(const ms_file* f) {
    return ms::guarded_value<int>(0, [&] {
        FileRef ref(f, false);
        return ref ? (int)ref->events.size() : 0;
    });
}

MS_API int ms_file_event_polarity(const ms_file* f, int event) {
    return ms::guarded_value<int>(0, [&] {
        FileRef ref(f, false);
        if (!ref || event < 0 || event >= (int)ref->events.size()) return 0;
        return ref->events[event].polarity;
    });
}

MS_API int ms_file_event_mz_range(const ms_file* f, int event, double* lo, double* hi) {
    bool unknown = false;
    const int rc = ms::guarded([&] {
        FileRef ref(f, false);
        ref.require("ms_file_event_mz_range");
        const ms_file::Event& ev = event_of(ref.get(), event);
        if (lo) *lo = ev.lo;
        if (hi) *hi = ev.hi;
        unknown = ev.lo <= 0 && ev.hi <= 0;
    });
    if (rc == MS_OK && unknown) { ms::set_error("acquisition range unknown"); return MS_NOT_FOUND; }
    return rc;
}

MS_API const char* ms_file_recalibration(const ms_file* f) {
    return ms::guarded_value<const char*>("", [&]() -> const char* {
        FileRef ref(f, false);
        if (!ref) return "";
        Outputs& o = outputs_of(ref.get());
        o.recal = ref->reader->recalibration();
        return o.recal.c_str();
    });
}

MS_API long ms_event_scans(const ms_file* f, int event, const double** rt, const double** tic, const double** bpc) {
    if (rt) *rt = nullptr;
    if (tic) *tic = nullptr;
    if (bpc) *bpc = nullptr;
    return ms::guarded_value<long>(0, [&] {
        FileRef ref(f, false);
        if (!ref || event < 0 || event >= (int)ref->events.size()) return 0L;
        const ms_file::Event& ev = ref->events[event];
        Outputs& o = outputs_of(ref.get());
        const size_t ne = ref->events.size();
        if (o.ev_rt.size() < ne) { o.ev_rt.resize(ne); o.ev_tic.resize(ne); o.ev_bpc.resize(ne); }
        o.ev_rt[event] = ev.rt;
        o.ev_tic[event] = ev.tic;
        o.ev_bpc[event] = ev.bpc;
        if (rt) *rt = o.ev_rt[event].data();
        if (tic) *tic = o.ev_tic[event].data();
        if (bpc) *bpc = o.ev_bpc[event].data();
        return (long)ev.idx.size();
    });
}

// Since 2.8: flags of the scans of an event (bit 0: saturated readings left out),
// in the order of ms_event_scans; the array is owned by the calling thread.
MS_API long ms_event_flags(const ms_file* f, int event, const unsigned char** flags) {
    if (flags) *flags = nullptr;
    return ms::guarded_value<long>(0, [&] {
        FileRef ref(f, false);
        if (!ref || event < 0 || event >= (int)ref->events.size()) return 0L;
        static thread_local std::vector<unsigned char> s_flags;
        s_flags = ref->events[event].flags;
        if (flags) *flags = s_flags.data();
        return (long)s_flags.size();
    });
}

MS_API int ms_xic(ms_file* f, int event, double mz, double tol, const double** y, long* n) {
    if (y) *y = nullptr;
    if (n) *n = 0;
    return ms::guarded([&] {
        FileRef ref(f);
        ref.require("ms_xic");
        if (!y || !n) throw std::invalid_argument("ms_xic: null output");
        const ms_file::Event& ev = event_of(ref.get(), event);
        if (!std::isfinite(mz) || !std::isfinite(tol)) throw std::invalid_argument("ms_xic: the m/z and the tolerance must be finite numbers");
        tol = std::fabs(tol);
        ms::Reader& r = *ref->reader;
        const size_t ne = ev.idx.size(), nall = r.scans().size();
        std::vector<double> yy, out;
        bool own;
        {
            std::lock_guard<std::mutex> lk(ref->reader_mu);
            own = r.xic(ev.polarity, mz, tol, yy);
        }
        if (own) {
            if (yy.size() == nall && nall != ne) {
                out.resize(ne);
                for (size_t k = 0; k < ne; k++) out[k] = yy[ev.idx[k]];
            } else if (yy.size() == ne) {
                out.swap(yy);
            } else {
                throw ms::Error("reader chromatogram of the wrong length");
            }
        } else {
            // from the centroids: summed intensity within mz +- tol per scan
            out.assign(ne, 0.0);
            std::vector<double> cm, ci;
            for (size_t k = 0; k < ne; k++) {
                if (ref->closing.load()) throw ms::Cancelled("cancelled: the file was closed");
                {
                    std::lock_guard<std::mutex> lk(ref->reader_mu);
                    r.centroids(ev.idx[k], cm, ci);
                }
                double s = 0;
                const size_t m = std::min(cm.size(), ci.size());
                for (size_t j = 0; j < m; j++)
                    if (std::fabs(cm[j] - mz) <= tol) s += ci[j];
                out[k] = s;
            }
        }
        Outputs& o = outputs_of(ref.get());
        o.xic_y.swap(out);
        *y = o.xic_y.data();
        *n = (long)ne;
    });
}

// Mass chromatograms of several windows in one pass over the scans (hrms_calib.find_calibrant_segment
// asks for one per calibrant ion: one pass instead of one per window); the same values as n calls of
// ms_xic (each window sums the centroids in their order).
MS_API int ms_xic_multi(ms_file* f, int event, const double* mz, const double* tol, int nw, double* out, long n_out) {
    return ms::guarded([&] {
        FileRef ref(f);
        ref.require("ms_xic_multi");
        if (!mz || !tol || !out || nw < 0) throw std::invalid_argument("ms_xic_multi: null argument");
        const ms_file::Event& ev = event_of(ref.get(), event);
        const size_t ne = ev.idx.size();
        if (n_out != (long)ne) throw std::invalid_argument("ms_xic_multi: the output must hold n_scans values per window");
        std::vector<double> tl(nw);
        for (int w = 0; w < nw; w++) {
            if (!std::isfinite(mz[w]) || !std::isfinite(tol[w]))
                throw std::invalid_argument("ms_xic_multi: the m/z and the tolerances must be finite numbers");
            tl[w] = std::fabs(tol[w]);
        }
        if (nw == 0) return;
        ms::Reader& r = *ref->reader;
        const size_t nall = r.scans().size();
        std::vector<double> yy;
        bool own;
        {
            std::lock_guard<std::mutex> lk(ref->reader_mu);
            own = r.xic(ev.polarity, mz[0], tl[0], yy);
        }
        if (own) {   // readers with chromatograms of their own (timsTOF): one call per window, as ms_xic
            for (int w = 0; w < nw; w++) {
                if (w > 0) {
                    std::lock_guard<std::mutex> lk(ref->reader_mu);
                    if (!r.xic(ev.polarity, mz[w], tl[w], yy)) throw ms::Error("reader chromatogram not available");
                }
                double* o = out + (size_t)w * ne;
                if (yy.size() == nall && nall != ne) {
                    for (size_t k = 0; k < ne; k++) o[k] = yy[ev.idx[k]];
                } else if (yy.size() == ne) {
                    std::copy(yy.begin(), yy.end(), o);
                } else {
                    throw ms::Error("reader chromatogram of the wrong length");
                }
            }
            return;
        }
        std::vector<double> cm, ci;
        for (size_t k = 0; k < ne; k++) {
            if (ref->closing.load()) throw ms::Cancelled("cancelled: the file was closed");
            {
                std::lock_guard<std::mutex> lk(ref->reader_mu);
                r.centroids(ev.idx[k], cm, ci);
            }
            const size_t m = std::min(cm.size(), ci.size());
            bool sorted = true;   // centroids in increasing m/z (as the readers give them): binary search
            for (size_t j = 1; j < m && sorted; j++) sorted = cm[j - 1] <= cm[j];
            for (int w = 0; w < nw; w++) {
                double s = 0;
                size_t j0 = 0, j1 = m;
                if (sorted) {   // the points with |cm - mz| <= tol lie in [j0, j1): the same sum, in order
                    j0 = std::lower_bound(cm.begin(), cm.begin() + m, mz[w] - tl[w] * 1.0000001 - 1e-300) - cm.begin();
                    j1 = std::upper_bound(cm.begin(), cm.begin() + m, mz[w] + tl[w] * 1.0000001 + 1e-300) - cm.begin();
                }
                for (size_t j = j0; j < j1; j++)
                    if (std::fabs(cm[j] - mz[w]) <= tl[w]) s += ci[j];
                out[(size_t)w * ne + k] = s;
            }
        }
    });
}

MS_API int ms_average(ms_file* f, int event, double t0, double t1, const double* bg, int n_bg,
                      const double** mz, const double** it, long* n, long* n_scans,
                      ms_progress_fn cb, void* user) {
    if (mz) *mz = nullptr;
    if (it) *it = nullptr;
    if (n) *n = 0;
    if (n_scans) *n_scans = 0;
    return ms::guarded([&] {
        FileRef ref(f);
        ref.require("ms_average");
        if (!mz || !it || !n || !n_scans) throw std::invalid_argument("ms_average: null output");
        if (n_bg < 0) throw std::invalid_argument("ms_average: negative number of background ranges");
        if (n_bg > 0 && !bg) throw std::invalid_argument("ms_average: null background ranges");
        if (std::isnan(t0) || std::isnan(t1)) throw std::invalid_argument("ms_average: the time range is not a number");
        const ms_file::Event& ev = event_of(ref.get(), event);
        Progress pr;
        pr.cb = cb;
        pr.user = user;
        pr.f = ref.get();
        std::vector<double> amz, ait;
        const std::vector<long> scans = scans_in(ev, t0, t1);
        try {
            if (!scans.empty()) sum_scans(ref.get(), scans, progress_wrap, &pr, amz, ait);
            if (!amz.empty()) {
                const double ns = (double)scans.size();
                for (double& v : ait) v /= ns;
                if (n_bg > 0) {
                    std::vector<long> bscans;
                    for (int k = 0; k < n_bg; k++) {
                        const std::vector<long> s = scans_in(ev, bg[2 * k], bg[2 * k + 1]);
                        bscans.insert(bscans.end(), s.begin(), s.end());
                    }
                    std::sort(bscans.begin(), bscans.end());
                    bscans.erase(std::unique(bscans.begin(), bscans.end()), bscans.end());
                    if (!bscans.empty()) {
                        std::vector<double> bmz, bit;
                        sum_scans(ref.get(), bscans, progress_wrap, &pr, bmz, bit);
                        if (!bmz.empty() && bmz.size() == bit.size()) {
                            const double nb = (double)bscans.size();
                            for (double& v : bit) v /= nb;
                            const long m = (long)amz.size();
                            double* out = ait.data();
                            const double* x = amz.data();
                            double bw;
                            {
                                std::lock_guard<std::mutex> lk(ref->reader_mu);
                                bw = ref->reader->bin_width();
                            }
                            if (bw > 0) {
                                // binned sums (Shimadzu): subtract bin by bin, as LCDFile.average
                                // matches the bin keys (the bin m/z is the intensity weighted mean)
                                std::unordered_map<long long, double> bg_of;
                                bg_of.reserve(bmz.size() * 2);
                                for (size_t j = 0; j < bmz.size(); j++)
                                    if (bit[j] > 0 && std::fabs(bmz[j] / bw) < 9e18) bg_of[(long long)std::nearbyint(bmz[j] / bw)] = bit[j];
                                for (long j = 0; j < m; j++) {
                                    if (!(out[j] > 0) || !(std::fabs(x[j] / bw) < 9e18)) continue;
                                    auto it_ = bg_of.find((long long)std::nearbyint(x[j] / bw));
                                    if (it_ == bg_of.end()) continue;
                                    const double v = out[j] - it_->second;
                                    out[j] = v < 0 ? 0.0 : v;
                                }
                            } else {
                                const int nt = ms::threads();
                                (void)nt;
#ifdef _OPENMP
#pragma omp parallel for schedule(static) num_threads(nt) if (m > 100000)
#endif
                                for (long j = 0; j < m; j++) {
                                    const double v = out[j] - interp0(x[j], bmz, bit);
                                    out[j] = v < 0 ? 0.0 : v;
                                }
                            }
                        }
                    }
                }
            }
        } catch (const ms::Cancelled&) {
            closed_cancel(ref);
            throw;
        }
        Outputs& o = outputs_of(ref.get());
        o.avg_mz.swap(amz);
        o.avg_it.swap(ait);
        *mz = o.avg_mz.data();
        *it = o.avg_it.data();
        *n = (long)o.avg_mz.size();
        *n_scans = (long)scans.size();
    });
}

MS_API int ms_scan(ms_file* f, int event, double t, const double** mz, const double** it, long* n, long* index) {
    if (mz) *mz = nullptr;
    if (it) *it = nullptr;
    if (n) *n = 0;
    return ms::guarded([&] {
        FileRef ref(f);
        ref.require("ms_scan");
        if (!mz || !it || !n) throw std::invalid_argument("ms_scan: null output");
        if (std::isnan(t)) throw std::invalid_argument("ms_scan: the time is not a number");
        const ms_file::Event& ev = event_of(ref.get(), event);
        if (ev.idx.empty()) throw ms::Error("ms_scan: this scan event has no scans");
        size_t best = 0;
        for (size_t k = 1; k < ev.rt.size(); k++)
            if (std::fabs(ev.rt[k] - t) < std::fabs(ev.rt[best] - t)) best = k;
        const std::vector<long> one(1, ev.idx[best]);
        Progress pr;
        pr.f = ref.get();
        std::vector<double> smz, sit;
        try {
            sum_scans(ref.get(), one, progress_wrap, &pr, smz, sit);
        } catch (const ms::Cancelled&) {
            closed_cancel(ref);
            throw;
        }
        Outputs& o = outputs_of(ref.get());
        o.scan_mz.swap(smz);
        o.scan_it.swap(sit);
        *mz = o.scan_mz.data();
        *it = o.scan_it.data();
        *n = (long)o.scan_mz.size();
        if (index) *index = ev.idx[best];
    });
}

MS_API int ms_file_set_bin_width(ms_file* f, double binw) {
    bool unsupported = false;
    const int rc = ms::guarded([&] {
        FileRef ref(f);
        ref.require("ms_file_set_bin_width");
        // (a bin width below 1 ppm of an m/z unit would ask for billions of bins)
        if (!(binw >= 0) || binw > 10 || (binw > 0 && binw < 1e-6))
            throw std::invalid_argument("bin width must be 0 (the default 0.05 Da) or between 0.000001 and 10 Da");
        std::lock_guard<std::mutex> lk(ref->reader_mu);
        if (!ref->reader->set_bin_width(binw)) unsupported = true;
    });
    if (rc == MS_OK && unsupported) { ms::set_error("no binning for this format"); return MS_UNSUPPORTED; }
    return rc;
}

// Shimadzu PDA matrix and sample information: forwarded to the reader (the
// arrays belong to the calling thread, see the top of this file).
MS_API int ms_pda(ms_file* f, const double** wl, long* n_wl, const double** t, long* n_t, const double** a) {
    if (wl) *wl = nullptr;
    if (t) *t = nullptr;
    if (a) *a = nullptr;
    if (n_wl) *n_wl = 0;
    if (n_t) *n_t = 0;
    bool none = false;
    const int rc = ms::guarded([&] {
        FileRef ref(f);
        ref.require("ms_pda");
        if (!wl || !n_wl || !t || !n_t || !a) throw std::invalid_argument("ms_pda: null output");
        std::vector<double> w, tt, aa;
        bool have;
        {
            std::lock_guard<std::mutex> lk(ref->reader_mu);
            have = ref->reader->pda(w, tt, aa);
        }
        if (!have || w.empty() || tt.empty() || aa.size() != w.size() * tt.size()) { none = true; return; }
        Outputs& o = outputs_of(ref.get());
        o.pda_wl.swap(w);
        o.pda_t.swap(tt);
        o.pda_a.swap(aa);
        *wl = o.pda_wl.data(); *n_wl = (long)o.pda_wl.size();
        *t = o.pda_t.data();   *n_t = (long)o.pda_t.size();
        *a = o.pda_a.data();
    });
    if (rc == MS_OK && none) { ms::set_error("no PDA data in this file"); return MS_NOT_FOUND; }
    return rc;
}

MS_API const char* ms_sample_info(ms_file* f) {
    return ms::guarded_value<const char*>("", [&]() -> const char* {
        FileRef ref(f);
        if (!ref) return "";
        std::string s;
        {
            std::lock_guard<std::mutex> lk(ref->reader_mu);
            s = ref->reader->sample_info();
        }
        Outputs& o = outputs_of(ref.get());
        o.info.swap(s);
        return o.info.c_str();
    });
}

}

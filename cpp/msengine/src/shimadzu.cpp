// Shimadzu LabSolutions .lcd files (LCMS Postrun). Port of the Python reference
// lcms_data.LCDFile / read_sample_info and lcms_pda.PDAData, on top of a port
// of the scan decoding of OpenSZRaw (RawReader: Mass Raw Data single quad,
// QTFL RawData centroid, TTFL Raw Data IT-TOF variants of the .lcd container).
//
// The scan stream decoding (ShimadzuScans below) is derived from OpenSZRaw
// (https://github.com/Sigilweaver/OpenSZRaw), Copyright Nathan Riley,
// licensed under the Apache License, Version 2.0; see
// third_party/OPENSZRAW_LICENSE. The source files ported: crates/openszraw/src/
// reader.rs, raw/mod.rs, raw/mass_raw.rs, raw/qtfl.rs, raw/ttfl.rs.
//
// On top of the scans, as lcms_data does:
//   * interleaved scan events (polarity switching) are split by the event
//     number of the scan headers (single quad), else by the repeating pattern
//     of the retention time steps,
//   * polarity and m/z range of each event come from the MS method stream,
//   * the m/z scale is checked against the method's scan range (LCMS-2020
//     files decode with m/z twice too high) and corrected,
//   * saturated readings (>= 2.1e9) are left out,
//   * spectra are kept as float32 like the Python arrays so that all values
//     and the binned averages agree bit for bit.
#include "file.h"
#include "cfb.h"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <set>

namespace ms {
namespace {

// ------------------------------------------------------------- utilities
inline uint16_t rd16(const uint8_t* p) { return (uint16_t)(p[0] | (p[1] << 8)); }
inline uint32_t rd32(const uint8_t* p) { return (uint32_t)p[0] | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24); }
inline uint64_t rd64(const uint8_t* p) { return (uint64_t)rd32(p) | ((uint64_t)rd32(p + 4) << 32); }
inline int32_t rdi32(const uint8_t* p) { return (int32_t)rd32(p); }
inline double rdf64(const uint8_t* p) { uint64_t u = rd64(p); double d; std::memcpy(&d, &u, 8); return d; }

std::string upper(std::string s) { for (char& c : s) c = (char)std::toupper((unsigned char)c); return s; }

// numpy's float32 sum of a contiguous array (pairwise summation of the add
// reduction, loops_utils.h.src), so that the TIC of a scan with saturated
// readings equals float(seg.sum()) of the Python reference.
float np_pairwise_sum_f32(const float* a, long n) {
    if (n < 8) {
        float res = 0.f;
        for (long i = 0; i < n; i++) res += a[i];
        return res;
    } else if (n <= 128) {
        float r[8];
        for (int j = 0; j < 8; j++) r[j] = a[j];
        long i;
        for (i = 8; i < n - (n % 8); i += 8)
            for (int j = 0; j < 8; j++) r[j] += a[i + j];
        float res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
        for (; i < n; i++) res += a[i];
        return res;
    } else {
        long n2 = n / 2;
        n2 -= n2 % 8;
        return np_pairwise_sum_f32(a, n2) + np_pairwise_sum_f32(a + n2, n - n2);
    }
}
inline float np_sum_f32(const float* a, long n) { return 0.f + np_pairwise_sum_f32(a, n); }

// ======================================================================
// OpenSZRaw port: the MS scans of the .lcd container
// ======================================================================
struct RawScan {
    double rt_sec = 0;
    std::vector<double> mz;
    std::vector<float> it;
    bool ms1 = true;
    int event = -1;     // scan event of the header (single quad), -1 = unknown
    int polarity = 0;   // polarity bit of the header (single quad), 0 = unknown
};

struct ShimadzuScans {
    std::string variant;            // "singlequad", "qtfl", "ttfl"
    std::vector<RawScan> scans;     // in file order, failed scans skipped (as OpenSZRaw)

    // raw::detect_variant for .lcd files
    static std::string detect(const CompoundFile& cf) {
        if (cf.exists("TTFL Raw Data")) return "ttfl";
        if (cf.exists("QTFL RawData/Centroid Index")) return "qtfl";
        if (cf.exists("Mass Raw Data/MS Raw Data")) return "singlequad";
        if (cf.exists("TLM Raw Data"))
            throw Error("'TLM Raw Data' storage found - this is a QQQ/triple-quadrupole .lcd file, which is not decoded yet");
        throw Error("none of 'TTFL Raw Data', 'QTFL RawData', 'Mass Raw Data', or 'TLM Raw Data' storage found in .lcd file");
    }

    static std::vector<uint32_t> u32_array(const std::vector<uint8_t>& d) {
        std::vector<uint32_t> v(d.size() / 4);
        for (size_t i = 0; i < v.size(); i++) v[i] = rd32(d.data() + 4 * i);
        return v;
    }

    void load(const CompoundFile& cf, ms_progress_fn cb, void* user) {
        variant = detect(cf);
        if (variant == "singlequad") load_single_quad(cf, cb, user);
        else if (variant == "qtfl") load_qtfl(cf, cb, user);
        else load_ttfl(cf, cb, user);
    }

    // ---- raw/mass_raw.rs: 64 byte header (rt ms at 0x04, u16 peak count at
    // 0x36), then n records of [u16 mz*10][intensity, 1..4 bytes LE]
    static bool parse_single_quad(const uint8_t* b, size_t len, RawScan& out) {
        if (len < 64) return false;
        const uint32_t rt_ms = rd32(b + 4);
        const size_t n = rd16(b + 0x36);
        const size_t payload = len - 64;
        if (n == 0 || payload % n != 0) return false;
        const size_t width = payload / n;
        if (width <= 2 || width > 6) return false;
        const size_t iw = width - 2;
        out.rt_sec = (double)rt_ms / 1000.0;
        // the header also names the scan event (u16 at 0x08: 0, 1, ... in the
        // order of the method) and its polarity (bit 0x04 of the byte at 0x2C,
        // set for negative ions, as in the event table of the MS method)
        out.event = rd16(b + 8);
        out.polarity = (b[0x2C] & 0x04) ? -1 : 1;
        out.mz.resize(n);
        out.it.resize(n);
        const uint8_t* p = b + 64;
        for (size_t i = 0; i < n; i++, p += width) {
            out.mz[i] = (double)rd16(p) / 10.0;
            uint32_t v = 0;
            for (size_t k = 0; k < iw; k++) v |= (uint32_t)p[2 + k] << (8 * k);
            out.it[i] = (float)v;
        }
        return true;
    }

    void load_single_quad(const CompoundFile& cf, ms_progress_fn cb, void* user) {
        const std::vector<uint8_t> si = cf.read("Mass Raw Data/Spectrum Index");
        const std::vector<uint8_t> raw = cf.read("Mass Raw Data/MS Raw Data");
        const std::vector<uint32_t> offsets = u32_array(si);
        const size_t n = offsets.size();
        scans.reserve(n);
        for (size_t i = 0; i < n; i++) {
            if (i % 500 == 0) tick(cb, user, (long)i, (long)n);
            const size_t start = offsets[i], end = i + 1 < n ? offsets[i + 1] : raw.size();
            if (start > end || end > raw.size()) continue;
            RawScan s;
            if (!parse_single_quad(raw.data() + start, end - start, s)) continue;
            scans.push_back(std::move(s));
        }
    }

    // ---- raw/qtfl.rs: Centroid Index records of 24 bytes (offset, cycle, event id),
    // scan: 64 byte header (bpi 0x10, payload size 0x18, intensity width 0x24),
    // n x u64 m/z * 1e12 then n intensities
    void load_qtfl(const CompoundFile& cf, ms_progress_fn cb, void* user) {
        const std::vector<uint8_t> ci = cf.read("QTFL RawData/Centroid Index");
        const std::vector<uint8_t> data = cf.read("QTFL RawData/Centroid Data");
        const std::vector<uint8_t> rt = cf.read("QTFL RawData/Retention Time");
        if (ci.size() % 24) throw Error("Centroid Index stream size is not a multiple of 24");
        if (rt.size() % 12) throw Error("Retention Time stream size is not a multiple of 12");
        const size_t n = ci.size() / 24, nrt = rt.size() / 12;
        for (size_t i = 0; i < n; i++) {
            if (i % 500 == 0) tick(cb, user, (long)i, (long)n);
            const uint8_t* rec = ci.data() + 24 * i;
            const size_t start = rd32(rec), end = i + 1 < n ? rd32(ci.data() + 24 * (i + 1)) : data.size();
            const uint32_t event_id = rd32(rec + 20);
            if (start > end || end > data.size()) continue;
            const uint8_t* b = data.data() + start;
            const size_t len = end - start;
            if (len < 64) continue;
            const size_t payload_size = rd32(b + 0x18), iw = rd32(b + 0x24);
            RawScan s;
            s.rt_sec = (i < nrt ? (double)rd32(rt.data() + 12 * i) : 0.0) / 1000.0;
            s.ms1 = event_id <= 1;
            if (payload_size > 0) {
                if (len - 64 < payload_size) continue;
                if (!(iw == 1 || iw == 2 || iw == 4)) continue;
                const size_t rec_size = 8 + iw;
                if (payload_size % rec_size) continue;
                const size_t np = payload_size / rec_size;
                const uint8_t* pm = b + 64;
                const uint8_t* pi = pm + 8 * np;
                s.mz.resize(np);
                s.it.resize(np);
                for (size_t k = 0; k < np; k++) s.mz[k] = (double)rd64(pm + 8 * k) / 1e12;
                for (size_t k = 0; k < np; k++) {
                    uint32_t v = iw == 1 ? pi[k] : iw == 2 ? rd16(pi + 2 * k) : rd32(pi + 4 * k);
                    s.it[k] = (float)v;
                }
            }
            scans.push_back(std::move(s));
        }
    }

    // ---- raw/ttfl.rs: Data Index of 16 byte subsets (offset, -, entry index, -),
    // scans RLE coded on a time bin index axis, calibrated with the tuning result
    struct Run { uint16_t skip; std::vector<uint16_t> values; };
    static bool decode_rle(const std::vector<uint16_t>& w, size_t start, uint16_t max_run, std::vector<Run>& runs, size_t& end) {
        size_t i = start;
        const size_t n = w.size();
        runs.clear();
        while (i < n) {
            const uint16_t v = w[i];
            if (v == 0x8000) { end = i + 1; return true; }
            if (v > 0x8000 && v <= 0x8000 + max_run) {
                const size_t run_len = v - 0x8000;
                if (i + 2 + run_len > n) return false;
                Run r;
                r.skip = w[i + 1];
                r.values.assign(w.begin() + i + 2, w.begin() + i + 2 + run_len);
                runs.push_back(std::move(r));
                i += 2 + run_len;
            } else {
                return false;
            }
        }
        return false;
    }
    static std::vector<uint16_t> words(const uint8_t* b, size_t len) {
        std::vector<uint16_t> w(len / 2);
        for (size_t i = 0; i < w.size(); i++) w[i] = rd16(b + 2 * i);
        return w;
    }
    static bool find_prefix_end(const uint8_t* p, size_t len, size_t& at) {
        if (len < 4) return false;
        std::vector<Run> runs;
        for (size_t i = 0; i + 1 < len; i += 2) {
            if (p[i + 1] != 0x80) continue;
            if (p[i] > 64) continue;
            const size_t tail = len - i;
            if (tail % 2) continue;
            const std::vector<uint16_t> w = words(p + i, tail);
            size_t end = 0;
            if (decode_rle(w, 0, 4096, runs, end) && end == w.size()) { at = i; return true; }
        }
        return false;
    }
    static bool decode_ttfl(const uint8_t* b, size_t len, std::vector<double>& index_axis, std::vector<float>& it) {
        index_axis.clear();
        it.clear();
        if (len < 64) return false;
        const uint8_t* p = b + 64;
        const size_t plen = len - 64;
        if (plen == 0) return true;
        size_t pe = 0;
        if (!find_prefix_end(p, plen, pe)) return false;
        if ((plen - pe) % 2) return false;
        const std::vector<uint16_t> w = words(p + pe, plen - pe);
        std::vector<Run> runs;
        size_t end = 0;
        if (!decode_rle(w, 0, 4096, runs, end) || end != w.size()) return false;
        uint32_t pos = 0;
        for (const Run& r : runs) {
            pos += r.skip;
            for (uint16_t v : r.values) {
                index_axis.push_back((double)pos);
                it.push_back((float)v);
                pos++;
            }
        }
        return true;
    }
    static bool parse_calibration(const std::vector<uint8_t>& d, double& a, double& b) {
        std::vector<double> masses, times;
        for (size_t i = 0; i < 9; i++) {
            const size_t off = 3022 + 4 * i;
            if (off + 4 > d.size()) break;
            const uint32_t raw = rd32(d.data() + off);
            if (raw == 0) break;
            masses.push_back((double)raw * 1.0e-4);
        }
        if (masses.size() < 3) return false;
        for (size_t i = 0; i < masses.size(); i++) {
            const size_t off = 3150 + 8 * i;
            if (off + 8 > d.size()) return false;
            times.push_back(rdf64(d.data() + off));
        }
        const double n = (double)masses.size();
        double sx = 0, sy = 0, sxx = 0, sxy = 0;
        for (size_t i = 0; i < masses.size(); i++) {
            const double x = std::sqrt(masses[i]);
            sx += x; sy += times[i]; sxx += x * x; sxy += x * times[i];
        }
        const double denom = n * sxx - sx * sx;
        if (std::fabs(denom) < 1e-9) return false;
        a = (n * sxy - sx * sy) / denom;
        b = (sy - a * sx) / n;
        if (!std::isfinite(a) || !std::isfinite(b) || a <= 0.0) return false;
        return true;
    }
    void load_ttfl(const CompoundFile& cf, ms_progress_fn cb, void* user) {
        const std::vector<uint8_t> di = cf.read("TTFL Raw Data/Data Index");
        const std::vector<uint8_t> raw = cf.read("TTFL Raw Data/MS Raw Data");
        const std::vector<uint8_t> rt = cf.read("TTFL Raw Data/Retention Time");
        const std::vector<uint32_t> rt_ms = u32_array(rt);
        const size_t n_full = di.size() / 64, rem = di.size() % 64;
        if (rem != 0 && rem % 16 != 0) throw Error("Data Index stream size is not a whole number of 16-byte subsets");
        struct Subset { size_t entry_i, sub_i; uint32_t offset; };
        std::vector<Subset> subsets;
        auto parse_block = [&](const uint8_t* blk, size_t ns) {
            for (size_t s = 0; s < ns; s++) {
                const uint8_t* sub = blk + 16 * s;
                subsets.push_back(Subset{rd32(sub + 8), s, rd32(sub)});
            }
        };
        for (size_t k = 0; k < n_full; k++) parse_block(di.data() + 64 * k, 4);
        if (rem) parse_block(di.data() + 64 * n_full, rem / 16);
        std::set<uint32_t> offs;
        for (const Subset& s : subsets) offs.insert(s.offset);
        const std::vector<uint32_t> sorted(offs.begin(), offs.end());
        double ca = 0, cbb = 0;
        bool have_cal = false;
        for (const char* path : {"TTFL Tuning/Tuning Result 00", "TTFL Tuning/Tuning Result 01", "TTFL Tuning/Tuning Result 02"}) {
            std::vector<uint8_t> d;
            if (cf.read_opt(path, d) && parse_calibration(d, ca, cbb)) { have_cal = true; break; }
        }
        std::vector<double> axis;
        std::vector<float> it;
        for (size_t i = 0; i < subsets.size(); i++) {
            if (i % 500 == 0) tick(cb, user, (long)i, (long)subsets.size());
            const Subset& s = subsets[i];
            const size_t start = s.offset;
            auto pos = std::lower_bound(sorted.begin(), sorted.end(), s.offset);
            const size_t end = (pos != sorted.end() && pos + 1 != sorted.end()) ? *(pos + 1) : raw.size();
            if (start > end || end > raw.size()) continue;
            if (!decode_ttfl(raw.data() + start, end - start, axis, it)) continue;
            RawScan sc;
            sc.rt_sec = (s.entry_i < rt_ms.size() ? (double)rt_ms[s.entry_i] : 0.0) / 1000.0;
            sc.mz.resize(axis.size());
            for (size_t k = 0; k < axis.size(); k++) {
                const double x = have_cal ? (axis[k] - cbb) / ca : axis[k];
                sc.mz[k] = have_cal ? x * x : axis[k];
            }
            sc.it = it;
            scans.push_back(std::move(sc));
        }
    }
};

// ======================================================================
// lcms_data helpers
// ======================================================================
struct MethodEvent { int polarity = 0; double lo = 0, hi = 0; };   // polarity 0 = unknown

// lcms_data._method_events: the MS method stream, UPD ID="Mass" block, hex blob
// with the event table
std::vector<MethodEvent> method_events(const CompoundFile& cf) {
    std::vector<MethodEvent> events;
    const CompoundFile::Entry* found = nullptr;
    for (const CompoundFile::Entry& e : cf.entries()) {
        if (e.type != 2 || e.name != "GUC.1.METHOD") continue;
        const size_t slash = e.path.rfind('/');
        if (slash == std::string::npos) continue;
        const std::string parent = e.path.substr(0, slash);
        const size_t s2 = parent.rfind('/');
        const std::string pname = s2 == std::string::npos ? parent : parent.substr(s2 + 1);
        if (upper(pname).find("LCMS") == std::string::npos) continue;
        found = &e;
        break;
    }
    if (!found) return events;
    const std::vector<uint8_t> raw = cf.read(*found);
    // decode("utf-16-le", "ignore"): one character per code unit (pairs count once,
    // lone surrogates dropped); non ASCII characters cannot match the patterns
    std::string text;
    text.reserve(raw.size() / 2);
    for (size_t i = 0; i + 1 < raw.size(); i += 2) {
        const uint16_t u = rd16(raw.data() + i);
        if (u >= 0xD800 && u <= 0xDBFF) {
            if (i + 3 < raw.size()) {
                const uint16_t v = rd16(raw.data() + i + 2);
                if (v >= 0xDC00 && v <= 0xDFFF) { text += '\x80'; i += 2; }
            }
            continue;
        }
        if (u >= 0xDC00 && u <= 0xDFFF) continue;
        text += u < 0x80 ? (char)u : '\x80';
    }
    const size_t i0 = text.find("UPD ID=\"Mass\"");
    if (i0 == std::string::npos) return events;
    const std::string win = text.substr(i0, 20000);
    // first run of at least 64 hex digits
    size_t hs = std::string::npos, he = 0;
    for (size_t k = 0; k < win.size();) {
        if (!std::isxdigit((unsigned char)win[k])) { k++; continue; }
        size_t j = k;
        while (j < win.size() && std::isxdigit((unsigned char)win[j])) j++;
        if (j - k >= 64) { hs = k; he = j; break; }
        k = j;
    }
    if (hs == std::string::npos) return events;
    std::vector<uint8_t> b;
    for (size_t k = hs; k + 1 < he; k += 2) {
        auto hv = [](char c) { return c <= '9' ? c - '0' : (std::tolower((unsigned char)c) - 'a' + 10); };
        b.push_back((uint8_t)(hv(win[k]) * 16 + hv(win[k + 1])));
    }
    size_t off = 0;
    while (off + 20 < b.size()) {
        const uint32_t t_end = rd32(&b[off]), t_ev = rd32(&b[off + 4]), zero = rd32(&b[off + 8]);
        const uint32_t a = rd32(&b[off + 12]), c = rd32(&b[off + 16]);
        off += 1;
        if (zero != 0 || !(0 < t_ev && t_ev < 100000) || !(1000 <= t_end && t_end <= 100000000)) continue;
        if (!(0 < a && a < c)) continue;
        const double lo = a / 10000.0, hi = c / 10000.0;
        if (!(5 <= lo && lo < hi && hi <= 10000)) continue;
        const size_t start = off - 1;
        MethodEvent ev;
        ev.lo = lo;
        ev.hi = hi;
        // b.find(b"\xa0\xa5", start + 20, start + 260)
        const size_t s0 = start + 20, s1 = std::min(start + 260, b.size());
        for (size_t j = s0; j + 1 < s1; j++) {
            if (b[j] == 0xA0 && b[j + 1] == 0xA5) {
                if (j > 0) ev.polarity = (b[j - 1] & 0x04) ? -1 : 1;
                break;
            }
        }
        events.push_back(ev);
        off = start + 20;
    }
    return events;
}

// lcms_data._event_period: number of interleaved scan events from the pattern of
// the retention time steps (seconds)
int event_period(const std::vector<double>& rts, int max_period = 6) {
    std::vector<double> d;
    for (size_t i = 1; i < rts.size(); i++) d.push_back(rts[i] - rts[i - 1]);
    if (d.size() < 20) return 1;
    if (d.size() > 2000) d.resize(2000);
    const size_t n = d.size();
    for (int p = 1; p <= max_period; p++) {
        bool all = true;
        for (size_t i = 0; i + p < n; i++) {
            const double a = d[i], bb = d[i + p];
            if (!(std::fabs(a - bb) <= 0.002 + 0.02 * std::fabs(a))) { all = false; break; }
        }
        if (!all) continue;
        if (p == 1) return p;
        double mean = 0;
        for (int k = 0; k < p; k++) mean += d[k];
        mean /= p;
        double var = 0;
        for (int k = 0; k < p; k++) var += (d[k] - mean) * (d[k] - mean);
        if (std::sqrt(var / p) > 1e-4) return p;
    }
    return 1;
}

// cp1252 byte -> UTF-8
void cp1252_to_utf8(const std::string& s, std::string& out) {
    static const uint16_t hi[32] = {0x20AC, 0xFFFD, 0x201A, 0x0192, 0x201E, 0x2026, 0x2020, 0x2021, 0x02C6, 0x2030, 0x0160,
                                    0x2039, 0x0152, 0xFFFD, 0x017D, 0xFFFD, 0xFFFD, 0x2018, 0x2019, 0x201C, 0x201D, 0x2022,
                                    0x2013, 0x2014, 0x02DC, 0x2122, 0x0161, 0x203A, 0x0153, 0xFFFD, 0x017E, 0x0178};
    for (unsigned char c : s) {
        uint32_t u = c < 0x80 ? c : c < 0xA0 ? hi[c - 0x80] : c;
        if (u < 0x80) out += (char)u;
        else if (u < 0x800) { out += (char)(0xC0 | (u >> 6)); out += (char)(0x80 | (u & 0x3F)); }
        else { out += (char)(0xE0 | (u >> 12)); out += (char)(0x80 | ((u >> 6) & 0x3F)); out += (char)(0x80 | (u & 0x3F)); }
    }
}

std::string strip(const std::string& s) {
    size_t a = 0, b = s.size();
    while (a < b && std::isspace((unsigned char)s[a])) a++;
    while (b > a && std::isspace((unsigned char)s[b - 1])) b--;
    return s.substr(a, b - a);
}

// lcms_data._stox
std::string stox(const std::string& v0) {
    const std::string v = strip(v0);
    auto unhex = [](const std::string& h, std::string& out) {
        if (h.size() % 2) return false;
        for (size_t i = 0; i < h.size(); i += 2) {
            if (!std::isxdigit((unsigned char)h[i]) || !std::isxdigit((unsigned char)h[i + 1])) return false;
            auto hv = [](char c) { return c <= '9' ? c - '0' : (std::tolower((unsigned char)c) - 'a' + 10); };
            out += (char)(hv(h[i]) * 16 + hv(h[i + 1]));
        }
        return true;
    };
    if (v.compare(0, 6, "@StoX@") == 0) {
        std::string bytes, out;
        if (!unhex(v.substr(6), bytes)) return "";
        cp1252_to_utf8(bytes, out);
        return strip(out);
    }
    if (v.compare(0, 6, "@FtoX@") == 0) {
        const std::string h = v.substr(6);
        if (h.size() == 8) {
            std::string bytes;
            if (!unhex(h, bytes)) return "";
            uint32_t u = ((uint32_t)(uint8_t)bytes[0] << 24) | ((uint32_t)(uint8_t)bytes[1] << 16) |
                         ((uint32_t)(uint8_t)bytes[2] << 8) | (uint32_t)(uint8_t)bytes[3];
            float f;
            std::memcpy(&f, &u, 4);
            char buf[64];
            std::snprintf(buf, sizeof buf, "%g", (double)f);
            return buf;
        }
        return h;
    }
    return v;
}

// <name>...</name> inside block (first match), as the Python regex with re.S
bool tag(const std::string& block, const char* name, std::string& out) {
    const std::string open = std::string("<") + name + ">", close = std::string("</") + name + ">";
    const size_t a = block.find(open);
    if (a == std::string::npos) return false;
    const size_t b = block.find(close, a + open.size());
    if (b == std::string::npos) return false;
    out = block.substr(a + open.size(), b - a - open.size());
    return true;
}

// lcms_data.read_sample_info as "key\tvalue\n" lines (Python dict order)
std::string sample_info_text(const CompoundFile& cf) {
    std::string out;
    std::vector<uint8_t> raw;
    if (!cf.read_opt("File Property", raw)) return out;
    const std::string txt(raw.begin(), raw.end());   // latin-1: one char per byte
    auto add = [&](const char* key, const std::string& v) { out += key; out += '\t'; out += v; out += '\n'; };
    size_t a = txt.find("<SampleInfo>");
    if (a != std::string::npos) {
        size_t b = txt.find("</SampleInfo>", a);
        if (b != std::string::npos) {
            const std::string blk = txt.substr(a, b + 13 - a);
            static const char* keys[][2] = {{"sample_name", "smpl_name"}, {"sample_id", "smpl_id"}, {"sample_type", "smpl_type"},
                                            {"operator", "operator_name"}, {"vial", "szVialNum"}, {"tray", "tray_name"},
                                            {"inj_vol", "inj_vol"}};
            for (auto& k : keys) {
                std::string v;
                if (tag(blk, k[1], v)) { const std::string s = stox(v); if (!s.empty()) add(k[0], s); }
            }
            std::string slo, shi;
            if (tag(blk, "dwLowDateTime", slo) && tag(blk, "dwHighDateTime", shi)) {
                // (strtoll's end pointers point into these strings: they must live until the checks)
                const std::string s_lo = strip(slo), s_hi = strip(shi);
                char* e1 = nullptr; char* e2 = nullptr;
                const long long lo = std::strtoll(s_lo.c_str(), &e1, 10), hi = std::strtoll(s_hi.c_str(), &e2, 10);
                if (e1 && *e1 == 0 && e2 && *e2 == 0 && !s_lo.empty() && !s_hi.empty() && hi >= 0 && hi <= 0xFFFFFFFFLL) {
                    const long long ft = (long long)(((unsigned long long)hi << 32) | ((unsigned long long)lo & 0xFFFFFFFFULL));
                    // years 1970 to 9999 (outside, datetime.fromtimestamp fails and the key is left out)
                    if (ft >= 116444736000000000LL && ft < 2650467744000000000LL) {
                        // datetime.fromtimestamp((ft - 116444736000000000) / 1e7): local time
                        const double secs = (double)(ft - 116444736000000000LL) / 1e7;
                        double ip = std::floor(secs);
                        long us = (long)std::llround((secs - ip) * 1e6);
                        if (us >= 1000000) { us -= 1000000; ip += 1; }
                        const time_t t = (time_t)ip;
                        std::tm tmv = std::tm();
#ifdef _WIN32
                        const bool tm_ok = localtime_s(&tmv, &t) == 0;
#else
                        const bool tm_ok = localtime_r(&t, &tmv) != nullptr;
#endif
                        char buf[64];
                        std::snprintf(buf, sizeof buf, "%04d-%02d-%02d %02d:%02d:%02d", tmv.tm_year + 1900, tmv.tm_mon + 1,
                                      tmv.tm_mday, tmv.tm_hour, tmv.tm_min, tmv.tm_sec);
                        std::string s = buf;
                        if (us) { std::snprintf(buf, sizeof buf, ".%06ld", us); s += buf; }
                        if (tm_ok) add("acquired", s);
                    }
                }
            }
        }
    }
    a = txt.find("<SampleInfoFile>");
    if (a != std::string::npos) {
        size_t b = txt.find("</SampleInfoFile>", a);
        if (b != std::string::npos) {
            const std::string blk = txt.substr(a, b + 17 - a);
            static const char* keys[][2] = {{"method_file", "methodfile"}, {"batch_file", "batchfile"}, {"tune_file", "tunefile"}};
            for (auto& k : keys) {
                std::string v;
                if (tag(blk, k[1], v)) {
                    std::string s = stox(v);
                    if (s.empty()) continue;
                    for (char& c : s) if (c == '/') c = '\\';
                    const size_t p = s.rfind('\\');
                    add(k[0], p == std::string::npos ? s : s.substr(p + 1));
                }
            }
        }
    }
    std::vector<uint8_t> c;
    if (cf.read_opt("File Comment", c) && !c.empty()) {
        std::string s;
        if (c.size() > 1 && c[1] == 0) {
            for (size_t i = 0; i + 1 < c.size(); i += 2) {
                const uint16_t u = rd16(c.data() + i);
                if (u < 0x80) s += (char)u;
                else if (u < 0x800) { s += (char)(0xC0 | (u >> 6)); s += (char)(0x80 | (u & 0x3F)); }
                else if (u >= 0xD800 && u <= 0xDFFF) { s.clear(); break; }   // UnicodeDecodeError -> ""
                else { s += (char)(0xE0 | (u >> 12)); s += (char)(0x80 | ((u >> 6) & 0x3F)); s += (char)(0x80 | (u & 0x3F)); }
            }
        } else {
            cp1252_to_utf8(std::string(c.begin(), c.end()), s);
        }
        // strip("\x00 \r\n")
        size_t p0 = 0, p1 = s.size();
        auto st = [](char ch) { return ch == 0 || ch == ' ' || ch == '\r' || ch == '\n'; };
        while (p0 < p1 && st(s[p0])) p0++;
        while (p1 > p0 && st(s[p1 - 1])) p1--;
        s = s.substr(p0, p1 - p0);
        if (!s.empty()) add("comment", s);
    }
    return out;
}

// ======================================================================
// lcms_pda.PDAData
// ======================================================================
// variable length signed integer code (top three bits of the first byte = length)
void pda_decode(const std::vector<uint8_t>& raw, size_t nrec, size_t nwl, std::vector<float>& A, size_t& nrows) {
    A.assign(nrec * nwl, 0.f);
    size_t off = 0;
    const size_t n_raw = raw.size();
    std::vector<long long> vals;
    nrows = nrec;
    for (size_t r = 0; r < nrec; r++) {
        if (off + 24 > n_raw || raw[off] != 0x52 || raw[off + 1] != 0x43) { nrows = r; break; }
        const size_t reclen = rd32(&raw[off + 12]);
        size_t pos = off + 24;
        const size_t end = std::min(off + reclen, n_raw);
        vals.clear();
        while (pos + 2 <= end && vals.size() < nwl) {
            const size_t n = rd16(&raw[pos]);
            size_t i = pos + 2;
            const size_t stop = std::min(i + n, n_raw);
            long long acc = 0;
            bool first = true;
            while (i < stop) {
                const unsigned b0 = raw[i];
                const unsigned pre = b0 >> 5;
                long long v;
                if (pre == 0) {
                    v = b0 & 0x1F;
                    if (v & 0x10) v -= 0x20;
                    i += 1;
                } else if (pre == 1) {
                    if (i + 1 >= n_raw) break;
                    v = ((long long)(b0 & 0x1F) << 8) | raw[i + 1];
                    if (v & 0x1000) v -= 0x2000;
                    i += 2;
                } else if (pre == 2) {
                    if (i + 2 >= n_raw) break;
                    v = ((long long)(b0 & 0x1F) << 16) | ((long long)raw[i + 1] << 8) | raw[i + 2];
                    if (v & 0x100000) v -= 0x200000;
                    i += 3;
                } else {
                    if (i + 3 >= n_raw) break;
                    v = ((long long)(b0 & 0x1F) << 24) | ((long long)raw[i + 1] << 16) | ((long long)raw[i + 2] << 8) | raw[i + 3];
                    if (v & 0x10000000) v -= 0x20000000;
                    i += 4;
                }
                acc = first ? v : acc + v;
                first = false;
                vals.push_back(acc);
            }
            pos = stop + 2;
        }
        const size_t m = std::min(vals.size(), nwl);
        float* row = A.data() + r * nwl;
        for (size_t k = 0; k < m; k++) row[k] = (float)vals[k];
        if (reclen == 0) { nrows = r + 1; break; }
        off += reclen;
    }
}

struct PDAMatrix {
    std::vector<double> wl, t, A;   // A: n_t x n_wl row major, mAU
    bool loaded = false;
};

bool read_pda(const CompoundFile& cf, PDAMatrix& out) {
    const std::string base = "PDA 3D Raw Data/";
    const CompoundFile::Entry* e = cf.find(base + "3D Raw Data");
    if (!e || e->type != 2 || e->size == 0) return false;
    const std::vector<uint8_t> raw = cf.read(*e);
    const std::vector<uint8_t> wt = cf.read(base + "Wavelength Table");
    std::vector<uint8_t> st;
    cf.read_opt(base + "Status", st);
    if (wt.size() < 4) return false;
    const size_t nwl = (size_t)std::max(0, rdi32(wt.data()));
    if (wt.size() < 4 + 4 * nwl) throw Error("PDA wavelength table too short");
    std::vector<double> wl(nwl);
    for (size_t i = 0; i < nwl; i++) wl[i] = rdi32(wt.data() + 4 + 4 * i) / 100.0;
    size_t nrec = 0;
    long long t0_ms = 0, interval = 0;
    double wl_lo = 0, wl_hi = 0;
    bool have_range = false;
    if (st.size() >= 28) {
        nrec = (size_t)std::max(0, rdi32(st.data() + 4));
        wl_lo = rdi32(st.data() + 8) / 100.0;
        wl_hi = rdi32(st.data() + 12) / 100.0;
        t0_ms = rdi32(st.data() + 16);
        interval = rdi32(st.data() + 24);
        have_range = true;
    }
    if (!nrec) {
        size_t off = 0;
        while (off + 16 <= raw.size() && raw[off] == 'R' && raw[off + 1] == 'C') {
            nrec++;
            const uint32_t len = rd32(&raw[off + 12]);
            if (len == 0) break;
            off += len;
        }
    }
    if (!interval) interval = 160;
    nrec = std::min(nrec, raw.size() / 24);   // every record has a 24 byte header: a damaged Status cannot ask for more
    if ((double)nrec * (double)nwl > 2e8) throw Error("PDA data too large (damaged file?)");
    std::vector<float> A;
    size_t nrows = 0;
    pda_decode(raw, nrec, nwl, A, nrows);
    std::vector<char> keep(nwl, 1);
    if (have_range && wl_hi > wl_lo)
        for (size_t i = 0; i < nwl; i++) keep[i] = (wl[i] >= wl_lo - 0.01) && (wl[i] <= wl_hi + 0.61);
    out.wl.clear();
    for (size_t i = 0; i < nwl; i++) if (keep[i]) out.wl.push_back(wl[i]);
    const size_t nk = out.wl.size();
    out.A.resize(nrows * nk);
    for (size_t r = 0; r < nrows; r++) {
        const float* row = A.data() + r * nwl;
        double* dst = out.A.data() + r * nk;
        size_t j = 0;
        for (size_t i = 0; i < nwl; i++) if (keep[i]) dst[j++] = (double)(row[i] / 1000.0f);   // micro-AU -> mAU, float32 as numpy
    }
    out.t.resize(nrows);
    for (size_t r = 0; r < nrows; r++) out.t[r] = (double)(t0_ms + (long long)r * interval) / 60000.0;
    out.loaded = true;
    return true;
}

// ======================================================================
// the reader
// ======================================================================
const float SATURATION = 2.1e9f;   // lcms_data.SATURATION on float32 intensities
const long long kMaxBins = 20000000;   // bins of one summed spectrum (lcms_data.MAX_BINS)

class ShimadzuReader : public Reader {
public:
    ShimadzuReader(const std::string& path, ms_progress_fn cb, void* user) : path_(path) {
        CompoundFile cf(path);
        ShimadzuScans sz;
        sz.load(cf, cb, user);
        variant_ = sz.variant;
        if (variant_ == "ttfl") instrument_ = "LCMS-IT-TOF";
        else if (variant_ == "qtfl") instrument_ = "LCMS-9030";
        // ---- LCDFile.__init__: flatten the scans as float32 arrays
        std::vector<RawScan>& raws = sz.scans;
        std::vector<double> rts;
        std::vector<long> ms1;
        for (size_t i = 0; i < raws.size(); i++) {
            if (!raws[i].ms1) { n_msms_++; continue; }
            ms1.push_back((long)i);
        }
        const size_t n = ms1.size();
        rts.resize(n);
        std::vector<double> tic(n), bpc(n);
        offsets_.assign(n + 1, 0);
        size_t total = 0;
        for (size_t k = 0; k < n; k++) total += raws[ms1[k]].mz.size();
        mz_.resize(total);
        it_.resize(total);
        for (size_t k = 0; k < n; k++) {
            const RawScan& s = raws[ms1[k]];
            rts[k] = s.rt_sec;
            const size_t a = offsets_[k], m = s.mz.size();
            offsets_[k + 1] = a + m;
            double sum = 0;
            double top = 0;
            bool any = false;
            for (size_t j = 0; j < m; j++) {
                mz_[a + j] = (float)s.mz[j];
                it_[a + j] = s.it[j];
                sum += (double)s.it[j];
                if (!any || (double)s.it[j] > top) { top = s.it[j]; any = true; }
            }
            tic[k] = sum;      // openszraw effective_tic: the sum of the intensities
            bpc[k] = any ? top : 0.0;
        }
        // scan events: the event number of the scan headers (single quad) when
        // they name 2 to 8 events, else the repeating pattern of the time steps
        // (scan k in event k % n). A scan that cannot be decoded is left out (as
        // OpenSZRaw does): with the pattern alone, one such scan made the whole
        // file a single event of mixed positive and negative scans
        std::vector<int> hev(n), hpol(n);
        for (size_t k = 0; k < n; k++) {
            hev[k] = raws[ms1[k]].event;
            hpol[k] = raws[ms1[k]].polarity;
        }
        raws.clear();
        raws.shrink_to_fit();
        std::vector<int> ids;
        bool by_header = n > 0;
        for (size_t k = 0; k < n && by_header; k++) {
            if (hev[k] < 0) by_header = false;
            else if (std::find(ids.begin(), ids.end(), hev[k]) == ids.end()) {
                ids.push_back(hev[k]);
                if (ids.size() > 8) by_header = false;
            }
        }
        std::vector<int> scan_event(n, 0);
        if (by_header && ids.size() >= 2) {
            std::sort(ids.begin(), ids.end());
            n_events_ = (int)ids.size();
            for (size_t k = 0; k < n; k++)
                scan_event[k] = (int)(std::lower_bound(ids.begin(), ids.end(), hev[k]) - ids.begin());
        } else {
            n_events_ = event_period(rts);
            for (size_t k = 0; k < n; k++) scan_event[k] = (int)(k % n_events_);
        }
        const std::vector<MethodEvent> meth = method_events(cf);
        // m/z scale check against the scan range of the method
        mz_scale_ = 1.0;
        if (!meth.empty() && total > 0) {
            double hi = 0;
            for (const MethodEvent& m : meth) hi = std::max(hi, m.hi);
            float mx = mz_[0];
            for (size_t j = 1; j < total; j++) if (mz_[j] > mx) mx = mz_[j];
            const double k = hi ? (double)mx / hi : 1.0;
            const double rk = std::round(k);
            if (rk >= 2 && std::fabs(k - rk) < 0.05) {
                mz_scale_ = 1.0 / rk;
                const float f = (float)mz_scale_;
                for (size_t j = 0; j < total; j++) mz_[j] *= f;
            }
        }
        // saturated readings: left out
        saturated_.assign(n, false);
        for (size_t k = 0; k < n; k++) {
            const size_t a = offsets_[k], b = offsets_[k + 1];
            bool sat = false;
            for (size_t j = a; j < b; j++) if (it_[j] >= SATURATION) { it_[j] = 0.f; sat = true; }
            if (!sat) continue;
            saturated_[k] = true;
            tic[k] = (double)np_sum_f32(it_.data() + a, (long)(b - a));
            if (b > a) {
                size_t kk = a;
                for (size_t j = a + 1; j < b; j++) if (it_[j] > it_[kk]) kk = j;
                bpc[k] = it_[kk];
            }
        }
        // events: polarity and m/z range of the method; a polarity the method
        // does not give comes from the headers of the event's scans (majority)
        const bool use_meth = (int)meth.size() == n_events_;
        std::vector<long> n_pos(n_events_, 0), n_neg(n_events_, 0);
        for (size_t k = 0; k < n; k++) {
            if (hpol[k] > 0) n_pos[scan_event[k]]++;
            else if (hpol[k] < 0) n_neg[scan_event[k]]++;
        }
        for (int e = 0; e < n_events_; e++) {
            ScanInfo ev;
            ev.polarity = 0;
            if (use_meth) {
                ev.polarity = meth[e].polarity;   // 0 = unknown
                ev.mz_lo = meth[e].lo;
                ev.mz_hi = meth[e].hi;
            }
            if (ev.polarity == 0 && n_pos[e] + n_neg[e] > 0) ev.polarity = n_neg[e] > n_pos[e] ? -1 : 1;
            events_.push_back(ev);
        }
        scans_.resize(n);
        for (size_t k = 0; k < n; k++) {
            ScanInfo& s = scans_[k];
            const int e = scan_event[k];
            s.rt = rts[k] / 60.0;
            s.polarity = events_[e].polarity < 0 ? -1 : 1;
            s.profile = false;
            s.mz_lo = events_[e].mz_lo;
            s.mz_hi = events_[e].mz_hi;
            s.tic = tic[k];
            s.bpc = bpc[k];
            s.event = e;
            s.saturated = saturated_[k];
        }
        // PDA and sample information (small streams; read now, the container is in memory)
        try { read_pda(cf, pda_); } catch (const std::exception&) { pda_.loaded = false; }
        sample_info_ = sample_info_text(cf);
    }

    const char* kind() const override { return "Shimadzu .lcd"; }
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
        const size_t a = offsets_[i], b = offsets_[i + 1];
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
            const size_t a = offsets_[k], b = offsets_[k + 1];
            double s = 0;
            for (size_t j = a; j < b; j++)
                if (std::fabs(mz_[j] - m) <= t) s += (double)it_[j];
            y[k] = s;
        }
        return true;
    }

    // LCDFile._binned + _fill on the 0.05 Da grid: the plain sum per bin
    // (file.cpp divides by the scan count), bin m/z = intensity weighted mean,
    // empty bins of the grid at zero
    bool sum_profiles(const std::vector<long>& idx, Spectrum& out, ms_progress_fn cb, void* user) override {
#ifdef MS_SHIMADZU_LOGBINS   // (test build: the generic 1 ppm LogBins path of file.cpp instead)
        (void)idx; (void)out; (void)cb; (void)user;
        return false;
#endif
        out.profile = false;
        out.mz.clear();
        out.it.clear();
        if (idx.empty()) return true;
        const double binw = binw_ > 0 ? binw_ : 0.05;
        // keys of every point in order, then accumulate per unique key in order of occurrence
        long long kmin = 0, kmax = 0;
        bool any = false;
        for (long i : idx) {
            if (i < 0 || (size_t)i + 1 >= offsets_.size()) throw std::invalid_argument("scan index out of range");
            for (size_t j = offsets_[i]; j < offsets_[i + 1]; j++) {
                const long long key = (long long)std::nearbyint((double)mz_[j] / binw);
                if (!any) { kmin = kmax = key; any = true; }
                else { kmin = std::min(kmin, key); kmax = std::max(kmax, key); }
            }
        }
        if (!any) return true;
        if (kmax - kmin + 1 > kMaxBins) {
            // a tiny bin width over the whole m/z range: gigabytes for one spectrum
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
            for (size_t j = offsets_[i]; j < offsets_[i + 1]; j++) {
                const double m = mz_[j], w = it_[j];
                const size_t b = (size_t)((long long)std::nearbyint(m / binw) - kmin);
                s_i[b] += w;
                s_m[b] += m * w;
                used[b] = 1;
            }
        }
        // _fill: grid from kmin-1 to kmax+1
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

    bool pda(std::vector<double>& wl, std::vector<double>& t, std::vector<double>& A) override {
        if (!pda_.loaded) return false;
        wl = pda_.wl;
        t = pda_.t;
        A = pda_.A;
        return true;
    }
    std::string sample_info() override { return sample_info_; }

private:
    std::string path_, variant_, instrument_, empty_, sample_info_;
    long n_msms_ = 0;
    int n_events_ = 1;
    double mz_scale_ = 1.0;
    double binw_ = 0;                      // bin width of sum_profiles (0 = 0.05 Da)
    std::vector<float> mz_, it_;
    std::vector<size_t> offsets_;
    std::vector<bool> saturated_;
    std::vector<ScanInfo> scans_, events_;
    PDAMatrix pda_;
};

}  // namespace

std::unique_ptr<Reader> open_shimadzu(const std::string& path, ms_progress_fn cb, void* user) {
    return std::unique_ptr<Reader>(new ShimadzuReader(path, cb, user));
}

}  // namespace ms

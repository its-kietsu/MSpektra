// mzML reader (ProteoWizard output): MS1 spectra, 32/64 bit float arrays,
// zlib or no compression, with or without the index of an indexed mzML.
// A small streaming tokenizer of its own (no XML library): the file is read
// in chunks and only the tags and the base64 payloads are looked at.
//
// Memory model (as hrms_data.MzMLFile): the centroids of every scan plus
// TIC/BPC are always kept; profile scans are kept as float32 intensities
// with a shared m/z axis up to kCacheBytes, further scans are read again
// from the file by their byte offset from the index.
#include "file.h"
#include "fileio.h"
#include "miniz.h"
#include <algorithm>
#include <cctype>
#include <cstring>
#include <cstdlib>
#include <memory>
#include <numeric>

namespace ms {
namespace {

// profile scans kept in memory (MSENGINE_MZML_CACHE_MB overrides, for tests)
size_t cache_bytes() {
    const char* e = std::getenv("MSENGINE_MZML_CACHE_MB");
    if (e && *e) return (size_t)std::strtoull(e, nullptr, 10) * 1024 * 1024;
    return size_t(400) * 1024 * 1024;
}

// -------------------------------------------------------------- base64
struct Base64 {
    uint32_t acc = 0;
    int bits = 0;
    static const signed char* table() {
        static const std::vector<signed char> t = [] {
            std::vector<signed char> v(256, -1);
            const char* a = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
            for (int i = 0; i < 64; i++) v[(unsigned char)a[i]] = (signed char)i;
            v[(unsigned char)'-'] = 62;   // url safe variant, harmless
            v[(unsigned char)'_'] = 63;
            return v;
        }();
        return t.data();
    }
    void feed(const char* p, size_t n, std::vector<unsigned char>& out) {
        const signed char* t = table();
        for (size_t i = 0; i < n; i++) {
            const signed char v = t[(unsigned char)p[i]];
            if (v < 0) continue;   // whitespace, '=' padding, anything else
            acc = (acc << 6) | (uint32_t)v;
            bits += 6;
            if (bits >= 8) {
                bits -= 8;
                out.push_back((unsigned char)((acc >> bits) & 0xff));
            }
        }
    }
};

// ----------------------------------------------------------- tokenizer
class Tokenizer {
public:
    Tokenizer(std::FILE* f, int64_t start, size_t chunk = size_t(4) << 20) : f_(f) {
        buf_.resize(chunk);
        seek(start);
    }
    void seek(int64_t off) {
        fseek64(f_, off, SEEK_SET);
        base_ = off;
        pos_ = end_ = 0;
        eof_ = false;
    }
    int64_t offset() const { return base_ + (int64_t)pos_; }   // file position of the next byte

    // the next tag's text between '<' and '>' (comments, CDATA and
    // processing instructions skipped); false at the end of the file
    bool next_tag(std::string& tag) {
        for (;;) {
            // find '<'
            while (pos_ < end_) {
                const char* p = (const char*)std::memchr(buf_.data() + pos_, '<', end_ - pos_);
                if (!p) { pos_ = end_; break; }
                pos_ = (size_t)(p - buf_.data());
                if (complete_tag(tag)) {
                    if (tag.empty()) continue;   // comment, CDATA or processing instruction skipped
                    return true;
                }
                if (eof_) return false;
                break;   // buffer refilled by complete_tag; search again
            }
            if (pos_ >= end_ && !refill()) return false;
        }
    }

    // base64 text up to the next '<' (the </binary> tag), decoded
    void read_binary(std::vector<unsigned char>& out) {
        Base64 b;
        for (;;) {
            if (pos_ >= end_ && !refill()) return;
            const char* p = (const char*)std::memchr(buf_.data() + pos_, '<', end_ - pos_);
            const size_t stop = p ? (size_t)(p - buf_.data()) : end_;
            b.feed(buf_.data() + pos_, stop - pos_, out);
            pos_ = stop;
            if (p) return;
        }
    }

private:
    std::FILE* f_;
    std::vector<char> buf_;
    size_t pos_ = 0, end_ = 0;
    int64_t base_ = 0;
    bool eof_ = false;

    // moves the unread bytes to the front and reads more; false at eof
    bool refill() {
        if (pos_ > 0) {
            std::memmove(buf_.data(), buf_.data() + pos_, end_ - pos_);
            base_ += (int64_t)pos_;
            end_ -= pos_;
            pos_ = 0;
        }
        if (end_ == buf_.size()) buf_.resize(buf_.size() * 2);
        if (eof_) return false;
        const size_t n = std::fread(buf_.data() + end_, 1, buf_.size() - end_, f_);
        if (n == 0) { eof_ = true; return false; }
        end_ += n;
        return true;
    }

    // buf_[pos_] == '<': extracts the tag if it is complete in the buffer,
    // reading more when it is not. Returns true with the tag (pos_ after
    // its '>'); false when the buffer had to be refilled (call again) or
    // the file ended inside the tag.
    bool complete_tag(std::string& tag) {
        for (;;) {
            const size_t start = pos_;
            size_t i = start + 1;
            const char* d = buf_.data();
            // comments, CDATA, processing instructions
            const char* close = nullptr;
            size_t clen = 0;
            if (i + 3 <= end_ && std::memcmp(d + i, "!--", 3) == 0) { close = "-->"; clen = 3; }
            else if (i + 8 <= end_ && std::memcmp(d + i, "![CDATA[", 8) == 0) { close = "]]>"; clen = 3; }
            else if (i < end_ && d[i] == '?') { close = "?>"; clen = 2; }
            if (close) {
                bool found = false;
                for (size_t j = i; j + clen <= end_; j++)
                    if (std::memcmp(d + j, close, clen) == 0) { pos_ = j + clen; found = true; break; }
                if (found) { tag.clear(); return true; }  // skipped: caller loops (empty tag)
            } else {
                char quote = 0;
                for (size_t j = i; j < end_; j++) {
                    const char c = d[j];
                    if (quote) { if (c == quote) quote = 0; }
                    else if (c == '"' || c == '\'') quote = c;
                    else if (c == '>') {
                        tag.assign(d + i, j - i);
                        pos_ = j + 1;
                        return true;
                    }
                }
            }
            // incomplete: read more (keeps the tag start)
            if (!refill()) { pos_ = end_; return false; }
        }
    }
};

// attribute value of a tag, "" when absent
std::string attr(const std::string& tag, const char* name) {
    const size_t ln = std::strlen(name);
    size_t p = 0;
    while ((p = tag.find(name, p)) != std::string::npos) {
        const bool boundary = p > 0 && (tag[p - 1] == ' ' || tag[p - 1] == '\t' || tag[p - 1] == '\n' || tag[p - 1] == '\r');
        size_t q = p + ln;
        while (q < tag.size() && (tag[q] == ' ' || tag[q] == '\t')) q++;
        if (boundary && q < tag.size() && tag[q] == '=') {
            q++;
            while (q < tag.size() && (tag[q] == ' ' || tag[q] == '\t')) q++;
            if (q < tag.size() && (tag[q] == '"' || tag[q] == '\'')) {
                const char quote = tag[q];
                const size_t e = tag.find(quote, q + 1);
                return tag.substr(q + 1, e == std::string::npos ? std::string::npos : e - q - 1);
            }
            return "";
        }
        p += ln;
    }
    return "";
}

// element name of a tag ("spectrum", "/spectrum", ...)
std::string tag_name(const std::string& tag) {
    size_t e = 0;
    while (e < tag.size() && tag[e] != ' ' && tag[e] != '\t' && tag[e] != '\n' && tag[e] != '\r' && !(e > 0 && tag[e] == '/')) e++;
    std::string name = tag.substr(0, e);
    const size_t c = name.find(':');   // namespace prefix
    return c == std::string::npos ? name : name.substr(c + 1);
}

bool starts_with(const std::string& s, const char* p) { return s.compare(0, std::strlen(p), p) == 0; }

// ------------------------------------------------------ one <spectrum>
struct Spec {
    int level = 1;
    double rt = 0;
    int pol = 0;        // 1 = negative (as the Python reference)
    bool prof = false;
    double lo = 0, hi = 0;
    bool has_lo = false, has_hi = false;
    std::vector<double> mz, it;
    bool has_mz = false, has_it = false;
};

std::vector<unsigned char> inflate(const std::vector<unsigned char>& in, size_t guess) {
    // the array length of the file is a hint only (a damaged value must not ask for terabytes):
    // zlib expands at most about 1032 times
    guess = std::min<size_t>(guess, in.size() * 1040 + 1024);
    std::vector<unsigned char> out(std::max<size_t>(guess, in.size() * 4 + 64));
    mz_stream s;
    std::memset(&s, 0, sizeof s);
    if (mz_inflateInit(&s) != MZ_OK) throw Error("zlib: init failed");
    s.next_in = in.data();
    s.avail_in = (unsigned)in.size();
    size_t done = 0;
    for (;;) {
        if (done == out.size()) out.resize(out.size() * 2);
        s.next_out = out.data() + done;
        s.avail_out = (unsigned)std::min<size_t>(out.size() - done, 1u << 30);
        const int rc = mz_inflate(&s, MZ_NO_FLUSH);
        done = (size_t)s.total_out;
        if (rc == MZ_STREAM_END) break;
        if (rc != MZ_OK && !(rc == MZ_BUF_ERROR && s.avail_out == 0)) { mz_inflateEnd(&s); throw Error("zlib: corrupt binary array"); }
        if (s.avail_in == 0 && s.avail_out != 0) break;   // truncated input: keep what we have
    }
    mz_inflateEnd(&s);
    out.resize(done);
    return out;
}

// compressions other than zlib (hrms_data._MZML_OTHER_COMPRESSION): MS-Numpress
// (msconvert --numpressLinear / --numpressPic / --numpressSlof) and the
// truncation and prediction schemes of msconvert --mzTruncation
const char* other_compression(const std::string& acc) {
    static const char* const tab[][2] = {
        {"MS:1002312", "MS-Numpress linear prediction compression"},
        {"MS:1002313", "MS-Numpress positive integer compression"},
        {"MS:1002314", "MS-Numpress short logged float compression"},
        {"MS:1002746", "MS-Numpress linear prediction compression followed by zlib compression"},
        {"MS:1002747", "MS-Numpress positive integer compression followed by zlib compression"},
        {"MS:1002748", "MS-Numpress short logged float compression followed by zlib compression"},
        {"MS:1003088", "truncation, delta prediction and zlib compression"},
        {"MS:1003089", "truncation, linear prediction and zlib compression"},
        {"MS:1003090", "delta prediction and zlib compression"},
        {"MS:1003091", "linear prediction and zlib compression"}};
    for (const auto& t : tab) if (acc == t[0]) return t[1];
    return nullptr;
}

// parses the body of a <spectrum> whose start tag was just read, up to </spectrum>
void parse_spectrum(Tokenizer& tok, const std::string& start_tag, Spec& sp) {
    sp = Spec();
    if (!start_tag.empty() && start_tag.back() == '/') return;   // <spectrum ... /> (no content)
    const std::string dl = attr(start_tag, "defaultArrayLength");
    const size_t default_len = dl.empty() ? 0 : (size_t)std::strtoull(dl.c_str(), nullptr, 10);
    std::string tag;
    bool in_bda = false;
    int kind = 0;        // 1 m/z, 2 intensity
    // value type of the array: 'f' 32 bit float, 'd' 64 bit float, 'i' 32 bit integer, 'l' 64 bit integer
    char type = 'd';
    bool zlib = false, have = false;
    std::string other;   // a compression that is not supported
    size_t arr_len = default_len;
    std::vector<unsigned char> raw;
    while (tok.next_tag(tag)) {
        if (tag.empty()) continue;
        if (tag[0] == '/') {
            const std::string name = tag_name(tag.substr(1));
            if (name == "spectrum") return;
            if (name == "binaryDataArray" && in_bda) {
                if (kind && have) {
                    if (!other.empty())
                        throw Error("This mzML file stores its spectra with " + other + ", which MSpektra cannot read. "
                                    "Convert the data again with ProteoWizard msconvert without that option (zlib "
                                    "compression is fine).");
                    std::vector<double>& dst = kind == 1 ? sp.mz : sp.it;
                    const size_t es = (type == 'f' || type == 'i') ? 4 : 8;
                    if (zlib && !raw.empty()) raw = inflate(raw, arr_len * es);
                    const size_t n = raw.size() / es;
                    dst.resize(n);
                    const unsigned char* src = raw.data();
                    if (type == 'f') {
                        for (size_t i = 0; i < n; i++) { float v; std::memcpy(&v, src + 4 * i, 4); dst[i] = (double)v; }
                    } else if (type == 'i') {
                        for (size_t i = 0; i < n; i++) { int32_t v; std::memcpy(&v, src + 4 * i, 4); dst[i] = (double)v; }
                    } else if (type == 'l') {
                        for (size_t i = 0; i < n; i++) { int64_t v; std::memcpy(&v, src + 8 * i, 8); dst[i] = (double)v; }
                    } else if (n) {
                        if (n) std::memcpy(dst.data(), src, 8 * n);
                    }
                    (kind == 1 ? sp.has_mz : sp.has_it) = true;
                }
                in_bda = false;
            }
            continue;
        }
        const std::string name = tag_name(tag);
        if (name == "cvParam") {
            const std::string acc = attr(tag, "accession");
            if (in_bda) {
                if (acc == "MS:1000514") kind = 1;
                else if (acc == "MS:1000515") kind = 2;
                else if (acc == "MS:1000521") type = 'f';
                else if (acc == "MS:1000523") type = 'd';
                else if (acc == "MS:1000519") type = 'i';
                else if (acc == "MS:1000522") type = 'l';
                else if (acc == "MS:1000574") zlib = true;
                else if (const char* oc = other_compression(acc)) {
                    const std::string nm = attr(tag, "name");
                    other = nm.empty() ? oc : nm;
                }
            } else if (acc == "MS:1000511") {
                const double lv = std::strtod(attr(tag, "value").c_str(), nullptr);
                sp.level = (lv >= -1e9 && lv <= 1e9) ? (int)lv : 0;   // ((int)NaN is undefined: odd levels count as MS/MS)
            } else if (acc == "MS:1000129") {
                sp.pol = 1;
            } else if (acc == "MS:1000128") {
                sp.prof = true;
            } else if (acc == "MS:1000016") {
                const double t = std::strtod(attr(tag, "value").c_str(), nullptr);
                // unit by name or accession (UO:0000010 second, UO:0000031 minute, UO:0000028 millisecond)
                std::string unit = attr(tag, "unitName");
                const std::string ua = attr(tag, "unitAccession");
                std::transform(unit.begin(), unit.end(), unit.begin(), [](unsigned char c) { return (char)std::tolower(c); });
                if (starts_with(unit, "second") || ua == "UO:0000010") sp.rt = t / 60.0;
                else if (starts_with(unit, "millisecond") || ua == "UO:0000028") sp.rt = t / 60000.0;
                else sp.rt = t;
            } else if (acc == "MS:1000501") {
                sp.lo = std::strtod(attr(tag, "value").c_str(), nullptr);
                sp.has_lo = true;
            } else if (acc == "MS:1000500") {
                sp.hi = std::strtod(attr(tag, "value").c_str(), nullptr);
                sp.has_hi = true;
            }
        } else if (name == "binaryDataArray") {
            in_bda = true;
            kind = 0;
            type = 'd';
            zlib = false;
            have = false;
            other.clear();
            raw.clear();
            const std::string al = attr(tag, "arrayLength");
            arr_len = al.empty() ? default_len : (size_t)std::strtoull(al.c_str(), nullptr, 10);
            if (tag.back() == '/') in_bda = false;   // empty element
        } else if (name == "binary" && in_bda) {
            raw.clear();
            if (tag.back() != '/') tok.read_binary(raw);
            have = true;
        }
    }
    throw Error("mzML: file ends inside a <spectrum>");
}

// byte offsets of the spectra of an indexed mzML (file order), empty if none
std::vector<int64_t> index_offsets(std::FILE* f) {
    std::vector<int64_t> out;
    const int64_t size = file_size(f);
    const int64_t tail_len = std::min<int64_t>(size, 4096);
    std::string tail((size_t)tail_len, '\0');
    fseek64(f, size - tail_len, SEEK_SET);
    if (std::fread(&tail[0], 1, (size_t)tail_len, f) != (size_t)tail_len) return out;
    size_t p = tail.find("<indexListOffset>");
    if (p == std::string::npos) return out;
    const int64_t off = (int64_t)std::strtoll(tail.c_str() + p + 17, nullptr, 10);
    if (off <= 0 || off >= size || size - off > ((int64_t)256 << 20)) return out;   // (an index is a few MB; a damaged offset must not load the file)
    std::string idx((size_t)(size - off), '\0');
    fseek64(f, off, SEEK_SET);
    if (std::fread(&idx[0], 1, idx.size(), f) != idx.size()) return out;
    size_t s = idx.find("<index name=\"spectrum\">");
    if (s == std::string::npos) return out;
    size_t e = idx.find("</index>", s);
    if (e == std::string::npos) e = idx.size();
    p = s;
    while ((p = idx.find("<offset", p)) != std::string::npos && p < e) {
        const size_t gt = idx.find('>', p);
        if (gt == std::string::npos || gt >= e) break;
        out.push_back((int64_t)std::strtoll(idx.c_str() + gt + 1, nullptr, 10));
        p = gt + 1;
    }
    return out;
}

// ---------------------------------------------------------------- reader
class MzmlReader : public Reader {
public:
    MzmlReader(const std::string& path, ms_progress_fn cb, void* user) : path_(path) { load(cb, user); }

    const char* kind() const override { return "mzML"; }
    const std::string& instrument() const override { return instrument_; }
    const std::string& recalibration() const override { return empty_; }
    long n_msms() const override { return n_msms_; }
    const std::vector<ScanInfo>& scans() const override { return scans_; }
    bool common_axis() const override { return has_profile_ && common_; }

    void spectrum(long i, Spectrum& out) override {
        if (i < 0 || i >= (long)scans_.size()) throw std::invalid_argument("mzML: scan index out of range");
        if (!scans_[i].profile) {   // centroid scans: the centroids are the spectrum
            centroids(i, out.mz, out.it);
            out.profile = false;
            return;
        }
        const Cached& c = cache_[i];
        if (c.cached) {
            out.mz.assign(c.mz->begin(), c.mz->end());
            out.it.resize(c.it.size());
            for (size_t j = 0; j < c.it.size(); j++) out.it[j] = (double)c.it[j];
            out.profile = scans_[i].profile;
            return;
        }
        if (offsets_[i] >= 0) {
            // no silent switch to the centroids when the scan cannot be read again: a
            // spectrum summed partly from profiles and partly from centroids would be wrong
            try {
                reread(offsets_[i], out);
            } catch (const std::exception& e) {
                throw Error("Scan " + std::to_string(i + 1) + " of " + path_ + " could not be read again from the file "
                            "(moved, changed or deleted while open?): " + e.what());
            }
            const size_t n = std::min(out.mz.size(), out.it.size());
            out.mz.resize(n);
            out.it.resize(n);
            out.profile = scans_[i].profile;
            return;
        }
        centroids(i, out.mz, out.it);
        out.profile = false;
    }

    void centroids(long i, std::vector<double>& mz, std::vector<double>& it) override {
        if (i < 0 || i >= (long)scans_.size()) throw std::invalid_argument("mzML: scan index out of range");
        const size_t a = coff_[i], b = coff_[i + 1];
        mz.assign(cmz_.begin() + a, cmz_.begin() + b);
        it.assign(cit_.begin() + a, cit_.begin() + b);
    }

private:
    struct Cached {
        bool cached = false;
        std::shared_ptr<std::vector<double>> mz;
        std::vector<float> it;
    };
    std::string path_, instrument_, empty_;
    long n_msms_ = 0;
    bool has_profile_ = false, common_ = true;
    std::vector<ScanInfo> scans_;
    std::vector<int64_t> offsets_;        // -1 = unknown
    std::vector<Cached> cache_;
    std::vector<double> cmz_, cit_;       // centroids of all scans, concatenated
    std::vector<size_t> coff_;            // scan i: [coff_[i], coff_[i+1])

    // A scan read again from the file: every call has its own file handle and
    // tokenizer (no state shared between calls or threads; the file is not kept open).
    void reread(int64_t off, Spectrum& out) const {
        std::unique_ptr<std::FILE, FileCloser> f(fopen_utf8(path_, "rb"));
        if (!f) throw Error("cannot open " + path_);
        Tokenizer tok(f.get(), off, size_t(1) << 20);
        std::string tag;
        while (tok.next_tag(tag)) {
            if (tag.empty()) continue;
            if (tag_name(tag) == "spectrum" && tag[0] != '/') {
                Spec sp;
                parse_spectrum(tok, tag, sp);
                if (!sp.has_it) sp.it.assign(sp.mz.size(), 0.0);
                out.mz.swap(sp.mz);
                out.it.swap(sp.it);
                return;
            }
            throw Error("mzML: no <spectrum> at the indexed offset");
        }
        throw Error("mzML: no <spectrum> at the indexed offset");
    }

    void load(ms_progress_fn cb, void* user) {
        std::unique_ptr<std::FILE, FileCloser> f(fopen_utf8(path_, "rb"));
        if (!f) throw Error("cannot open " + path_);
        const int64_t size = file_size(f.get());
        const std::vector<int64_t> index = index_offsets(f.get());
        Tokenizer tok(f.get(), 0);
        std::string tag;
        Spec sp;
        std::vector<double> cm, ci;
        std::shared_ptr<std::vector<double>> shared;   // the m/z axis of the previous scan
        std::shared_ptr<std::vector<double>> last_profile_axis;
        size_t kept = 0;
        const size_t cache_limit = cache_bytes();
        long k = 0;                                    // spectra in file order
        std::vector<double> rts;
        tick(cb, user, 0, (long)(size >> 10));
        while (tok.next_tag(tag)) {
            if (tag.empty() || tag[0] == '/') continue;
            const std::string name = tag_name(tag);
            if (name != "spectrum") continue;
            parse_spectrum(tok, tag, sp);
            const long file_index = k++;
            if (cb && k % 200 == 0) tick(cb, user, (long)(tok.offset() >> 10), (long)(size >> 10));   // kB: long is 32 bit on Windows
            if (sp.level != 1) { n_msms_++; continue; }
            if (!sp.has_it) sp.it.assign(sp.mz.size(), 0.0);
            if (sp.it.size() != sp.mz.size()) {
                const size_t n = std::min(sp.it.size(), sp.mz.size());
                sp.mz.resize(n);
                sp.it.resize(n);
            }
            ScanInfo si;
            si.rt = sp.rt;
            si.polarity = sp.pol ? -1 : 1;
            si.profile = sp.prof;
            si.mz_lo = sp.has_lo ? sp.lo : 0;
            si.mz_hi = sp.has_hi ? sp.hi : 0;
            // centroids (TIC and BPC from them, as the reference does)
            if (sp.prof) {
                centroid(sp.mz.data(), sp.it.data(), (long)sp.mz.size(), 0.002, cm, ci);
            } else {
                cm = sp.mz;
                ci = sp.it;
            }
            double tic = 0, bpc = 0;
            for (double v : ci) { tic += v; if (v > bpc) bpc = v; }
            si.tic = tic;
            si.bpc = bpc;
            scans_.push_back(si);
            offsets_.push_back(file_index < (long)index.size() ? index[file_index] : -1);
            coff_.push_back(cmz_.size());
            cmz_.insert(cmz_.end(), cm.begin(), cm.end());
            cit_.insert(cit_.end(), ci.begin(), ci.end());
            // profile cache (float32 intensities, the m/z axis shared when it repeats)
            Cached c;
            if (shared && shared->size() == sp.mz.size() && std::equal(shared->begin(), shared->end(), sp.mz.begin())) {
                c.mz = shared;
            } else {
                c.mz = std::make_shared<std::vector<double>>(std::move(sp.mz));
                shared = c.mz;
                kept += c.mz->size() * sizeof(double);
            }
            kept += sp.it.size() * sizeof(float);
            if (sp.prof) {
                has_profile_ = true;
                if (last_profile_axis && last_profile_axis != c.mz &&
                    !(last_profile_axis->size() == c.mz->size() && std::equal(last_profile_axis->begin(), last_profile_axis->end(), c.mz->begin())))
                    common_ = false;
                last_profile_axis = c.mz;
            }
            if (offsets_.back() >= 0 && kept > cache_limit) {
                c.cached = false;
                c.mz.reset();
            } else if (sp.prof) {
                c.cached = true;
                c.it.resize(sp.it.size());
                for (size_t j = 0; j < sp.it.size(); j++) c.it[j] = (float)sp.it[j];
            } else {
                c.cached = false;   // centroid scans are served from the centroids
                c.mz.reset();
            }
            cache_.push_back(std::move(c));
        }
        coff_.push_back(cmz_.size());
        if (scans_.empty()) throw Error("No MS1 spectra in " + path_);
        // time order (stable)
        const size_t n = scans_.size();
        std::vector<size_t> order(n);
        std::iota(order.begin(), order.end(), size_t(0));
        std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) { return scans_[a].rt < scans_[b].rt; });
        bool sorted = true;
        for (size_t i = 0; i < n; i++) if (order[i] != i) { sorted = false; break; }
        if (!sorted) {
            std::vector<ScanInfo> s2(n);
            std::vector<int64_t> o2(n);
            std::vector<Cached> c2(n);
            std::vector<double> m2, i2;
            std::vector<size_t> co2(n + 1);
            m2.reserve(cmz_.size());
            i2.reserve(cit_.size());
            for (size_t i = 0; i < n; i++) {
                const size_t j = order[i];
                s2[i] = scans_[j];
                o2[i] = offsets_[j];
                c2[i] = std::move(cache_[j]);
                co2[i] = m2.size();
                m2.insert(m2.end(), cmz_.begin() + coff_[j], cmz_.begin() + coff_[j + 1]);
                i2.insert(i2.end(), cit_.begin() + coff_[j], cit_.begin() + coff_[j + 1]);
            }
            co2[n] = m2.size();
            scans_.swap(s2);
            offsets_.swap(o2);
            cache_.swap(c2);
            cmz_.swap(m2);
            cit_.swap(i2);
            coff_.swap(co2);
        }
        tick(cb, user, (long)(size >> 10), (long)(size >> 10));
    }
};

}  // namespace

std::unique_ptr<Reader> open_mzml(const std::string& path, ms_progress_fn cb, void* user) {
    return std::unique_ptr<Reader>(new MzmlReader(path, cb, user));
}

}  // namespace ms

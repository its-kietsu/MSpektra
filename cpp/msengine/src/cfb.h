// Internal: a small reader of OLE2 / Compound File Binary (CFB) containers,
// enough for Shimadzu LabSolutions .lcd / .qgd files: header, DIFAT, FAT,
// mini FAT, directory tree, stream extraction (version 3 and 4 files).
#pragma once
#include "common.h"
#include "fileio.h"
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <memory>
#include <string>
#include <vector>

namespace ms {

class CompoundFile {
public:
    struct Entry {
        std::string path;     // "Storage/Sub/Stream" (UTF-8), root not included
        std::string name;
        int type = 0;         // 1 storage, 2 stream, 5 root
        uint32_t start = 0;   // first sector (or mini sector)
        uint64_t size = 0;
        uint64_t created = 0, modified = 0;   // FILETIME, 0 = unset
        uint32_t left = 0xFFFFFFFF, right = 0xFFFFFFFF, child = 0xFFFFFFFF;
    };

    // reads the whole file into memory and parses the structures; throws ms::Error
    explicit CompoundFile(const std::string& path) {
        std::unique_ptr<std::FILE, FileCloser> fp(fopen_utf8(path, "rb"));
        if (!fp) throw Error("Cannot open " + path);
        const int64_t n = file_size(fp.get());
        if (n < 512) throw Error("Not an OLE2 compound file: " + path);
        data_.resize((size_t)n);
        if (std::fread(data_.data(), 1, (size_t)n, fp.get()) != (size_t)n) throw Error("Cannot read " + path);
        parse();
    }

    const std::vector<Entry>& entries() const { return entries_; }

    // the entry of a path ("A/B/C", case insensitive as olefile does), or null
    const Entry* find(const std::string& path) const {
        for (const Entry& e : entries_)
            if (iequal(e.path, path)) return &e;
        return nullptr;
    }
    bool exists(const std::string& path) const { return find(path) != nullptr; }
    bool is_stream(const std::string& path) const { const Entry* e = find(path); return e && e->type == 2; }

    // contents of a stream; throws if missing or not a stream
    std::vector<uint8_t> read(const std::string& path) const {
        const Entry* e = find(path);
        if (!e || e->type != 2) throw Error("stream '" + path + "' not found");
        return read(*e);
    }
    // contents of a stream; false when it does not exist
    bool read_opt(const std::string& path, std::vector<uint8_t>& out) const {
        const Entry* e = find(path);
        if (!e || e->type != 2) return false;
        out = read(*e);
        return true;
    }

    std::vector<uint8_t> read(const Entry& e) const {
        std::vector<uint8_t> out;
        if (e.size == 0) return out;
        out.reserve((size_t)std::min<uint64_t>(e.size, data_.size()));   // (a damaged size must not ask for terabytes)
        if (e.type != 5 && e.size < mini_cutoff_) {
            // mini stream: chain in the mini FAT, data inside the root entry's stream
            uint32_t s = e.start;
            uint64_t left = e.size;
            size_t guard = 0;
            while (s < 0xFFFFFFFA && left > 0) {
                if (++guard > minifat_.size() + 1) throw Error("corrupt mini FAT chain");
                const uint64_t off = (uint64_t)s * mini_sector_;
                if (off + mini_sector_ > mini_stream_.size()) throw Error("mini sector out of range");
                const size_t k = (size_t)std::min<uint64_t>(left, mini_sector_);
                out.insert(out.end(), mini_stream_.begin() + off, mini_stream_.begin() + off + k);
                left -= k;
                s = s < minifat_.size() ? minifat_[s] : 0xFFFFFFFE;
            }
        } else {
            read_chain(e.start, e.size, out);
        }
        return out;
    }

private:
    std::vector<uint8_t> data_;
    uint32_t sector_ = 512, mini_sector_ = 64, mini_cutoff_ = 4096;
    std::vector<uint32_t> fat_, minifat_;
    std::vector<uint8_t> mini_stream_;
    std::vector<Entry> entries_;

    static bool iequal(const std::string& a, const std::string& b) {
        if (a.size() != b.size()) return false;
        for (size_t i = 0; i < a.size(); i++)
            if (std::tolower((unsigned char)a[i]) != std::tolower((unsigned char)b[i])) return false;
        return true;
    }
    uint16_t u16(size_t off) const { return (uint16_t)(data_[off] | (data_[off + 1] << 8)); }
    uint32_t u32(size_t off) const {
        return (uint32_t)data_[off] | ((uint32_t)data_[off + 1] << 8) | ((uint32_t)data_[off + 2] << 16) | ((uint32_t)data_[off + 3] << 24);
    }
    uint64_t u64(size_t off) const { return (uint64_t)u32(off) | ((uint64_t)u32(off + 4) << 32); }

    size_t sector_offset(uint32_t s) const { return ((size_t)s + 1) * sector_; }

    void read_chain(uint32_t start, uint64_t size, std::vector<uint8_t>& out) const {
        uint32_t s = start;
        uint64_t left = size;
        size_t guard = 0;
        while (s < 0xFFFFFFFA && left > 0) {
            if (++guard > fat_.size() + 1) throw Error("corrupt FAT chain");
            const size_t off = sector_offset(s);
            if (off + sector_ > data_.size()) throw Error("sector out of range");
            const size_t k = (size_t)std::min<uint64_t>(left, sector_);
            out.insert(out.end(), data_.begin() + off, data_.begin() + off + k);
            left -= k;
            s = s < fat_.size() ? fat_[s] : 0xFFFFFFFE;
        }
    }

    void parse() {
        static const uint8_t magic[8] = {0xD0, 0xCF, 0x11, 0xE0, 0xA1, 0xB1, 0x1A, 0xE1};
        if (std::memcmp(data_.data(), magic, 8) != 0) throw Error("not a valid OLE2/CFBF container");
        const int major = u16(0x1A);
        const int sshift = u16(0x1E), mshift = u16(0x20);
        if (sshift < 7 || sshift > 20 || mshift < 2 || mshift > sshift) throw Error("bad CFB sector size");
        sector_ = 1u << sshift;
        mini_sector_ = 1u << mshift;
        const uint32_t n_dir_sectors = u32(0x28);
        const uint32_t n_fat = u32(0x2C);
        const uint32_t first_dir = u32(0x30);
        mini_cutoff_ = u32(0x38);
        const uint32_t first_minifat = u32(0x3C), n_minifat = u32(0x40);
        const uint32_t first_difat = u32(0x44), n_difat = u32(0x48);
        const uint32_t per_sector = sector_ / 4;
        // DIFAT: 109 entries in the header, then a chain of DIFAT sectors
        std::vector<uint32_t> fat_sectors;
        for (uint32_t i = 0; i < 109 && fat_sectors.size() < n_fat; i++) {
            const uint32_t s = u32(0x4C + 4 * i);
            if (s >= 0xFFFFFFFA) break;
            fat_sectors.push_back(s);
        }
        uint32_t d = first_difat;
        for (uint32_t k = 0; k < n_difat && d < 0xFFFFFFFA && k <= data_.size() / sector_; k++) {
            const size_t off = sector_offset(d);
            if (off + sector_ > data_.size()) throw Error("DIFAT sector out of range");
            for (uint32_t i = 0; i + 1 < per_sector && fat_sectors.size() < n_fat; i++) {
                const uint32_t s = u32(off + 4 * i);
                if (s >= 0xFFFFFFFA) break;
                fat_sectors.push_back(s);
            }
            d = u32(off + 4 * (per_sector - 1));
        }
        fat_.reserve((size_t)fat_sectors.size() * per_sector);
        for (uint32_t s : fat_sectors) {
            const size_t off = sector_offset(s);
            if (off + sector_ > data_.size()) throw Error("FAT sector out of range");
            for (uint32_t i = 0; i < per_sector; i++) fat_.push_back(u32(off + 4 * i));
        }
        // mini FAT
        if (n_minifat > 0 && first_minifat < 0xFFFFFFFA) {
            std::vector<uint8_t> mf;
            read_chain(first_minifat, (uint64_t)n_minifat * sector_, mf);
            minifat_.resize(mf.size() / 4);
            for (size_t i = 0; i < minifat_.size(); i++)
                minifat_[i] = (uint32_t)mf[4 * i] | ((uint32_t)mf[4 * i + 1] << 8) | ((uint32_t)mf[4 * i + 2] << 16) | ((uint32_t)mf[4 * i + 3] << 24);
        }
        // directory
        std::vector<uint8_t> dir;
        const uint64_t dir_size = (major == 4 && n_dir_sectors > 0) ? (uint64_t)n_dir_sectors * sector_
                                                                     : (uint64_t)fat_.size() * sector_;  // v3: unknown, read the chain
        read_chain(first_dir, dir_size, dir);
        const size_t n_entries = dir.size() / 128;
        std::vector<Entry> raw(n_entries);
        for (size_t i = 0; i < n_entries; i++) {
            const uint8_t* p = dir.data() + 128 * i;
            Entry& e = raw[i];
            auto r16 = [&](size_t o) { return (uint16_t)(p[o] | (p[o + 1] << 8)); };
            auto r32 = [&](size_t o) { return (uint32_t)p[o] | ((uint32_t)p[o + 1] << 8) | ((uint32_t)p[o + 2] << 16) | ((uint32_t)p[o + 3] << 24); };
            auto r64 = [&](size_t o) { return (uint64_t)r32(o) | ((uint64_t)r32(o + 4) << 32); };
            const int nlen = r16(0x40);
            e.type = p[0x42];
            // name: UTF-16LE, nlen bytes including the terminator
            int nchars = nlen >= 2 ? (nlen - 2) / 2 : 0;
            if (nchars > 31) nchars = 31;
            for (int c = 0; c < nchars; c++) append_utf8(e.name, r16(2 * c));
            e.left = r32(0x44);
            e.right = r32(0x48);
            e.child = r32(0x4C);
            e.created = r64(0x64);
            e.modified = r64(0x6C);
            e.start = r32(0x74);
            e.size = major == 3 ? (uint64_t)r32(0x78) : r64(0x78);
        }
        if (raw.empty() || raw[0].type != 5) throw Error("CFB: no root directory entry");
        // the root's stream holds the mini stream
        read_chain(raw[0].start, raw[0].size, mini_stream_);
        // tree walk from the root's child
        entries_.clear();
        std::vector<char> seen(n_entries, 0);
        walk(raw, raw[0].child, "", seen, 0);
    }

    void walk(const std::vector<Entry>& raw, uint32_t id, const std::string& prefix, std::vector<char>& seen, int depth) {
        if (id >= raw.size() || seen[id] || depth > 64) return;
        seen[id] = 1;
        const Entry& e = raw[id];
        walk(raw, e.left, prefix, seen, depth + 1);
        if (e.type == 1 || e.type == 2) {
            Entry out = e;
            out.path = prefix.empty() ? e.name : prefix + "/" + e.name;
            entries_.push_back(out);
            if (e.type == 1) walk(raw, e.child, out.path, seen, depth + 1);
        }
        walk(raw, e.right, prefix, seen, depth + 1);
    }

    static void append_utf8(std::string& s, uint32_t c) {
        if (c < 0x80) s += (char)c;
        else if (c < 0x800) { s += (char)(0xC0 | (c >> 6)); s += (char)(0x80 | (c & 0x3F)); }
        else { s += (char)(0xE0 | (c >> 12)); s += (char)(0x80 | ((c >> 6) & 0x3F)); s += (char)(0x80 | (c & 0x3F)); }
    }
};

}  // namespace ms

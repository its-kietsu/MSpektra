// Private helpers of bruker.cpp: UTF-8 paths on both platforms, dynamic
// libraries, SQLite (read only), md5 (for the temp copy name).
#pragma once
#include "common.h"
#include <atomic>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
#include <map>
#include <functional>

#ifdef _WIN32
#  ifndef NOMINMAX
#    define NOMINMAX
#  endif
#  ifndef WIN32_LEAN_AND_MEAN
#    define WIN32_LEAN_AND_MEAN
#  endif
#  include <windows.h>
#  include <sys/types.h>
#  include <sys/stat.h>
#else
#  include <dlfcn.h>
#  include <sys/stat.h>
#  include <dirent.h>
#  include <unistd.h>
#endif

#include "sqlite3.h"

namespace ms {
namespace sys {

#ifdef _WIN32
inline std::wstring to_wide(const std::string& s) {
    if (s.empty()) return L"";
    int n = MultiByteToWideChar(CP_UTF8, 0, s.c_str(), (int)s.size(), nullptr, 0);
    std::wstring w(n > 0 ? n : 0, L'\0');
    if (n > 0) MultiByteToWideChar(CP_UTF8, 0, s.c_str(), (int)s.size(), &w[0], n);
    return w;
}
inline std::string to_utf8(const std::wstring& w) {
    if (w.empty()) return "";
    int n = WideCharToMultiByte(CP_UTF8, 0, w.c_str(), (int)w.size(), nullptr, 0, nullptr, nullptr);
    std::string s(n > 0 ? n : 0, '\0');
    if (n > 0) WideCharToMultiByte(CP_UTF8, 0, w.c_str(), (int)w.size(), &s[0], n, nullptr, nullptr);
    return s;
}
#endif

inline bool is_sep(char c) { return c == '/' || c == '\\'; }

inline std::string join(const std::string& a, const std::string& b) {
    if (a.empty()) return b;
    if (is_sep(a.back())) return a + b;
#ifdef _WIN32
    return a + "\\" + b;
#else
    return a + "/" + b;
#endif
}

inline std::string strip_sep(std::string p) {
    while (p.size() > 1 && is_sep(p.back())) {
#ifdef _WIN32
        if (p.size() == 3 && p[1] == ':') break;   // "C:\"
#endif
        p.pop_back();
    }
    return p;
}

inline std::string dirname(const std::string& p0) {
    std::string p = strip_sep(p0);
    size_t k = p.find_last_of("/\\");
    if (k == std::string::npos) return "";
    if (k == 0) return p.substr(0, 1);
    return p.substr(0, k);
}

inline std::string basename(const std::string& p0) {
    std::string p = strip_sep(p0);
    size_t k = p.find_last_of("/\\");
    return k == std::string::npos ? p : p.substr(k + 1);
}

inline std::string lower(std::string s) {
    for (char& c : s) if (c >= 'A' && c <= 'Z') c = (char)(c - 'A' + 'a');
    return s;
}

inline bool ends_with(const std::string& s, const std::string& suf) {
    return s.size() >= suf.size() && s.compare(s.size() - suf.size(), suf.size(), suf) == 0;
}

struct Stat { bool exists = false, is_dir = false; int64_t size = 0; int64_t mtime = 0; };

inline Stat stat_path(const std::string& p) {
    Stat s;
#ifdef _WIN32
    struct _stat64 st;
    if (_wstat64(to_wide(strip_sep(p)).c_str(), &st) == 0) {
        s.exists = true; s.is_dir = (st.st_mode & _S_IFDIR) != 0; s.size = st.st_size; s.mtime = (int64_t)st.st_mtime;
    }
#else
    struct stat st;
    if (::stat(p.c_str(), &st) == 0) {
        s.exists = true; s.is_dir = S_ISDIR(st.st_mode); s.size = st.st_size; s.mtime = (int64_t)st.st_mtime;
    }
#endif
    return s;
}
inline bool is_file(const std::string& p) { Stat s = stat_path(p); return s.exists && !s.is_dir; }
inline bool is_dir(const std::string& p) { Stat s = stat_path(p); return s.exists && s.is_dir; }

inline std::vector<std::string> listdir(const std::string& d) {
    std::vector<std::string> out;
#ifdef _WIN32
    WIN32_FIND_DATAW fd;
    HANDLE h = FindFirstFileW(to_wide(join(d, "*")).c_str(), &fd);
    if (h == INVALID_HANDLE_VALUE) return out;
    do {
        std::string n = to_utf8(fd.cFileName);
        if (n != "." && n != "..") out.push_back(n);
    } while (FindNextFileW(h, &fd));
    FindClose(h);
#else
    DIR* dir = opendir(d.c_str());
    if (!dir) return out;
    while (dirent* e = readdir(dir)) {
        std::string n = e->d_name;
        if (n != "." && n != "..") out.push_back(n);
    }
    closedir(dir);
#endif
    return out;
}

// removes "/./" segments (the Windows API normalizes itself)
inline std::string normalize_dots(std::string p) {
    for (size_t k; (k = p.find("/./")) != std::string::npos;) p.erase(k, 2);
    return p;
}

inline std::string abspath(const std::string& p) {
#ifdef _WIN32
    std::wstring w = to_wide(p);
    DWORD n = GetFullPathNameW(w.c_str(), 0, nullptr, nullptr);
    if (!n) return p;
    std::wstring buf(n, L'\0');
    DWORD m = GetFullPathNameW(w.c_str(), n, &buf[0], nullptr);
    buf.resize(m);
    return strip_sep(to_utf8(buf));
#else
    if (!p.empty() && p[0] == '/') return strip_sep(normalize_dots(p));
    char cwd[4096];
    if (!getcwd(cwd, sizeof cwd)) return p;
    return strip_sep(normalize_dots(join(cwd, p)));
#endif
}

inline std::string temp_dir() {
#ifdef _WIN32
    wchar_t buf[MAX_PATH + 2];
    DWORD n = GetTempPathW(MAX_PATH + 1, buf);
    if (!n) return "C:\\Temp";
    return strip_sep(to_utf8(std::wstring(buf, n)));
#else
    for (const char* v : {"TMPDIR", "TEMP", "TMP"}) {
        const char* e = getenv(v);
        if (e && *e) return strip_sep(e);
    }
    return "/tmp";
#endif
}

inline bool mkdir_p(const std::string& p) {
    if (p.empty() || is_dir(p)) return true;
    std::string parent = dirname(p);
    if (!parent.empty() && parent != p && !mkdir_p(parent)) return false;
#ifdef _WIN32
    return CreateDirectoryW(to_wide(p).c_str(), nullptr) || GetLastError() == ERROR_ALREADY_EXISTS;
#else
    return ::mkdir(p.c_str(), 0777) == 0 || errno == EEXIST;
#endif
}

inline bool copy_file(const std::string& a, const std::string& b) {
#ifdef _WIN32
    return CopyFileW(to_wide(a).c_str(), to_wide(b).c_str(), FALSE) != 0;
#else
    FILE* in = fopen(a.c_str(), "rb");
    if (!in) return false;
    FILE* out = fopen(b.c_str(), "wb");
    if (!out) { fclose(in); return false; }
    std::vector<char> buf(1 << 20);
    size_t n;
    bool ok = true;
    while ((n = fread(buf.data(), 1, buf.size(), in)) > 0)
        if (fwrite(buf.data(), 1, n, out) != n) { ok = false; break; }
    fclose(in);
    if (fclose(out) != 0) ok = false;
    return ok;
#endif
}

// copies a folder (one level of files plus sub folders, recursively)
inline bool copy_tree(const std::string& a, const std::string& b) {
    if (!mkdir_p(b)) return false;
    for (const std::string& n : listdir(a)) {
        std::string pa = join(a, n), pb = join(b, n);
        if (is_dir(pa)) { if (!copy_tree(pa, pb)) return false; }
        else if (!copy_file(pa, pb)) return false;
    }
    return true;
}

// deletes a folder and everything in it (best effort)
inline void remove_tree(const std::string& p) {
    for (const std::string& n : listdir(p)) {
        std::string q = join(p, n);
        Stat s = stat_path(q);
        if (s.is_dir) {
            remove_tree(q);
        } else {
#ifdef _WIN32
            std::wstring w = to_wide(q);
            SetFileAttributesW(w.c_str(), FILE_ATTRIBUTE_NORMAL);   // (copies of read only files are read only)
            DeleteFileW(w.c_str());
#else
            ::unlink(q.c_str());
#endif
        }
    }
#ifdef _WIN32
    RemoveDirectoryW(to_wide(p).c_str());
#else
    ::rmdir(p.c_str());
#endif
}

// copies a folder to b through a temporary folder renamed to b when the copy is
// complete: b never exists half copied (an interrupted copy would else be taken
// for the copy forever). True when b exists afterwards.
inline bool copy_tree_complete(const std::string& a, const std::string& b) {
    // a name of its own for every copy (two programs or threads may copy at once)
    static std::atomic<unsigned> counter{0};
#ifdef _WIN32
    const unsigned long pid = (unsigned long)GetCurrentProcessId();
#else
    const unsigned long pid = (unsigned long)getpid();
#endif
    const std::string part = b + ".part" + std::to_string(pid) + "_" + std::to_string(counter++);
    if (is_dir(part)) remove_tree(part);
    if (!copy_tree(a, part)) { remove_tree(part); return false; }
#ifdef _WIN32
    const bool moved = MoveFileExW(to_wide(part).c_str(), to_wide(b).c_str(), 0) != 0;
#else
    const bool moved = ::rename(part.c_str(), b.c_str()) == 0;
#endif
    if (!moved) remove_tree(part);   // (another reader may have made the copy meanwhile)
    return is_dir(b);
}

// folder holding this library (dladdr / GetModuleHandleEx on an address in it)
inline std::string module_dir() {
#ifdef _WIN32
    HMODULE h = nullptr;
    if (!GetModuleHandleExW(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS | GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
                            (LPCWSTR)&module_dir, &h)) return "";
    wchar_t buf[4096];
    DWORD n = GetModuleFileNameW(h, buf, 4096);
    if (!n || n >= 4096) return "";
    return dirname(to_utf8(std::wstring(buf, n)));
#else
    Dl_info info;
    if (dladdr((void*)&module_dir, &info) && info.dli_fname && *info.dli_fname) return dirname(abspath(info.dli_fname));
    return "";
#endif
}

// ------------------------------------------------------------ dynamic library
class DynLib {
public:
    DynLib() {}
    explicit DynLib(const std::string& path) { open(path); }
    ~DynLib() { close(); }
    DynLib(const DynLib&) = delete;
    DynLib& operator=(const DynLib&) = delete;
    void open(const std::string& path) {
        close();
        path_ = path;
#ifdef _WIN32
        // the folder of the DLL first, so that its own dependencies are found
        h_ = (void*)LoadLibraryExW(to_wide(path).c_str(), nullptr, LOAD_WITH_ALTERED_SEARCH_PATH);
        if (!h_) {
            DWORD e = GetLastError();
            throw Error("could not load " + path + " (Windows error " + std::to_string((unsigned long)e) + ")");
        }
#else
        h_ = dlopen(path.c_str(), RTLD_NOW | RTLD_LOCAL);
        if (!h_) {
            const char* e = dlerror();
            throw Error("could not load " + path + (e ? std::string(": ") + e : ""));
        }
#endif
    }
    void close() {
        if (!h_) return;
#ifdef _WIN32
        FreeLibrary((HMODULE)h_);
#else
        dlclose(h_);
#endif
        h_ = nullptr;
    }
    void* sym(const char* name) const {
#ifdef _WIN32
        void* p = (void*)GetProcAddress((HMODULE)h_, name);
#else
        void* p = dlsym(h_, name);
#endif
        if (!p) throw Error(std::string("function ") + name + " not found in " + path_);
        return p;
    }
    template <class F> void get(F& fn, const char* name) const { fn = reinterpret_cast<F>(sym(name)); }
    bool loaded() const { return h_ != nullptr; }
    const std::string& path() const { return path_; }
private:
    void* h_ = nullptr;
    std::string path_;
};

// ------------------------------------------------------------------- SQLite
inline std::string file_uri(const std::string& path) {
    std::string a = abspath(path), out = "file://";
#ifdef _WIN32
    for (char& c : a) if (c == '\\') c = '/';
    if (a.size() >= 2 && a[1] == ':') out += "/";    // file:///C:/...
#endif
    static const char* hex = "0123456789ABCDEF";
    for (unsigned char c : a) {
        if (isalnum(c) || c == '-' || c == '.' || c == '_' || c == '~' || c == '/' || c == ':') out += (char)c;
        else { out += '%'; out += hex[c >> 4]; out += hex[c & 15]; }
    }
    return out;
}

class Db {
public:
    explicit Db(const std::string& path) {
        // read only, as a URI (mode=ro); plain read only open as the fallback
        int rc = sqlite3_open_v2((file_uri(path) + "?mode=ro").c_str(), &db_, SQLITE_OPEN_READONLY | SQLITE_OPEN_URI, nullptr);
        if (rc != SQLITE_OK) {
            if (db_) sqlite3_close(db_);
            db_ = nullptr;
            rc = sqlite3_open_v2(path.c_str(), &db_, SQLITE_OPEN_READONLY, nullptr);
        }
        if (rc != SQLITE_OK) {
            std::string msg = db_ ? sqlite3_errmsg(db_) : "out of memory";
            if (db_) sqlite3_close(db_);
            db_ = nullptr;
            throw Error("could not open " + path + ": " + msg);
        }
    }
    ~Db() { if (db_) sqlite3_close(db_); }
    Db(const Db&) = delete;
    Db& operator=(const Db&) = delete;

    // runs sql; row(stmt) for every row
    void query(const std::string& sql, const std::function<void(sqlite3_stmt*)>& row) const {
        sqlite3_stmt* st = nullptr;
        if (sqlite3_prepare_v2(db_, sql.c_str(), -1, &st, nullptr) != SQLITE_OK)
            throw Error(std::string("SQLite: ") + sqlite3_errmsg(db_));
        int rc;
        try {
            while ((rc = sqlite3_step(st)) == SQLITE_ROW) row(st);
        } catch (...) {
            sqlite3_finalize(st);   // (else the connection cannot close and keeps the file open)
            throw;
        }
        std::string err = rc == SQLITE_DONE ? "" : sqlite3_errmsg(db_);
        sqlite3_finalize(st);
        if (!err.empty()) throw Error("SQLite: " + err);
    }
    std::map<std::string, std::string> key_values(const std::string& sql) const {
        std::map<std::string, std::string> m;
        query(sql, [&](sqlite3_stmt* st) {
            m[text(st, 0)] = text(st, 1);
        });
        return m;
    }
    static std::string text(sqlite3_stmt* st, int col) {
        if (sqlite3_column_type(st, col) == SQLITE_NULL) return "";
        if (sqlite3_column_type(st, col) == SQLITE_BLOB) {
            const char* p = (const char*)sqlite3_column_blob(st, col);
            int n = sqlite3_column_bytes(st, col);
            return std::string(p ? p : "", p ? n : 0);
        }
        const unsigned char* p = sqlite3_column_text(st, col);
        return p ? std::string((const char*)p) : "";
    }
    static bool is_null(sqlite3_stmt* st, int col) { return sqlite3_column_type(st, col) == SQLITE_NULL; }
    static double num(sqlite3_stmt* st, int col, double dflt = 0.0) {
        return is_null(st, col) ? dflt : sqlite3_column_double(st, col);
    }
    static int64_t integer(sqlite3_stmt* st, int col, int64_t dflt = 0) {
        return is_null(st, col) ? dflt : (int64_t)sqlite3_column_int64(st, col);
    }
private:
    sqlite3* db_ = nullptr;
};

// ----------------------------------------------------------------------- md5
// (RFC 1321, compact; only used to name the temporary copy of a write
// protected .d folder exactly as the Python version does)
inline std::string md5_hex(const std::string& msg) {
    static const uint32_t K[64] = {
        0xd76aa478, 0xe8c7b756, 0x242070db, 0xc1bdceee, 0xf57c0faf, 0x4787c62a, 0xa8304613, 0xfd469501,
        0x698098d8, 0x8b44f7af, 0xffff5bb1, 0x895cd7be, 0x6b901122, 0xfd987193, 0xa679438e, 0x49b40821,
        0xf61e2562, 0xc040b340, 0x265e5a51, 0xe9b6c7aa, 0xd62f105d, 0x02441453, 0xd8a1e681, 0xe7d3fbc8,
        0x21e1cde6, 0xc33707d6, 0xf4d50d87, 0x455a14ed, 0xa9e3e905, 0xfcefa3f8, 0x676f02d9, 0x8d2a4c8a,
        0xfffa3942, 0x8771f681, 0x6d9d6122, 0xfde5380c, 0xa4beea44, 0x4bdecfa9, 0xf6bb4b60, 0xbebfbc70,
        0x289b7ec6, 0xeaa127fa, 0xd4ef3085, 0x04881d05, 0xd9d4d039, 0xe6db99e5, 0x1fa27cf8, 0xc4ac5665,
        0xf4292244, 0x432aff97, 0xab9423a7, 0xfc93a039, 0x655b59c3, 0x8f0ccc92, 0xffeff47d, 0x85845dd1,
        0x6fa87e4f, 0xfe2ce6e0, 0xa3014314, 0x4e0811a1, 0xf7537e82, 0xbd3af235, 0x2ad7d2bb, 0xeb86d391};
    static const uint32_t R[64] = {7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22, 7, 12, 17, 22,
                                   5, 9, 14, 20, 5, 9, 14, 20, 5, 9, 14, 20, 5, 9, 14, 20,
                                   4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23, 4, 11, 16, 23,
                                   6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21, 6, 10, 15, 21};
    uint32_t h0 = 0x67452301, h1 = 0xefcdab89, h2 = 0x98badcfe, h3 = 0x10325476;
    std::vector<unsigned char> m(msg.begin(), msg.end());
    uint64_t bits = (uint64_t)m.size() * 8;
    m.push_back(0x80);
    while (m.size() % 64 != 56) m.push_back(0);
    for (int i = 0; i < 8; i++) m.push_back((unsigned char)(bits >> (8 * i)));
    for (size_t off = 0; off < m.size(); off += 64) {
        uint32_t w[16];
        for (int i = 0; i < 16; i++)
            w[i] = (uint32_t)m[off + 4 * i] | ((uint32_t)m[off + 4 * i + 1] << 8) | ((uint32_t)m[off + 4 * i + 2] << 16) |
                   ((uint32_t)m[off + 4 * i + 3] << 24);
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
    std::string out;
    static const char* hex = "0123456789abcdef";
    for (uint32_t h : {h0, h1, h2, h3})
        for (int i = 0; i < 4; i++) { unsigned char byte = (unsigned char)(h >> (8 * i)); out += hex[byte >> 4]; out += hex[byte & 15]; }
    return out;
}

}  // namespace sys
}  // namespace ms

// Internal: stdio files opened from UTF-8 paths (wide API on Windows),
// 64 bit seeking.
#pragma once
#include <cstdio>
#include <cstdint>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef _WIN32
#include <windows.h>
#endif

namespace ms {

#ifdef _WIN32
inline std::wstring utf8_to_wide(const std::string& s) {
    int n = MultiByteToWideChar(CP_UTF8, 0, s.c_str(), -1, nullptr, 0);
    std::wstring w(n > 0 ? n : 1, L'\0');
    if (n > 0) MultiByteToWideChar(CP_UTF8, 0, s.c_str(), -1, &w[0], n);
    w.resize(n > 0 ? n - 1 : 0);
    return w;
}
#endif

inline std::FILE* fopen_utf8(const std::string& path, const char* mode) {
#ifdef _WIN32
    std::wstring wm(mode, mode + std::char_traits<char>::length(mode));
    return _wfopen(utf8_to_wide(path).c_str(), wm.c_str());
#else
    return std::fopen(path.c_str(), mode);
#endif
}

// 0: no such path, 1: a regular file, 2: a folder, 3: something else (device, pipe)
inline int path_type(const std::string& path) {
#ifdef _WIN32
    struct _stat64 st;
    if (_wstat64(utf8_to_wide(path).c_str(), &st) != 0) return 0;
    if (st.st_mode & _S_IFDIR) return 2;
    return (st.st_mode & _S_IFREG) ? 1 : 3;
#else
    struct stat st;
    if (::stat(path.c_str(), &st) != 0) return 0;
    if (S_ISDIR(st.st_mode)) return 2;
    return S_ISREG(st.st_mode) ? 1 : 3;
#endif
}

inline int fseek64(std::FILE* f, int64_t off, int whence) {
#ifdef _WIN32
    return _fseeki64(f, off, whence);
#else
    return fseeko(f, (off_t)off, whence);
#endif
}

inline int64_t ftell64(std::FILE* f) {
#ifdef _WIN32
    return _ftelli64(f);
#else
    return (int64_t)ftello(f);
#endif
}

inline int64_t file_size(std::FILE* f) {
    const int64_t cur = ftell64(f);
    fseek64(f, 0, SEEK_END);
    const int64_t n = ftell64(f);
    fseek64(f, cur, SEEK_SET);
    return n;
}

struct FileCloser {
    void operator()(std::FILE* f) const { if (f) std::fclose(f); }
};

}  // namespace ms

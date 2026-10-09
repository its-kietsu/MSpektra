#include "common.h"
#include <atomic>
#include <cstring>
#include <thread>
#include <cstdio>
#ifdef _WIN32
#include <windows.h>
#else
#include <sys/sysinfo.h>
#endif
#ifdef _OPENMP
#include <omp.h>
#endif

namespace ms {
// The message of the last failure in this thread. A plain character array: setting
// it never allocates or throws (it runs inside catch blocks), and a thread_local of
// a trivial type needs no destructor at thread exit (threads of the host program).
static thread_local char g_error[2048];

void set_error(const char* msg) noexcept {
    if (!msg) msg = "";
    std::strncpy(g_error, msg, sizeof g_error - 1);
    g_error[sizeof g_error - 1] = '\0';
}
void set_error(const std::string& msg) noexcept { set_error(msg.c_str()); }
void clear_error() noexcept { g_error[0] = '\0'; }

int error_code_of_current_exception() noexcept {
    try {
        throw;
    } catch (const Cancelled& e) {
        set_error(e.what());
        return MS_CANCELLED;
    } catch (const std::invalid_argument& e) {
        set_error(e.what());
        return MS_BAD_ARG;
    } catch (const std::bad_alloc&) {
        set_error("Not enough memory for this operation (the data or the settings need more memory than this computer can give)");
        return MS_ERROR;
    } catch (const std::length_error&) {
        set_error("Not enough memory for this operation (a size in the data or the settings is beyond what this computer can hold)");
        return MS_ERROR;
    } catch (const std::exception& e) {
        set_error(e.what());
        return MS_ERROR;
    } catch (...) {
        set_error("unknown error");
        return MS_ERROR;
    }
}

// thread count set by ms_set_threads (0 = automatic); read by every OpenMP region,
// written from any thread: atomic
static std::atomic<int> g_threads(0);

static int cores() {
    int n = (int)std::thread::hardware_concurrency();
    return n > 0 ? n : 1;
}
int max_threads() {
    const int m = 4 * cores();
    return m < 8 ? 8 : (m > 256 ? 256 : m);
}
int threads() {
    const int t = g_threads.load(std::memory_order_relaxed);
    if (t > 0) return t;
    const int n = cores();
    return n > 4 ? n - 1 : n;
}
double memory_total() {
#ifdef _WIN32
    MEMORYSTATUSEX st;
    st.dwLength = sizeof st;
    if (GlobalMemoryStatusEx(&st)) return (double)st.ullTotalPhys;
    return 0;
#else
    struct sysinfo si;
    if (sysinfo(&si) == 0) return (double)si.totalram * (double)si.mem_unit;
    return 0;
#endif
}

double memory_budget() {
    const double t = memory_total();
    const double half = t > 0 ? 0.5 * t : 0.0;
    return half > 1e9 ? half : 1e9;
}

std::string fmt_bytes(double b) {
    char buf[64];
    if (b >= 1e9) snprintf(buf, sizeof buf, "%.1f GB", b / 1e9);
    else snprintf(buf, sizeof buf, "%.0f MB", b / 1e6);
    return buf;
}
}  // namespace ms

extern "C" {
MS_API const char* ms_version(void) { return "msengine 0.1"; }
MS_API const char* ms_last_error(void) { return ms::g_error; }
MS_API int ms_set_threads(int n) {
    // more threads than processors only cost memory and time; far more make the
    // OpenMP runtime fail (its team arrays live on the stack): clamped
    const int m = ms::max_threads();
    ms::clear_error();
    ms::g_threads.store(n > 0 ? (n < m ? n : m) : 0, std::memory_order_relaxed);
#ifdef _OPENMP
    omp_set_num_threads(ms::threads());
#endif
    return MS_OK;
}
MS_API int ms_get_threads(void) { return ms::threads(); }
MS_API double ms_memory_budget(void) { return ms::memory_budget(); }
}

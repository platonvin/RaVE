#include <assert.h>
#include <stdio.h>

// Define platform-specific macros
#if defined(_WIN32) || defined(_WIN64)
    #include <windows.h>
    #define WINDOWS
#elif defined(__linux__) || defined(__unix__) || defined(__APPLE__) && defined(__MACH__)
    #include <unistd.h>
    #define UNIX
#else
    #error "Unsupported platform"
#endif

int get_core_count() {
    int num_cores = 0;

#ifdef WINDOWS
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    num_cores = sysinfo.dwNumberOfProcessors;
#elif defined(UNIX)
    num_cores = sysconf(_SC_NPROCESSORS_ONLN);
#else
#error define sysconf(_SC_NPROCESSORS_ONLN) alternative
#endif

    assert(num_cores != 0);

    return num_cores;
}
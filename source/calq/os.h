#ifndef CALQ_OS_H_
#define CALQ_OS_H_

#ifdef _WIN32
    // Defines for Windows (32-bit and 64-bit, this part is common)
    #define OS_WINDOWS
    #ifdef _WIN64
        // Defines for Windows (64-bit)
        #define OS_WINDOWS_64
    #else
        // Defines for Windows (32-bit)
        #define OS_WINDOWS_32
    #endif
#elif __APPLE__
    // Defines for Apple (all Mac OSes and iOSes)
    #define OS_APPLE
#elif __linux__
    // Defines for Linux
    #define OS_LINUX
#else
    #error "Operating system not supported"
#endif

#endif  // CALQ_OS_H_

/** @file os_config.h
 *  @brief This file detects the operation system we are working on and
 *         does the necessary defines.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: what (who)
 */

#ifndef OS_CONFIG_H
#define OS_CONFIG_H

#ifdef _WIN32
    // defines for Windows (32-bit and 64-bit, this part is common)
    #define OS_WINDOWS
    #ifdef _WIN64
        // defines for Windows (64-bit)
        #define OS_WINDOWS_64
    #else
        // defines for Windows (32-bit)
        #define OS_WINDOWS_32
    #endif
#elif __APPLE__
    // defines for Apple (all Mac OSs and iOSs)
    #define OS_APPLE
#elif __linux__
    // defines for Linux
    #define OS_LINUX
#else
    #error "Operating system not supported"
#endif

#endif // OS_CONFIG_H


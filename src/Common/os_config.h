/** @file os_config.h
 *  @brief This file detects the operation system we are working on and
 *         does the necessary defines.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef OS_CONFIG_H
#define OS_CONFIG_H

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

#endif // OS_CONFIG_H


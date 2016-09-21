/** @file log.h
 *  @brief This file contains useful log/debug macros.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef LOG_H
#define LOG_H

#include "Common/helpers.h"

// C-style debug macro
#define DBG
#ifdef DBG
    #define DEBUG(c,...) \
        do { \
            fflush(stdout); \
            fprintf(stderr, \
                    "DEBUG  %s %s:%s:%d: " c"\n", \
                    currentDateAndTime().c_str(), \
                    __FILE__, \
                    __FUNCTION__, \
                    __LINE__, \
                    ##__VA_ARGS__); \
        } while (false)
#else
    #define DEBUG(c,...) do { } while (false)
#endif

// C-style log macro
#define LOG(c,...) \
    do {\
        fflush(stdout); \
        fprintf(stdout, \
                "LOG  %s %s: " c"\n", \
                currentDateAndTime().c_str(), \
                removeFileNameExtension(fileBaseName(std::string(__FILE__))).c_str(), \
                ##__VA_ARGS__); \
    } while (false)

// Usage example: std::cout << ME << "Hello";
#define MODULE_NAME removeFileNameExtension(fileBaseName(std::string(__FILE__)))
#define CURRENT_DATE_AND_TIME currentDateAndTime()
#define COUT_PREFIX std::string("[" + MODULE_NAME + " @ " + CURRENT_DATE_AND_TIME + "] ")
#define ME COUT_PREFIX // for brevity

#endif // LOG_H


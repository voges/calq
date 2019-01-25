//
// Created by fabian on 1/25/19.
//

#ifndef CALQ_LOG_HPP
#define CALQ_LOG_HPP

#include "calqapp/helpers.h"

// C-style debug macro
// #define DBG
#ifdef DBG
#define CALQ_DEBUG(c, ...) \
        do { \
            fflush(stderr); \
            fprintf(stderr, \
                    "DEBUG  %s %s:%s:%d: " c"\n", \
                    calq::currentDateAndTime().c_str(), \
                    __FILE__, \
                    __FUNCTION__, \
                    __LINE__, \
                    ##__VA_ARGS__); \
        } while (false)
#else
    #define CALQ_DEBUG(c, ...) do { } while (false)
#endif

// -----------------------------------------------------------------------------

// C-style log macro
#define CALQ_LOG(c, ...) \
    do {\
        fflush(stdout); \
        fprintf(stdout, \
                "LOG  %s %s: " c"\n", \
                calqapp::currentDateAndTime().c_str(), \
                calqapp::removeFileNameExtension( \
                calqapp::fileBaseName(std::string(__FILE__))).c_str(), \
                ##__VA_ARGS__); \
    } while (false)

// -----------------------------------------------------------------------------

// C-style error macro
#define CALQ_ERROR(c, ...) \
    do {\
        fflush(stderr); \
        fprintf(stderr, \
                "ERROR  %s " c"\n", \
                calqapp::currentDateAndTime().c_str(), \
                ##__VA_ARGS__); \
    } while (false)


#endif //CALQ_LOG_HPP

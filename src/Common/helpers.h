/** @file helpers.h
 *  @brief This file contains generic helper functions.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  2016-09-22: Added LOG and DEBUG macros (voges)
 *  2016-09-21: Added function currentDateAndTime (voges)
 *  2016-07-14: Added functions fileBaseName and removeFileNameExtension (voges)
 */

#ifndef CQ_HELPERS_H
#define CQ_HELPERS_H

#include <string>

namespace cq {

std::string currentDateAndTime(void);
bool fileExists(const std::string &path);
std::string fileBaseName(const std::string &path);
std::string fileNameExtension(const std::string &path);
std::string removeFileNameExtension(const std::string &path);

}

// C-style debug macro
#define DBG
#ifdef DBG
    #define DEBUG(c,...) \
        do { \
            fflush(stderr); \
            fprintf(stderr, \
                    "DEBUG  %s %s:%s:%d: " c"\n", \
                    cq::currentDateAndTime().c_str(), \
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
                cq::currentDateAndTime().c_str(), \
                cq::removeFileNameExtension(cq::fileBaseName(std::string(__FILE__))).c_str(), \
                ##__VA_ARGS__); \
    } while (false)

// C-style error macro
#define ERROR(c,...) \
    do {\
        fflush(stderr); \
        fprintf(stderr, \
                "ERROR  %s " c"\n", \
                cq::currentDateAndTime().c_str(), \
                ##__VA_ARGS__); \
    } while (false)

#endif // CQ_HELPERS_H


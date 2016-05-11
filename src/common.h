/** @file common.h
 *  @brief This file contains common global definitions.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  2016-05-11: added safe debug macro (voges)
 */

#ifndef COMMON_H
#define COMMON_H

#include <string>

// safe debug macro
#define DBG
#ifdef DBG
    #define DEBUG(c,...)\
        do {\
            fprintf(stderr, "%s:%s:%d: " c, __FILE__, __PRETTY_FUNCTION__, \
                    __LINE__, ##__VA_ARGS__); \
        } while (false)
#else
    #define DEBUG(c,...) do { } while (false)
#endif

struct CLIOptions {
    size_t blockSize;
    bool force;
    std::string infileName;
    enum {
        MODE_COMPRESS,
        MODE_DECOMPRESS,
        MODE_INFO 
    } mode;
    std::string outfileName;
};

extern CLIOptions cliOptions;

const char * filenameExtension(const std::string &filename);
bool fileExists(const std::string &filename);
bool yesno(void);

#endif // COMMON_H


/** @file common.h
 *  @brief This file contains common global definitions.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#ifndef COMMON_H
#define COMMON_H

#include <string>

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


/** @file common.cc
 *  @brief This file contains the implementations of common functions.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  2016-05-29: Moved the fileExists function from access() (which
 *              is defined in unistd.h, which is Linux-specific) to
 *              std::ifstream (which is portable). (voges)
 */

#include "common.h"
#include <fstream>
#include <string.h>

const char * filenameExtension(const std::string &filename)
{
    const char *dot = strrchr(filename.c_str(), '.');
    if (!dot || dot == filename.c_str()) { return ""; }
    return (dot + 1);
}

bool fileExists(const std::string &filename)
{
    std::ifstream ifs(filename.c_str());
    return ifs.good();
}

bool yesno(void)
{
    int c = getchar();
    bool yes = c == 'y' || c == 'Y';
    while (c != '\n' && c != EOF) {
        c = getchar();
    }
    return yes;
}


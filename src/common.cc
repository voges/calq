/** @file common.cc
 *  @brief This file contains the implementations of the functions defined
 *         in common.h.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  2016-05-29: Moved the implementation of the fileExists function from
 *              access() (which is defined in unistd.h, which is
 *              Linux-specific) to std::ifstream (which is portable).
 *             (voges)
 */

#include "common.h"
#include <fstream>
#include <string.h>

const char * fileNameExtension(const std::string &fileName)
{
    const char *dot = strrchr(fileName.c_str(), '.');
    if (!dot || dot == fileName.c_str()) { return ""; }
    return (dot + 1);
}

bool fileExists(const std::string &fileName)
{
    std::ifstream ifs(fileName.c_str());
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


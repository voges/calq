/** @file common.cc
 *  @brief This file contains the implementations of common functions.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#include "common.h"
#include <string.h>
#include <unistd.h>

const char * filenameExtension(const std::string &filename)
{
    const char *dot = strrchr(filename.c_str(), '.');
    if (!dot || dot == filename.c_str()) { return ""; }
    return (dot + 1);
}

bool fileExists(const std::string &filename)
{
    return !access(filename.c_str(), F_OK | R_OK);
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


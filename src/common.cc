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

bool fileExists(const std::string &fileName)
{
    std::ifstream ifs(fileName.c_str());
    return ifs.good();
}

std::string fileBaseName(const std::string &path)
{
    std::string const &delims = "/\\";
  return path.substr(path.find_last_of(delims) + 1);
}

std::string fileNameExtension(const std::string &fileName)
{
    if (fileName.find_last_of(".") != std::string::npos) {
        return fileName.substr(fileName.find_last_of(".")+1);
    }
    return "";
}

std::string removeFileNameExtension(const std::string &fileName)
{
  std::string::size_type const p(fileName.find_last_of('.'));
  return p > 0 && p != std::string::npos ? fileName.substr(0, p) : fileName;
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


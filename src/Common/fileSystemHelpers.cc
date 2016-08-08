/** @file fileSystemHelpers.cc
 *  @brief This file contains the implementations of the functions defined
 *         in fileSystemHelpers.h.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  2017-07-14: Added functions fileBaseName and removeFileNameExtension (voges)
 *  2016-05-29: Moved the implementation of the fileExists function from
 *              access() (which is defined in unistd.h, which is
 *              Linux-specific) to std::ifstream (which is portable). (voges)
 */

#include "fileSystemHelpers.h"
#include <fstream>

bool fileExists(const std::string &path)
{
    std::ifstream ifs(path.c_str());
    return ifs.good();
}

std::string fileBaseName(const std::string &path)
{
    std::string const &delims = "/\\";
    return path.substr(path.find_last_of(delims) + 1);
}

std::string fileNameExtension(const std::string &path)
{
    if (path.find_last_of(".") != std::string::npos) {
        return path.substr(path.find_last_of(".")+1);
    }
    return "";
}

std::string removeFileNameExtension(const std::string &path)
{
    std::string::size_type const p(path.find_last_of('.'));
    return p > 0 && p != std::string::npos ? path.substr(0, p) : path;
}


/** @file helpers.cc
 *  @brief This file contains the implementations of the functions defined
 *         in helpers.h.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  2016-09-21: Added function currentDateAndTime (voges)
 *  2016-07-14: Added functions fileBaseName and removeFileNameExtension (voges)
 *  2016-05-29: Moved the implementation of the fileExists function from
 *              access() (which is defined in unistd.h, which is
 *              Linux-specific) to std::ifstream (which is portable). (voges)
 */

#include "Common/helpers.h"
#include "Common/Exceptions.h"
#include "Common/os_config.h"
#include <fstream>
#include <time.h>

std::string cq::currentDateAndTime(void)
{
    time_t rawtime = time(NULL);
    struct tm timeinfo;

    // Convert to UTC ('Zulu') time
#ifdef CQ_OS_WINDOWS
    errno_t err;
    err = gmtime_s(&timeinfo, &rawtime); 
    if (err != 0) {
        throwErrorException("gmtime_s failed");
    }
#else
    timeinfo = gmtime(&rawtime);
    if (timeinfo == NULL) {
        throwErrorException("gmtime failed");
    }
#endif

    // ISO 8601: 2007-04-05T14:30:21Z
    char timeString[] = "yyyy-mm-ddTHH:MM:SSZ";
    if (strftime(timeString, sizeof(timeString), "%Y-%m-%dT%H:%M:%SZ", &timeinfo) == 0) {
        throwErrorException("strftime failed");
    }
    std::string result(timeString);

    return result;
}

bool cq::fileExists(const std::string &path)
{
    if (path.empty() == true) {
        throwErrorException("path is empty");
    }
    std::ifstream ifs(path.c_str());
    return ifs.good();
}

std::string cq::fileBaseName(const std::string &path)
{
    if (path.empty() == true) {
        throwErrorException("path is empty");
    }
    std::string const &delims = "/\\";
    return path.substr(path.find_last_of(delims) + 1);
}

std::string cq::fileNameExtension(const std::string &path)
{
    if (path.empty() == true) {
        throwErrorException("path is empty");
    }
    if (path.find_last_of(".") != std::string::npos) {
        return path.substr(path.find_last_of(".")+1);
    }
    return "";
}

std::string cq::removeFileNameExtension(const std::string &path)
{
    if (path.empty() == true) {
        throwErrorException("path is empty");
    }
    std::string::size_type const p(path.find_last_of('.'));
    return p > 0 && p != std::string::npos ? path.substr(0, p) : path;
}


/** @file misc.cc
 *  @brief This file contains the implementations of the functions defined
 *         in misc.h.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "misc.h"
#include "Common/Exceptions.h"
#include <ctime>

std::string currentDateAndTime(void)
{
    time_t t = std::time(NULL);
    std::tm *ttm = localtime(&t);
    if (ttm == NULL) {
        throwErrorException("localtime failed");
    }

    // ISO 8601: 2007-04-05T14:30:21
    char timeString[] = "yyyy-mm-ddTHH:MM:SS";
    if (strftime(timeString, sizeof(timeString), "%Y-%m-%dT%H:%M:%S", ttm) == 0) {
        throwErrorException("strftime failed");
    }
    std::string result(timeString);

    return result;
}


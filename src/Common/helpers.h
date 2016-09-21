/** @file helpers.h
 *  @brief This file contains generic helper functions.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  2016-09-21: Added function currentDateAndTime (voges)
 *  2016-07-14: Added functions fileBaseName and removeFileNameExtension (voges)
 */

#ifndef HELPERS_H
#define HELPERS_H

#include <string>

std::string currentDateAndTime(void);
bool fileExists(const std::string &path);
std::string fileBaseName(const std::string &path);
std::string fileNameExtension(const std::string &path);
std::string removeFileNameExtension(const std::string &path);

#endif // HELPERS_H


/** @file fileSystemHelpers.h
 *  @brief This file contains helpers related to the file sytem.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  2017-07-14: Added functions fileBaseName and removeFileNameExtension (voges)
 */

#ifndef FILESYSTEMHELPERS_H
#define FILESYSTEMHELPERS_H

#include <string>

bool fileExists(const std::string &path);
std::string fileBaseName(const std::string &path);
std::string fileNameExtension(const std::string &path);
std::string removeFileNameExtension(const std::string &path);

#endif // FILESYSTEMHELPERS_H


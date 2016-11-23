/** @file FASTAFile.h
 *  @brief This file contains the definition of the FASTAFile class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef FASTAFILE_H
#define FASTAFILE_H

#include "IO/File.h"
#include <map>

class FASTAFile : public File {
public:
    FASTAFile(const std::string &path,
              const FASTAFile::Mode &mode = FASTAFile::MODE_READ);
    ~FASTAFile(void);

public:
    std::map<std::string, std::string> references;

private:
    void parse(void);

private:
    char *line;
    size_t lineSize;
};

#endif // FASTAFILE_H


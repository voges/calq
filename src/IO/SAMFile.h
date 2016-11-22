/** @file SAMFile.h
 *  @brief This file contains the definition of the SAMFile class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef SAMFILE_H
#define SAMFILE_H

#include "IO/File.h"
#include <map>

class SAMFile : public File {
public:
    SAMFile(const std::string &path, const char *mode, const size_t &blockSize);
    ~SAMFile(void);

    void readBlock(void);

    std::string header;

private:
    size_t blockSize;
    char *line;
    size_t lineSize;
};

#endif // SAMFILE_H


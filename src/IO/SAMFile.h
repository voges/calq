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
#include "IO/SAMRecord.h"
#include <deque>

class SAMFile : public File {
public:
    SAMFile(const std::string &path,
            const SAMFile::Mode &mode = SAMFile::MODE_READ);
    ~SAMFile(void);

    size_t readBlock(const size_t &blockSize);

public:
    std::deque<SAMRecord> block;
    std::string header;

private:
    char *line;
    size_t lineSize;
};

#endif // SAMFILE_H


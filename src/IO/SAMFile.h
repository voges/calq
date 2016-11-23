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

    void readBlock(const size_t &blockSize);

public:
    class Block {
    public:
        Block(void)
            : records()
            , numMappedRecords(0)
            , numUnmappedRecords(0)
            , numRecords(0) { }
        ~Block(void) { };

        void reset(void)
        {
            records.clear();
            numUnmappedRecords = 0;
            numMappedRecords = 0;
            numRecords = 0;
        }

    public:
        std::deque<SAMRecord> records;
        size_t numMappedRecords;
        size_t numUnmappedRecords;
        size_t numRecords;
    };

    Block currentBlock;
    std::string header;

private:
    char *line;
    size_t lineSize;
};

#endif // SAMFILE_H


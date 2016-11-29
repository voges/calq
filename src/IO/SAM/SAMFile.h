/** @file SAMFile.h
 *  @brief This file contains the definition of the SAMFile class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef CQ_SAMFILE_H
#define CQ_SAMFILE_H

#include "Common/constants.h"
#include "IO/File.h"
#include "IO/SAM/SAMBlock.h"

namespace cq {

class SAMFile : public File {
public:
    SAMFile(const std::string &path, const Mode &mode = MODE_READ);
    ~SAMFile(void);

    size_t nrBlocksRead(void) const;
    size_t nrMappedRecordsRead(void) const;
    size_t nrUnmappedRecordsRead(void) const;
    size_t nrRecordsRead(void) const;
    size_t readBlock(const size_t &blockSize);

public:
    SAMBlock currentBlock;
    std::string header;

private:
    static const size_t LINE_SIZE = sizeof(char) * (1*MB);

private:
    char *line_;
    size_t nrBlocksRead_;
    size_t nrMappedRecordsRead_;
    size_t nrUnmappedRecordsRead_;
};

}

#endif // CQ_SAMFILE_H


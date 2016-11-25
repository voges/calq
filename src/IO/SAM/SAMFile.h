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

    size_t numBlocksRead(void) const;
    size_t numMappedRecordsRead(void) const;
    size_t numUnmappedRecordsRead(void) const;
    size_t numRecordsRead(void) const;
    size_t readBlock(const size_t &blockSize);

public:
    SAMBlock currentBlock;
    std::string header;

private:
    static const size_t LINE_SIZE = sizeof(char) * (1*MB);

private:
    char *m_line;
    size_t m_numBlocksRead;
    size_t m_numMappedRecordsRead;
    size_t m_numUnmappedRecordsRead;
};

}

#endif // CQ_SAMFILE_H


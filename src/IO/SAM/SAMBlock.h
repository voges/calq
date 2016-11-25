/** @file SAMBlock.h
 *  @brief This file contains the definition of the SAMBlock class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef CQ_SAMBLOCK_H
#define CQ_SAMBLOCK_H

#include "IO/SAM/SAMRecord.h"
#include <deque>

namespace cq {

class SAMBlock {
    friend class SAMFile;

public:
    SAMBlock(void);
    ~SAMBlock(void);

    size_t numMappedRecords(void) const;
    size_t numUnmappedRecords(void) const;
    size_t numRecords(void) const;
    void reset(void);

public:
    std::deque<SAMRecord> records;

private:
    size_t m_numMappedRecords;
    size_t m_numUnmappedRecords;
};

}

#endif // CQ_SAMBLOCK_H


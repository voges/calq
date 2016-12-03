/** @file SAMBlock.h
 *  @brief This file contains the definition of the SAMBlock class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#ifndef CALQ_IO_SAM_SAMBLOCK_H_
#define CALQ_IO_SAM_SAMBLOCK_H_

#include <deque>

#include "IO/SAM/SAMRecord.h"

namespace calq {

class SAMBlock {
    friend class SAMFile;

public:
    SAMBlock(void);
    ~SAMBlock(void);

    size_t nrMappedRecords(void) const;
    size_t nrUnmappedRecords(void) const;
    size_t nrRecords(void) const;
    void reset(void);

    std::deque<SAMRecord> records;

private:
    size_t nrMappedRecords_;
    size_t nrUnmappedRecords_;
};

} // namespace calq

#endif // CALQ_IO_SAM_SAMBLOCK_H_


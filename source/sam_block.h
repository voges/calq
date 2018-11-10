/** @file SAMBlock.h
 *  @brief This file contains the definition of the SAMBlock class.
 */

// Copyright 2015-2017 Leibniz Universitaet Hannover

#ifndef CALQ_IO_SAM_SAMBLOCK_H_
#define CALQ_IO_SAM_SAMBLOCK_H_

#include <deque>

#include "sam_record.h"

namespace calq {

class SAMBlock {
    friend class SAMFile;

 public:
    SAMBlock();
    ~SAMBlock();

    size_t nrMappedRecords() const;
    size_t nrUnmappedRecords() const;
    size_t nrRecords() const;
    void reset();

    std::deque<SAMRecord> records;

 private:
    size_t nrMappedRecords_;
    size_t nrUnmappedRecords_;
};

}  // namespace calq

#endif  // CALQ_IO_SAM_SAMBLOCK_H_

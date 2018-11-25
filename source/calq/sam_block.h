#ifndef CALQ_SAM_BLOCK_H_
#define CALQ_SAM_BLOCK_H_

#include <deque>

#include "calq/sam_record.h"

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

}  // namespace calq

#endif  // CALQ_SAM_BLOCK_H_

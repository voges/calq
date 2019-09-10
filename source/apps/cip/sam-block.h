#ifndef CIP_SAM_BLOCK_H_
#define CIP_SAM_BLOCK_H_

#include <deque>
#include "calq/sam-record.h"

namespace cip {

class SamBlock {
    friend class SAMFile;

   public:
    SamBlock();
    size_t nrMappedRecords() const;
    size_t nrUnmappedRecords() const;
    size_t nrRecords() const;
    void reset();

    std::deque<calq::SamRecord> records;

   private:
    size_t nrMappedRecords_;
    size_t nrUnmappedRecords_;
};

}  // namespace cip

#endif  // CIP_SAM_BLOCK_H_

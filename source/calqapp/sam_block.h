#ifndef CALQAPP_SAM_BLOCK_H_
#define CALQAPP_SAM_BLOCK_H_

// -----------------------------------------------------------------------------

#include <deque>

// -----------------------------------------------------------------------------

#include "calqapp/sam_record.h"

// -----------------------------------------------------------------------------

namespace calqapp {

// -----------------------------------------------------------------------------

class SAMBlock
{
    // -------------------------------------------------------------------------

    friend class SAMFile;

    // -------------------------------------------------------------------------

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

// -----------------------------------------------------------------------------

}  // namespace calqapp

// -----------------------------------------------------------------------------

#endif  // CALQAPP_SAM_BLOCK_H_

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

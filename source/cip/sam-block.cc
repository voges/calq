#include "sam-block.h"

// -----------------------------------------------------------------------------

namespace cip {

// -----------------------------------------------------------------------------

SAMBlock::SAMBlock() : records(), nrMappedRecords_(0), nrUnmappedRecords_(0) {}

// -----------------------------------------------------------------------------

SAMBlock::~SAMBlock() = default;

// -----------------------------------------------------------------------------

size_t SAMBlock::nrMappedRecords() const { return nrMappedRecords_; }

// -----------------------------------------------------------------------------

size_t SAMBlock::nrUnmappedRecords() const { return nrUnmappedRecords_; }

// -----------------------------------------------------------------------------

size_t SAMBlock::nrRecords() const {
    return (nrMappedRecords_ + nrUnmappedRecords_);
}

// -----------------------------------------------------------------------------

void SAMBlock::reset() {
    records.clear();
    nrMappedRecords_ = 0;
    nrUnmappedRecords_ = 0;
}

// -----------------------------------------------------------------------------

}  // namespace cip

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

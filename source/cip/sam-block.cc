#include "sam-block.h"

namespace cip {

SamBlock::SamBlock() : records(), nrMappedRecords_(0), nrUnmappedRecords_(0) {}

size_t SamBlock::nrMappedRecords() const { return nrMappedRecords_; }

size_t SamBlock::nrUnmappedRecords() const { return nrUnmappedRecords_; }

size_t SamBlock::nrRecords() const { return (nrMappedRecords_ + nrUnmappedRecords_); }

void SamBlock::reset() {
    records.clear();
    nrMappedRecords_ = 0;
    nrUnmappedRecords_ = 0;
}

}  // namespace cip

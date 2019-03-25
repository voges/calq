#include "calq/sam_block.h"

namespace calq {

SAMBlock::SAMBlock(void)
    : records(),
      nrMappedRecords_(0),
      nrUnmappedRecords_(0) {}

SAMBlock::~SAMBlock(void) {}

size_t SAMBlock::nrMappedRecords(void) const {
    return nrMappedRecords_;
}

size_t SAMBlock::nrUnmappedRecords(void) const {
    return nrUnmappedRecords_;
}

size_t SAMBlock::nrRecords(void) const {
    return (nrMappedRecords_ + nrUnmappedRecords_);
}

void SAMBlock::reset(void) {
    records.clear();
    nrMappedRecords_ = 0;
    nrUnmappedRecords_ = 0;
}

}  // namespace calq

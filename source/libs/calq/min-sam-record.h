/**
 * @file min-sam-record.h
 */

#ifndef CALQ_MIN_SAM_RECORD_H_
#define CALQ_MIN_SAM_RECORD_H_

namespace calq {

struct MinSamRecord {
    uint32_t posMin;
    uint32_t posMax;
    std::string ref;  // = RNAME
    std::string cigar;
    std::string seq;
    std::string qual;
};

}  // namespace calq

#endif  // CALQ_MIN_SAM_RECORD_H_

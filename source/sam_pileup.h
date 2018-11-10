/** @file SAMPileup.h
 *  @brief This file contains the definition of the SAMPileup class.
 */

// Copyright 2015-2017 Leibniz Universitaet Hannover

#ifndef CALQ_IO_SAM_SAMPILEUP_H_
#define CALQ_IO_SAM_SAMPILEUP_H_

#include <string>
#include <vector>

#include "constants.h"

namespace calq {

class SAMPileup {
 public:
    SAMPileup();
    ~SAMPileup();

    bool empty() const;
    void clear();
    void print() const;
    void printQual() const;
    void printSeq() const;

    uint32_t pos;  // 0-based position of this pileup
    std::string qual;
    std::string seq;

    uint16_t hq_softcounter;  // High quality softclips next to this position
};

}  // namespace calq

#endif  // CALQ_IO_SAM_SAMPILEUP_H_

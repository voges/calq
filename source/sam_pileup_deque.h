/** @file SAMPileupDeque.h
 *  @brief This file contains the definition of the SAMPileupDeque class.
 */

// Copyright 2015-2017 Leibniz Universitaet Hannover

#ifndef CALQ_IO_SAM_SAMPILEUPDEQUE_H_
#define CALQ_IO_SAM_SAMPILEUPDEQUE_H_

#include <deque>

#include "sam_pileup.h"

namespace calq {

class SAMPileupDeque {
    friend class SAMRecord;

 public:
    SAMPileupDeque();
    ~SAMPileupDeque();

    const SAMPileup &back() const;
    void clear();
    bool empty() const;
    const SAMPileup &front() const;
    size_t length() const;
    const SAMPileup &operator[](const size_t &n) const;
    void pop_back();
    void pop_front();
    size_t size() const;

    void print() const;

    uint32_t posMax() const;
    uint32_t posMin() const;

    void setPosMax(const uint32_t &posMax);
    void setPosMin(const uint32_t &posMin);

 private:
    std::deque<SAMPileup> pileups_;
    uint32_t posMax_;
    uint32_t posMin_;
};

}  // namespace calq

#endif  // CALQ_IO_SAM_SAMPILEUPDEQUE_H_

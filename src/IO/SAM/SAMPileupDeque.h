/** @file SAMPileupDeque.h
 *  @brief This file contains the definition of the SAMPileupDeque class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#ifndef CALQ_IO_SAM_SAMPILEUPDEQUE_H_
#define CALQ_IO_SAM_SAMPILEUPDEQUE_H_

#include "IO/SAM/SAMPileup.h"

#include <deque>

namespace calq {

class SAMPileupDeque {
    friend class SAMRecord;

public:
    explicit SAMPileupDeque(void);
    ~SAMPileupDeque(void);

    bool empty(void) const;
    void clear(void);
    const SAMPileup & back(void) const;
    const SAMPileup & front(void) const;
    size_t length(void) const;
    uint32_t posMax(void) const;
    uint32_t posMin(void) const;
    void pop_back(void);
    void pop_front(void);
    void print(void) const;
    void setPosMax(const uint32_t &posMax);
    void setPosMin(const uint32_t &posMin);

private:
    std::deque<SAMPileup> pileups_;
    uint32_t posMax_;
    uint32_t posMin_;
};

} // namespace calq

#endif // CALQ_IO_SAM_SAMPILEUPDEQUE_H_


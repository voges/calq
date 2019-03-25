#ifndef CALQ_SAM_PILEUP_DEQUE_H_
#define CALQ_SAM_PILEUP_DEQUE_H_

#include <deque>

#include "calq/sam_pileup.h"

namespace calq {

class SAMPileupDeque {
    friend class SAMRecord;

 public:
    SAMPileupDeque(void);
    ~SAMPileupDeque(void);

    const SAMPileup & back(void) const;
    void clear(void);
    bool empty(void) const;
    const SAMPileup & front(void) const;
    size_t length(void) const;
    const SAMPileup & operator[](const size_t &n) const;
    void pop_back(void);
    void pop_front(void);
    size_t size(void) const;

    void print(void) const;

    uint32_t posMax(void) const;
    uint32_t posMin(void) const;

    void setPosMax(const uint32_t &posMax);
    void setPosMin(const uint32_t &posMin);

 private:
    std::deque<SAMPileup> pileups_;
    uint32_t posMax_;
    uint32_t posMin_;
};

}  // namespace calq

#endif  // CALQ_SAM_PILEUP_DEQUE_H_

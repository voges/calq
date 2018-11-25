#ifndef CALQ_SAM_PILEUP_H_
#define CALQ_SAM_PILEUP_H_

#include <string>

namespace calq {

class SAMPileup {
 public:
    SAMPileup(void);
    ~SAMPileup(void);

    bool empty(void) const;
    void clear(void);
    void print(void) const;
    void printQual(void) const;
    void printSeq(void) const;

    uint32_t pos;  // 0-based position of this pileup
    std::string qual;
    std::string seq;
};

}  // namespace calq

#endif  // CALQ_SAM_PILEUP_H_

/**
 * @file sam-pileup.h
 */

#ifndef CALQ_SAM_PILEUP_H_
#define CALQ_SAM_PILEUP_H_

#include <deque>
#include <string>
#include <vector>
#include "encoding-read.h"

namespace calq {

class SAMPileup {
   public:
    SAMPileup();

    uint32_t pos;            // 0-based position of this pileup
    std::string qual;        // Quality values
    std::string seq;         // Nucleotides
    char ref;                // Reference base
    uint16_t hqSoftclipCnt;  // High quality softclips next to this position
};

class SamPileupDeque {
   public:
    SamPileupDeque();
    const SAMPileup& back() const;
    void clear();
    bool empty() const;
    const SAMPileup& front() const;
    size_t length() const;
    const SAMPileup& operator[](const size_t& n) const;
    void pop_front();
    size_t size() const;
    uint32_t posMax() const;
    uint32_t posMin() const;
    void setPosMax(const uint32_t& posMax);
    void setPosMin(const uint32_t& posMin);
    void add(const EncodingRead& r, uint8_t qvOffset, uint8_t hqSoftClipThreshold);

   private:
    std::deque<SAMPileup> pileups_{};
    uint32_t posMax_;
    uint32_t posMin_;
};

}  // namespace calq

#endif  // CALQ_SAM_PILEUP_H_

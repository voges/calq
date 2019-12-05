#ifndef CALQ_SAM_PILEUP_H_
#define CALQ_SAM_PILEUP_H_

#include <deque>
#include <string>
#include <vector>
#include "min-sam-record.h"

namespace calq {

class SamPileup {
   public:
    SamPileup();
    SamPileup(const SamPileup &) = delete;
    SamPileup &operator=(const SamPileup &) = delete;
    SamPileup(SamPileup &&) = delete;
    SamPileup &operator=(SamPileup &&) = delete;
    ~SamPileup() = default;

    uint32_t pos;            // 0-based position of this pileup
    std::string qual;        // Quality values
    std::string seq;         // Nucleotides
    char ref;                // Reference base
    uint16_t hqSoftclipCnt;  // High quality softclips next to this position
};

class SamPileupDeque {
   public:
    SamPileupDeque();
    SamPileupDeque(const SamPileupDeque &) = delete;
    SamPileupDeque &operator=(const SamPileupDeque &) = delete;
    SamPileupDeque(SamPileupDeque &&) = delete;
    SamPileupDeque &operator=(SamPileupDeque &&) = delete;
    ~SamPileupDeque() = default;

    const SamPileup &back() const;
    void clear();
    bool empty() const;
    const SamPileup &front() const;
    size_t length() const;
    const SamPileup &operator[](const size_t &n) const;
    void pop_front();
    size_t size() const;
    uint32_t posMax() const;
    uint32_t posMin() const;
    void setPosMax(uint32_t posMax);
    void setPosMin(uint32_t posMin);
    void add(const MinSamRecord &r, uint8_t qualOffset, uint8_t hqSoftClipThreshold);

   private:
    std::deque<SamPileup> pileups_{};
    uint32_t posMax_;
    uint32_t posMin_;
};

}  // namespace calq

#endif  // CALQ_SAM_PILEUP_H_

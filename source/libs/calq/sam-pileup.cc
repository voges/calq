/**
 * @file sam-pileup.cc
 */

#include "sam-pileup.h"
#include "errors.h"
#include "min-sam-record.h"

namespace calq {

SamPileup::SamPileup() : pos(0), qual(""), seq(""), ref('N'), hqSoftclipCnt(0) {}

SamPileupDeque::SamPileupDeque() : pileups_(), posMax_(0), posMin_(0) {}

const SamPileup& SamPileupDeque::back() const {
    if (pileups_.empty()) {
        throwErrorException("Cannot access back() of empty deque");
    }
    return pileups_.back();
}

void SamPileupDeque::clear() {
    pileups_.clear();
    posMax_ = 0;
    posMin_ = 0;
}

bool SamPileupDeque::empty() const { return pileups_.empty(); }

const SamPileup& SamPileupDeque::front() const {
    if (pileups_.empty()) {
        throwErrorException("Cannot access front() of empty deque");
    }
    return pileups_.front();
}

size_t SamPileupDeque::length() const { return posMax_ - posMin_ + 1; }

const SamPileup& SamPileupDeque::operator[](const size_t& n) const { return pileups_.at(n); }

void SamPileupDeque::pop_front() {
    if (pileups_.empty()) {
        throwErrorException("Deque is empty");
    }
    pileups_.pop_front();
    posMin_++;
}

size_t SamPileupDeque::size() const { return pileups_.size(); }

uint32_t SamPileupDeque::posMax() const { return posMax_; }

uint32_t SamPileupDeque::posMin() const { return posMin_; }

void SamPileupDeque::setPosMax(const uint32_t posMax) {
    if (posMax < posMax_) {
        throwErrorException("posMax out of range");
    }
    posMax_ = posMax;
    pileups_.resize(length());
}

void SamPileupDeque::setPosMin(const uint32_t posMin) {
    if (posMin < posMin_) {
        throwErrorException("posMin out of range");
    }

    if (empty()) {
        posMin_ = posMin;
    } else {
        for (uint32_t i = posMin_; i < posMin; i++) {
            pop_front();
        }
    }
}

void SamPileupDeque::add(const MinSamRecord& r, const uint8_t qualOffset, const uint8_t hqSoftClipThreshold) {
    if (this->empty()) {
        throwErrorException("Pileup queue is empty");
    }
    if ((this->posMin() > r.posMin) || (this->posMax() < r.posMax - 1)) {
        throwErrorException("samPileupQueue does not overlap record");
    }

    size_t cigarIdx = 0;
    size_t cigarLen = r.cigar.length();
    size_t opLen = 0;  // length of current CIGAR operation
    size_t idx = 0;
    size_t pileupIdx = r.posMin - this->posMin();

    size_t softclips = 0;

    for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++) {
        if (isdigit(r.cigar[cigarIdx])) {
            opLen = opLen * 10 + (size_t)r.cigar[cigarIdx] - (size_t)'0';
            continue;
        }

        const auto HQ_SOFTCLIP_THRESHOLD = static_cast<const char>(hqSoftClipThreshold + qualOffset);

        switch (r.cigar[cigarIdx]) {
            case 'M':
            case '=':
            case 'X':
                for (size_t i = 0; i < opLen; i++) {
                    this->pileups_[pileupIdx].pos = static_cast<uint32_t>(this->posMin() + pileupIdx);
                    this->pileups_[pileupIdx].seq += r.seq[idx];
                    this->pileups_[pileupIdx].qual += r.qual[idx];
                    if (!r.ref.empty()) {
                        if (this->pileups_[pileupIdx].ref != 'N' &&
                            this->pileups_[pileupIdx].ref != r.ref[pileupIdx + this->posMin() - r.posMin]) {
                            throwErrorException("Non-matching reference between reads!");
                        }
                        this->pileups_[pileupIdx].ref = r.ref[pileupIdx + this->posMin() - r.posMin];
                    }

                    idx++;
                    pileupIdx++;
                }
                break;

            case 'S':
                for (int l = 0; l < static_cast<int>(opLen); ++l) {
                    if (r.qual[idx + l] >= HQ_SOFTCLIP_THRESHOLD) {
                        ++softclips;
                    }
                }
                // Fall-through
            case 'I':
                idx += opLen;
                break;
            case 'D':
            case 'N':
                pileupIdx += opLen;
                break;
            case 'H':
            case 'P':
                break;  // These have been clipped
            default:
                throwErrorException("Bad CIGAR string");
        }

        opLen = 0;
    }

    // Write clips
    for (size_t i = r.posMin - this->posMin(); i < pileupIdx; ++i) {
        this->pileups_[i].hqSoftclipCnt += softclips;
    }
}

}  // namespace calq

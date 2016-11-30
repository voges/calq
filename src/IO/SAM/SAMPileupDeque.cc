/** @file SAMPileupDeque.cc
 *  @brief This file contains the implementation of the SAMPileupDeque class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#include "IO/SAM/SAMPileupDeque.h"

#include "Common/Exceptions.h"

namespace calq {

SAMPileupDeque::SAMPileupDeque(void) : pileups_() , posMax_(0) , posMin_(0) {}

SAMPileupDeque::~SAMPileupDeque(void) {}

bool SAMPileupDeque::empty(void) const
{
    if (pileups_.empty() == true)
        return true;
    return false;
}

void SAMPileupDeque::clear(void)
{
    pileups_.clear();
    posMax_ = 0;
    posMin_ = 0;
}

const SAMPileup& SAMPileupDeque::back(void) const
{
    if (pileups_.empty() == true) {
        throwErrorException("Deque is empty");
    }
    return pileups_.back();
}


const SAMPileup & SAMPileupDeque::front(void) const
{
    if (pileups_.empty() == true) {
        throwErrorException("Deque is empty");
    }
    return pileups_.front();
}

size_t SAMPileupDeque::length(void) const
{
    return posMax_ - posMin_ + 1;
}

uint32_t SAMPileupDeque::posMax(void) const
{
    return posMax_;
}

uint32_t SAMPileupDeque::posMin(void) const
{
    return posMin_;
}

void SAMPileupDeque::pop_back(void)
{
    if (pileups_.empty() == true) {
        throwErrorException("Deque is empty");
    }
    pileups_.pop_back();
    posMax_--;
}

void SAMPileupDeque::pop_front(void)
{
    if (pileups_.empty() == true) {
        throwErrorException("Deque is empty");
    }
    pileups_.pop_front();
    posMin_++;
}

void SAMPileupDeque::print(void) const
{
    if (pileups_.empty() == true) {
        throwErrorException("Deque is empty");
    }

    for (auto const &samPileup : pileups_) {
        samPileup.printSeq();
    }
}

void SAMPileupDeque::setPosMax(const uint32_t &posMax)
{
    if (posMax < posMax_) {
        throwErrorException("posMax range");
    }
    posMax_ = posMax;
    pileups_.resize(length());
}

void SAMPileupDeque::setPosMin(const uint32_t &posMin)
{
    if (posMin < posMin_) {
        throwErrorException("posMin range");
    }

    if (empty() == true) {
        posMin_ = posMin;
    } else {
        for (uint32_t i = posMin_; i < posMin ; i++) { pop_front(); } 
    }
}

} // namespace calq


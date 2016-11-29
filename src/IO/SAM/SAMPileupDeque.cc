/** @file SAMPileupDeque.cc
 *  @brief This file contains the implementation of the SAMPileupDeque class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "IO/SAM/SAMPileupDeque.h"
#include "Common/Exceptions.h"

cq::SAMPileupDeque::SAMPileupDeque(void)
    : pileups_()
    , posMax_(0)
    , posMin_(0)
{
    // empty
}

cq::SAMPileupDeque::~SAMPileupDeque(void)
{
    // empty
}

bool cq::SAMPileupDeque::empty(void) const
{
    if (pileups_.empty() == true)
        return true;
    return false;
}

void cq::SAMPileupDeque::clear(void)
{
    pileups_.clear();
    posMax_ = 0;
    posMin_ = 0;
}

const cq::SAMPileup& cq::SAMPileupDeque::back(void) const
{
    if (pileups_.empty() == true) {
        throwErrorException("Deque is empty");
    }
    return pileups_.back();
}


const cq::SAMPileup & cq::SAMPileupDeque::front(void) const
{
    if (pileups_.empty() == true) {
        throwErrorException("Deque is empty");
    }
    return pileups_.front();
}

size_t cq::SAMPileupDeque::length(void) const
{
    return posMax_ - posMin_ + 1;
}

uint32_t cq::SAMPileupDeque::posMax(void) const
{
    return posMax_;
}

uint32_t cq::SAMPileupDeque::posMin(void) const
{
    return posMin_;
}

void cq::SAMPileupDeque::pop_back(void)
{
    if (pileups_.empty() == true) {
        throwErrorException("Deque is empty");
    }
    pileups_.pop_back();
    posMax_--;
}

void cq::SAMPileupDeque::pop_front(void)
{
    if (pileups_.empty() == true) {
        throwErrorException("Deque is empty");
    }
    pileups_.pop_front();
    posMin_++;
}

void cq::SAMPileupDeque::print(void) const
{
    if (pileups_.empty() == true) {
        throwErrorException("Deque is empty");
    }

    for (auto const &samPileup : pileups_) {
        samPileup.printSeq();
    }
}

void cq::SAMPileupDeque::setPosMax(const uint32_t &posMax)
{
    if (posMax < posMax_) {
        throwErrorException("posMax range");
    }
    posMax_ = posMax;
    pileups_.resize(length());
}

void cq::SAMPileupDeque::setPosMin(const uint32_t &posMin)
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


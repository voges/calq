/** @file Pileup.cc
 *  @brief This file contains the implementations of the Pileup and PileupQueue
 *         classes.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "QualCodec/Pileup.h"
#include "Common/Exceptions.h"

cq::Pileup::Pileup(void)
    : seq("")
    , qual("")
{
    // empty
}

cq::Pileup::~Pileup(void)
{
    // empty
}

void cq::Pileup::print(void) const
{
    printf("seq: %s\n", seq.c_str());
    printf("qual: %s\n", qual.c_str());
}

cq::PileupQueue::PileupQueue(void)
    : pileups()
    , posMax_(0)
    , posMin_(0)
{
    // empty
}

cq::PileupQueue::~PileupQueue(void)
{
    // empty
}

void cq::PileupQueue::clear(void)
{
    pileups.clear();
    posMax_ = 0;
    posMin_ = 0;
}

size_t cq::PileupQueue::length(void) const
{
    return posMax_ - posMin_ + 1;
}

uint32_t cq::PileupQueue::posMax(void) const
{
    return posMax_;
}

uint32_t cq::PileupQueue::posMin(void) const
{
    return posMin_;
}

void cq::PileupQueue::pop_front(void)
{
    pileups.pop_front();
}

void cq::PileupQueue::setPosMax(const uint32_t &posMax)
{
    if (posMax < posMin_) {
        throwErrorException("posMax is smaller than posMin");
    }
}

void cq::PileupQueue::setPosMin(const uint32_t &posMin)
{

}


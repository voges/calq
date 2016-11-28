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

#include "IO/SAM/Pileup.h"
#include "Common/Exceptions.h"

cq::Pileup::Pileup(void)
    : pos(0)
    , qual("")
    , seq("")
{
    // empty
}

cq::Pileup::~Pileup(void)
{
    // empty
}

bool cq::Pileup::empty(void) const
{
    if (seq.empty() == true)
        return true;
    return false;
}

void cq::Pileup::clear(void)
{
    pos = 0;
    qual = "";
    seq = "";
}

void cq::Pileup::print(void) const
{
    printQual();
    printSeq();
}

void cq::Pileup::printQual(void) const
{
    printf("qual: %s\n", qual.c_str());
}

void cq::Pileup::printSeq(void) const
{
    printf("seq: %s\n", seq.c_str());
}

cq::PileupQueue::PileupQueue(void)
    : pileups_()
    , posMax_(0)
    , posMin_(0)
{
    // empty
}

cq::PileupQueue::~PileupQueue(void)
{
    // empty
}

bool cq::PileupQueue::empty(void) const
{
    if (pileups_.empty() == true)
        return true;
    return false;
}

void cq::PileupQueue::clear(void)
{
    pileups_.clear();
    posMax_ = 0;
    posMin_ = 0;
}

const cq::Pileup& cq::PileupQueue::front(void) const
{
    return pileups_.front();
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
    pileups_.pop_front();
}

void cq::PileupQueue::setPosMax(const uint32_t &posMax)
{
    if (posMax <= posMax_) {
        throwErrorException("posMax range");
    }

    posMax_ = posMax;
    pileups_.resize(length());
}

void cq::PileupQueue::setPosMin(const uint32_t &posMin)
{
    if (posMin <= posMin_) {
        throwErrorException("posMin range");
    }

    for (uint32_t i = posMin_; i < posMin ; i++) { pileups_.pop_front(); } 
    posMin_ = posMin;
}


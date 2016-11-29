/** @file SAMPileup.cc
 *  @brief This file contains the implementation of the SAMPileup class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "IO/SAM/SAMPileup.h"

cq::SAMPileup::SAMPileup(void)
    : pos(0)
    , qual("")
    , seq("")
{
    // empty
}

cq::SAMPileup::~SAMPileup(void)
{
    // empty
}

bool cq::SAMPileup::empty(void) const
{
    if (seq.empty() == true)
        return true;
    return false;
}

void cq::SAMPileup::clear(void)
{
    pos = 0;
    qual = "";
    seq = "";
}

void cq::SAMPileup::print(void) const
{
    printQual();
    printSeq();
}

void cq::SAMPileup::printQual(void) const
{
    printf("%6d: %s\n", pos, qual.c_str());
}

void cq::SAMPileup::printSeq(void) const
{
    printf("%6d: %s\n", pos, seq.c_str());
}


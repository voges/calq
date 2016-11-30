/** @file SAMPileup.cc
 *  @brief This file contains the implementation of the SAMPileup class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#include "IO/SAM/SAMPileup.h"

namespace calq {

SAMPileup::SAMPileup(void) : pos(0), qual(""), seq("") {}

SAMPileup::~SAMPileup(void) {}

bool SAMPileup::empty(void) const
{
    if (seq.empty() == true)
        return true;
    return false;
}

void SAMPileup::clear(void)
{
    pos = 0;
    qual = "";
    seq = "";
}

void SAMPileup::print(void) const
{
    printQual();
    printSeq();
}

void SAMPileup::printQual(void) const
{
    printf("%6d: %s\n", pos, qual.c_str());
}

void SAMPileup::printSeq(void) const
{
    printf("%6d: %s\n", pos, seq.c_str());
}

} // namespace calq


#include "sam-pileup.h"

namespace calq {

SAMPileup::SAMPileup() : pos(0), qual(""), seq(""), ref('N'), hq_softcounter(0) {}

bool SAMPileup::empty() const { return seq.empty(); }

void SAMPileup::clear() {
    pos = 0;
    qual = "";
    seq = "";
    ref = 'N';
}

void SAMPileup::print() const { printf("%6u: %s %s\n", pos, seq.c_str(), qual.c_str()); }

void SAMPileup::printQual() const { printf("%6u: %s\n", pos, qual.c_str()); }

void SAMPileup::printSeq() const { printf("%6u: %s\n", pos, seq.c_str()); }

}  // namespace calq

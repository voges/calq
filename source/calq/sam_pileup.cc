#include "calq/sam_pileup.h"

namespace calq {
SAMPileup::SAMPileup() : pos(0), qual(""), seq(""), ref('N'), hq_softcounter(0) {}

SAMPileup::~SAMPileup() = default;

bool SAMPileup::empty() const {
    return seq.empty();
}

void SAMPileup::clear() {
    pos = 0;
    qual = "";
    seq = "";
    ref = 'N';
}

void SAMPileup::print() const {
    printf("%6d: %s %s\n", pos, seq.c_str(), qual.c_str());
}

void SAMPileup::printQual() const {
    printf("%6d: %s\n", pos, qual.c_str());
}

void SAMPileup::printSeq() const {
    printf("%6d: %s\n", pos, seq.c_str());
}

}  // namespace calq

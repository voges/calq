#include "calq/sam-pileup.h"
#include <gtest/gtest.h>

TEST(SamPileupDeque, Generic) {  // NOLINT(cert-err58-cpp)
    calq::MinSamRecord record1;
    record1.posMin = 0;
    record1.posMax = 6;
    record1.ref = "chr1";
    record1.cigar = "7M";
    record1.seq = "GATTACA";
    record1.qual = "QQQQQQQ";

    calq::SamPileupDeque deque;

    // TODO(Jan): the function add() should work w/o requiring the user to call setPos[Min|Max]() first
    deque.setPosMin(record1.posMin);
    deque.setPosMax(record1.posMax);

    uint8_t qualOffset = 33;
    uint8_t hqSoftclipThreshold = 0;

    deque.add(record1, qualOffset, hqSoftclipThreshold);
}

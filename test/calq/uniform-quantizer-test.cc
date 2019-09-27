#include "calq/uniform-quantizer.h"
#include <gtest/gtest.h>

TEST(UniformQuantizer, Initialization) {  // NOLINT(cert-err58-cpp)
    calq::UniformQuantizer q3(0, 1, 2);
}

TEST(UniformQuantizer, PositiveRange) {  // NOLINT(cert-err58-cpp)
    calq::UniformQuantizer q(0, 4, 2);

    EXPECT_EQ(q.valueToIndex(0), 0);
    EXPECT_EQ(q.valueToIndex(1), 0);
    EXPECT_EQ(q.valueToIndex(2), 0);
    EXPECT_EQ(q.valueToIndex(3), 1);
    EXPECT_EQ(q.valueToIndex(4), 1);

    EXPECT_EQ(q.indexToReconstructedValue(0), 1);
    EXPECT_EQ(q.indexToReconstructedValue(1), 3);
}

TEST(UniformQuantizer, NegativeRange) {  // NOLINT(cert-err58-cpp)
    calq::UniformQuantizer q(-5, -1, 2);

    EXPECT_EQ(q.valueToIndex(-5), 0);
    EXPECT_EQ(q.valueToIndex(-4), 0);
    EXPECT_EQ(q.valueToIndex(-3), 0);
    EXPECT_EQ(q.valueToIndex(-2), 1);
    EXPECT_EQ(q.valueToIndex(-1), 1);

    EXPECT_EQ(q.indexToReconstructedValue(0), -4);
    EXPECT_EQ(q.indexToReconstructedValue(1), -2);
}

TEST(UniformQuantizer, EightBinning) {  // NOLINT(cert-err58-cpp)
    calq::UniformQuantizer q(0, 40, 8);

    EXPECT_EQ(q.valueToIndex(0), 0);
    EXPECT_EQ(q.valueToIndex(1), 0);
    EXPECT_EQ(q.valueToIndex(2), 0);
    EXPECT_EQ(q.valueToIndex(3), 0);
    EXPECT_EQ(q.valueToIndex(4), 0);
    EXPECT_EQ(q.valueToIndex(5), 0);
    EXPECT_EQ(q.valueToIndex(6), 1);
    EXPECT_EQ(q.valueToIndex(7), 1);
    EXPECT_EQ(q.valueToIndex(8), 1);
    EXPECT_EQ(q.valueToIndex(9), 1);
    EXPECT_EQ(q.valueToIndex(10), 1);
    EXPECT_EQ(q.valueToIndex(11), 2);
    EXPECT_EQ(q.valueToIndex(12), 2);
    EXPECT_EQ(q.valueToIndex(13), 2);
    EXPECT_EQ(q.valueToIndex(14), 2);
    EXPECT_EQ(q.valueToIndex(15), 2);
    EXPECT_EQ(q.valueToIndex(16), 3);
    EXPECT_EQ(q.valueToIndex(17), 3);
    EXPECT_EQ(q.valueToIndex(18), 3);
    EXPECT_EQ(q.valueToIndex(19), 3);
    EXPECT_EQ(q.valueToIndex(20), 3);
    EXPECT_EQ(q.valueToIndex(21), 4);
    EXPECT_EQ(q.valueToIndex(22), 4);
    EXPECT_EQ(q.valueToIndex(23), 4);
    EXPECT_EQ(q.valueToIndex(24), 4);
    EXPECT_EQ(q.valueToIndex(25), 4);
    EXPECT_EQ(q.valueToIndex(26), 5);
    EXPECT_EQ(q.valueToIndex(27), 5);
    EXPECT_EQ(q.valueToIndex(28), 5);
    EXPECT_EQ(q.valueToIndex(29), 5);
    EXPECT_EQ(q.valueToIndex(30), 5);
    EXPECT_EQ(q.valueToIndex(31), 6);
    EXPECT_EQ(q.valueToIndex(32), 6);
    EXPECT_EQ(q.valueToIndex(33), 6);
    EXPECT_EQ(q.valueToIndex(34), 6);
    EXPECT_EQ(q.valueToIndex(35), 6);
    EXPECT_EQ(q.valueToIndex(36), 7);
    EXPECT_EQ(q.valueToIndex(37), 7);
    EXPECT_EQ(q.valueToIndex(38), 7);
    EXPECT_EQ(q.valueToIndex(39), 7);
    EXPECT_EQ(q.valueToIndex(40), 7);

    EXPECT_EQ(q.indexToReconstructedValue(0), 3);
    EXPECT_EQ(q.indexToReconstructedValue(1), 8);
    EXPECT_EQ(q.indexToReconstructedValue(2), 13);
    EXPECT_EQ(q.indexToReconstructedValue(3), 18);
    EXPECT_EQ(q.indexToReconstructedValue(4), 23);
    EXPECT_EQ(q.indexToReconstructedValue(5), 28);
    EXPECT_EQ(q.indexToReconstructedValue(6), 33);
    EXPECT_EQ(q.indexToReconstructedValue(7), 38);
}

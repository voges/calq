#include "calq/quantizer.h"
#include <gtest/gtest.h>
#include "calq/exceptions.h"

TEST(QuantizerTest, EmptyLut) {  // NOLINT(cert-err58-cpp)
    calq::Quantizer q;

    // The is no table in the quantizer. Thus we expect these functions to
    // throw exceptions.
    EXPECT_THROW(q.valueToIndex(0), calq::ErrorException);
    EXPECT_THROW(q.indexToReconstructedValue(0), calq::ErrorException);

    // The inverse LUT should be empty.
    const std::map<int, int>& inverseLut = q.inverseLut();
    EXPECT_EQ(inverseLut.size(), 0);
}

TEST(QuantizerTest, PopulatedLut) {  // NOLINT(cert-err58-cpp)
    std::map<int, int> inverseLut = {{0, 0}, {1, 0}, {2, 1}, {3, 1}};
    calq::Quantizer q(inverseLut);

    // The function valueToIndex() performs a look-up in the LUT; this is
    // however not initialized within the Quantizer. Thus we expect
    // an exception.
    EXPECT_THROW(q.valueToIndex(0), calq::ErrorException);
    EXPECT_THROW(q.valueToIndex(1), calq::ErrorException);
    EXPECT_THROW(q.valueToIndex(2), calq::ErrorException);
    EXPECT_THROW(q.valueToIndex(3), calq::ErrorException);

    // An index out of range should lead to an exception.
    EXPECT_THROW(q.indexToReconstructedValue(-1), calq::ErrorException);
    EXPECT_THROW(q.indexToReconstructedValue(4), calq::ErrorException);

    // These should work.
    EXPECT_EQ(q.indexToReconstructedValue(0), 0);
    EXPECT_EQ(q.indexToReconstructedValue(1), 0);
    EXPECT_EQ(q.indexToReconstructedValue(2), 1);
    EXPECT_EQ(q.indexToReconstructedValue(3), 1);
}

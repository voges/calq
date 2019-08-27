#include "calq/filter-buffer.h"
#include <gtest/gtest.h>

TEST(FilterBufferTest, Convolution) { // NOLINT(cert-err58-cpp)
    calq::FilterBuffer buffer([](size_t pos, size_t size) -> double { return pos / static_cast<double>(size - 1); }, 3);

    EXPECT_EQ(buffer.filter(), 0);
    buffer.push(1);
    EXPECT_EQ(buffer.filter(), 1);
    buffer.push(2);
    EXPECT_EQ(buffer.filter(), 2.5);
    buffer.push(3);
    EXPECT_EQ(buffer.filter(), 4);
    buffer.push(4);
}

TEST(FilterBufferTest, Offset) { // NOLINT(cert-err58-cpp)
    calq::FilterBuffer buffer([](size_t pos, size_t size) -> double { return pos / static_cast<double>(size - 1); }, 3);
    EXPECT_EQ(buffer.getOffset(), 2);
}

// void gaussKernelTest() {
//     // Test variance 1
//     GaussKernel k(1.0);
//
//     double buffer[63];
//     for (size_t i = 0; i < 63; ++i) {
//         buffer[i] = k.calcValue(i, 63);
//     }
//
//     equals(buffer[31], 0.3989);
//     equals(buffer[32], 0.2419);
//     equals(buffer[30], 0.2419);
//     equals(buffer[33], 0.0539);
//     equals(buffer[29], 0.0539);
//
//     // Test variance 17
//     GaussKernel k2(17.0);
//     double buffer2[127];
//     for (size_t i = 0; i < 127; ++i) {
//         buffer2[i] = k2.calcValue(i, 127);
//     }
//
//     equals(buffer2[63], 0.0234);
//     equals(buffer2[62], 0.0234);
//     equals(buffer2[64], 0.0234);
//     equals(buffer2[61], 0.0233);
//     equals(buffer2[65], 0.0233);
//     equals(buffer2[41], 0.0101);
//     equals(buffer2[85], 0.0101);
//
//     // Test length limits
//     equals(k2.calcMinSize(0.01), 47);
//     equals(k.calcMinSize(0.01), 7);
// }

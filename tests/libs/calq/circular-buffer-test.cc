#include "calq/circular-buffer.h"
#include <gtest/gtest.h>

TEST(CircularBufferTest, Size) {
    calq::CircularBuffer<int> buffer(3, 0);

    EXPECT_EQ(buffer.size(), 3);
}

TEST(CircularBufferTest, Access) {
    calq::CircularBuffer<int> buffer(3, 0);

    buffer.push(1);
    buffer.push(2);
    buffer.push(3);

    EXPECT_EQ(buffer.front(), 3);
    EXPECT_EQ(buffer[2], 3);
    EXPECT_EQ(buffer.back(), 1);
    EXPECT_EQ(buffer[0], 1);
    EXPECT_EQ(buffer[1], 2);
}

#include "calq/haplotyper.h"
#include <calq/calq-codec.h>
#include <gtest/gtest.h>

TEST(HaplotyperTest, Everything) {  // NOLINT(cert-err58-cpp)
    calq::Haplotyper h(5, 2, 33, 8, 5, 3, 5, false, true, calq::FilterType::GAUSS);

    EXPECT_EQ(h.getOffset(), 10);

    // First 5 pushes inside basespreader
    for (int i = 0; i < 5; ++i) {
        EXPECT_EQ(h.push("C", "}", 0, 'A'), 0);
    }

    // next 11 pushes inside filterbuffer approaching 0.73
    for (int i = 0; i < 10; ++i) {
        h.push("C", "}", 0, 'A');
    }
    EXPECT_EQ(h.push("C", "}", 0, 'A'), 7);
    EXPECT_EQ(h.push("C", "}", 0, 'A'), 7);

    // reset
    calq::Haplotyper h2(5, 2, 33, 8, 5, 3, 5, false, true, calq::FilterType::GAUSS);

    h2.push("CCC", "}}}", 15, 'A');

    // push spreaded bases into filterbuffer
    for (int i = 0; i < 8; ++i) {
        h2.push("CCCCCCCCCCCCCC", "}}}}}}}}}}}}", 0, 'C');  // Activity close to 0
    }
    // Reach maximum
    EXPECT_EQ(h2.push("CCCCCCCCCCCCCC", "}}}}}}}}}}}}", 0, 'C'), 7);
    EXPECT_EQ(h2.push("CCCCCCCCCCCCCC", "}}}}}}}}}}}}", 0, 'C'), 7);
    EXPECT_EQ(h2.push("CCCCCCCCCCCCCC", "}}}}}}}}}}}}", 0, 'C'), 7);

    // Reach minimum
    for (int i = 0; i < 9; ++i) {
        h2.push("CCCCCCCCCCCCCC", "}}}}}}}}}}}}", 0, 'C');
    }

    EXPECT_EQ(h2.push("CCCCCCCCCCCCCC", "}}}}}}}}}}}}", 0, 'C'), 0);
}

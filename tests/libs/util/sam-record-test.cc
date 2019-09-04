#include "util/sam-record.h"
#include <gtest/gtest.h>
#include "helpers.h"

TEST(SamRecord, Minimal) {  // NOLINT(cert-err58-cpp)
    std::string gitRootDir = exec("git rev-parse --show-toplevel");
}

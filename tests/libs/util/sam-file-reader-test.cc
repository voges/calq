#include "util/sam-file-reader.h"
#include <gtest/gtest.h>
#include "helpers.h"

TEST(SamFileReader, Minimal) {  // NOLINT(cert-err58-cpp)
    std::string gitRootDir = exec("git rev-parse --show-toplevel");
}

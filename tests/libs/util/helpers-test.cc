#include "helpers.h"
#include <gtest/gtest.h>
#include "util/helpers.h"

TEST(Helpers, Minimal) {  // NOLINT(cert-err58-cpp)
    std::string gitRootDir = exec("git rev-parse --show-toplevel");
}

#include "util/file-writer.h"
#include <gtest/gtest.h>
#include "helpers.h"

TEST(FileWriter, Minimal) {  // NOLINT(cert-err58-cpp)
    std::string gitRootDir = exec("git rev-parse --show-toplevel");
}

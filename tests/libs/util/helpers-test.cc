#include "helpers.h"
#include <gtest/gtest.h>
#include "util/helpers.h"

TEST(Helpers, FileExists) {  // NOLINT(cert-err58-cpp)
    std::string gitRootDir = exec("git rev-parse --show-toplevel");

    EXPECT_EQ(util::fileExists(gitRootDir + "/README.md"), true);
    EXPECT_EQ(util::fileExists(gitRootDir + "/nonexistent"), false);
}

TEST(Helpers, FileNameExtension) {  // NOLINT(cert-err58-cpp)
    std::string gitRootDir = exec("git rev-parse --show-toplevel");

    EXPECT_EQ(util::fileNameExtension(gitRootDir + "/README.md"), "md");
    EXPECT_EQ(util::fileNameExtension(gitRootDir + "/no-extension"), "");
}

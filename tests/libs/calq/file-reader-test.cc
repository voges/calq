#include "calq/file-reader.h"
#include <gtest/gtest.h>
#include "helpers.h"

TEST(FileReader, Constructor) {  // NOLINT(cert-err58-cpp)
    std::string gitRootDir = calq_tests::exec("git rev-parse --show-toplevel");
    calq::FileReader fileReader(gitRootDir + "/resources/test-files/tests/0x041f5aac");
}

TEST(FileReader, ReadLine) {  // NOLINT(cert-err58-cpp)
    std::string gitRootDir = calq_tests::exec("git rev-parse --show-toplevel");
    calq::FileReader fileReader(gitRootDir + "/resources/test-files/tests/test-text.txt");

    std::string line;

    fileReader.readLine(&line);
    EXPECT_EQ(line, "test-line");

    fileReader.readLine(&line);
    EXPECT_EQ(line, "  test-line-leading-ws");

    fileReader.readLine(&line);
    EXPECT_EQ(line, "test-line-trailing-ws");
}

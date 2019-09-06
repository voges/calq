#include <gtest/gtest.h>
#include "helpers.h"
#include "util/file-reader.h"

TEST(File, Constructor) {  // NOLINT(cert-err58-cpp)
    std::string gitRootDir = exec("git rev-parse --show-toplevel");
    util::FileReader reader(gitRootDir + "/resources/test-files/tests/0x041f5aac");
}

TEST(File, AdvanceAndTell) {  // NOLINT(cert-err58-cpp)
    std::string gitRootDir = exec("git rev-parse --show-toplevel");
    util::FileReader reader(gitRootDir + "/resources/test-files/tests/0x041f5aac");

    EXPECT_EQ(reader.tell(), 0);
    EXPECT_NO_THROW(reader.advance(0));
    EXPECT_EQ(reader.tell(), 0);
    EXPECT_NO_THROW(reader.advance(1));
    EXPECT_EQ(reader.tell(), 1);
    EXPECT_NO_THROW(reader.advance(1));
    EXPECT_EQ(reader.tell(), 2);
    EXPECT_NO_THROW(reader.advance(1));
    EXPECT_EQ(reader.tell(), 3);
}

TEST(File, Eof) {  // NOLINT(cert-err58-cpp)
    std::string gitRootDir = exec("git rev-parse --show-toplevel");
    util::FileReader reader(gitRootDir + "/resources/test-files/tests/0x041f5aac");

    EXPECT_EQ(reader.eof(), false);
}

TEST(File, SeekFromCurAndTell) {  // NOLINT(cert-err58-cpp)
    std::string gitRootDir = exec("git rev-parse --show-toplevel");
    util::FileReader reader(gitRootDir + "/resources/test-files/tests/0x041f5aac");

    EXPECT_NO_THROW(reader.seekFromCur(0));
    EXPECT_EQ(reader.tell(), 0);
}

TEST(File, SeekFromEndAndTell) {  // NOLINT(cert-err58-cpp)
    std::string gitRootDir = exec("git rev-parse --show-toplevel");
    util::FileReader reader(gitRootDir + "/resources/test-files/tests/0x041f5aac");

    EXPECT_NO_THROW(reader.seekFromEnd(0));
    EXPECT_EQ(reader.tell(), 4);
}

TEST(File, SeekFromSetAndTell) {  // NOLINT(cert-err58-cpp)
    std::string gitRootDir = exec("git rev-parse --show-toplevel");
    util::FileReader reader(gitRootDir + "/resources/test-files/tests/0x041f5aac");

    EXPECT_NO_THROW(reader.seekFromSet(0));
    EXPECT_EQ(reader.tell(), 0);
}

TEST(File, Size) {  // NOLINT(cert-err58-cpp)
    std::string gitRootDir = exec("git rev-parse --show-toplevel");
    util::FileReader reader(gitRootDir + "/resources/test-files/tests/0x041f5aac");

    EXPECT_EQ(reader.size(), 4);
}

TEST(File, Tell) {  // NOLINT(cert-err58-cpp)
    std::string gitRootDir = exec("git rev-parse --show-toplevel");
    util::FileReader reader(gitRootDir + "/resources/test-files/tests/0x041f5aac");

    EXPECT_EQ(reader.tell(), 0);
    reader.advance(1);
    EXPECT_EQ(reader.tell(), 1);
}

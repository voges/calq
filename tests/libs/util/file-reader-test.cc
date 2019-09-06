#include "util/file-reader.h"
#include <gtest/gtest.h>
#include "helpers.h"

TEST(FileReader, Constructor) {  // NOLINT(cert-err58-cpp)
    std::string gitRootDir = exec("git rev-parse --show-toplevel");
    util::FileReader reader(gitRootDir + "/resources/test-files/tests/0x041f5aac");
}

TEST(FileReader, Advance) {  // NOLINT(cert-err58-cpp)
    std::string gitRootDir = exec("git rev-parse --show-toplevel");
    util::FileReader reader(gitRootDir + "/resources/test-files/tests/0x041f5aac");

    EXPECT_NO_THROW(reader.advance(0));
    EXPECT_NO_THROW(reader.advance(1));
    EXPECT_NO_THROW(reader.advance(2));
    EXPECT_NO_THROW(reader.advance(3));
}

TEST(FileReader, Eof) {  // NOLINT(cert-err58-cpp)
    std::string gitRootDir = exec("git rev-parse --show-toplevel");
    util::FileReader reader(gitRootDir + "/resources/test-files/tests/0x041f5aac");

    EXPECT_EQ(reader.eof(), false);
}

TEST(FileReader, Handle) {  // NOLINT(cert-err58-cpp)
    std::string gitRootDir = exec("git rev-parse --show-toplevel");
    util::FileReader reader(gitRootDir + "/resources/test-files/tests/0x041f5aac");

    EXPECT_NO_THROW(reader.handle());
}

TEST(FileReader, ReadLine) {  // NOLINT(cert-err58-cpp)
    std::string gitRootDir = exec("git rev-parse --show-toplevel");
    util::FileReader reader(gitRootDir + "/resources/test-files/tests/test-text.txt");

    std::string line;
    reader.readLine(&line);
    EXPECT_EQ(line, "test-text");
}

TEST(FileReader, SeekFromCur) {  // NOLINT(cert-err58-cpp)
    std::string gitRootDir = exec("git rev-parse --show-toplevel");
    util::FileReader reader(gitRootDir + "/resources/test-files/tests/0x041f5aac");

    EXPECT_NO_THROW(reader.seekFromCur(0));
}

TEST(FileReader, SeekFromEnd) {  // NOLINT(cert-err58-cpp)
    std::string gitRootDir = exec("git rev-parse --show-toplevel");
    util::FileReader reader(gitRootDir + "/resources/test-files/tests/0x041f5aac");

    EXPECT_NO_THROW(reader.seekFromEnd(0));
}

TEST(FileReader, SeekFromSet) {  // NOLINT(cert-err58-cpp)
    std::string gitRootDir = exec("git rev-parse --show-toplevel");
    util::FileReader reader(gitRootDir + "/resources/test-files/tests/0x041f5aac");

    EXPECT_NO_THROW(reader.seekFromSet(0));
}

TEST(FileReader, Size) {  // NOLINT(cert-err58-cpp)
    std::string gitRootDir = exec("git rev-parse --show-toplevel");
    util::FileReader reader(gitRootDir + "/resources/test-files/tests/0x041f5aac");

    EXPECT_EQ(reader.size(), 4);
}

TEST(FileReader, Tell) {  // NOLINT(cert-err58-cpp)
    std::string gitRootDir = exec("git rev-parse --show-toplevel");
    util::FileReader reader(gitRootDir + "/resources/test-files/tests/0x041f5aac");

    EXPECT_EQ(reader.tell(), 0);
    reader.advance(1);
    EXPECT_EQ(reader.tell(), 1);
}

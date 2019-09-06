#include "util/fasta-file-reader.h"
#include <gtest/gtest.h>
#include "helpers.h"
#include "util/fasta-record.h"

TEST(FastaFileReader, Minimal) {  // NOLINT(cert-err58-cpp)
    std::string gitRootDir = exec("git rev-parse --show-toplevel");

    util::FastaFileReader reader(gitRootDir + "/resources/test-files/fasta/minimal.fasta");

    std::vector<util::FastaRecord> records;
    reader.parse(&records);

    EXPECT_EQ(records.size(), 1);
    EXPECT_EQ(records.front().header, ">minimal");
    EXPECT_EQ(records.front().sequence, "GATTACA");
}

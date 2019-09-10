#include "calq/fasta-file-reader.h"
#include <gtest/gtest.h>
#include "calq/fasta-record.h"
#include "helpers.h"

TEST(FastaFileReader, Minimal) {  // NOLINT(cert-err58-cpp)
    std::string gitRootDir = calq_tests::exec("git rev-parse --show-toplevel");

    calq::FastaFileReader reader(gitRootDir + "/resources/test-files/fasta/minimal.fasta");

    std::vector<calq::FastaRecord> records;
    reader.parse(&records);

    EXPECT_EQ(records.size(), 1);
    EXPECT_EQ(records.front().header, ">minimal");
    EXPECT_EQ(records.front().sequence, "GATTACA");
}

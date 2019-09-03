#include "util/fasta-file-reader.h"
#include <gtest/gtest.h>
#include "util/fasta-record.h"
#include "util/string-helpers.h"

static std::string exec(const std::string& cmd) {
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) { return "ERROR"; }
    const int bufferSize = 256;
    char buffer[bufferSize];
    std::string result;
    while(!feof(pipe)) {
        if(fgets(buffer, bufferSize, pipe) != nullptr) {
            result += buffer;
        }
    }
    pclose(pipe);
    return result;
}

TEST(FastaFileReader, Minimal) {  // NOLINT(cert-err58-cpp)
    std::string gitRootDir = exec("git rev-parse --show-toplevel");
    util::rtrim(gitRootDir);

    util::FastaFileReader reader(gitRootDir + "/resources/test-files/fasta/minimal.fasta");

    std::vector<util::FastaRecord> records;
    reader.parse(&records);

    EXPECT_EQ(records.size(), 1);
    EXPECT_EQ(records.front().header, ">minimal");
    EXPECT_EQ(records.front().sequence, "GATTACA");
}

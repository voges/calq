/**
 * @file fasta-file-reader.h
 */

#ifndef CALQ_FASTA_FILE_READER_H_
#define CALQ_FASTA_FILE_READER_H_

#include <fstream>
#include <string>
#include <vector>
#include "fasta-record.h"

namespace calq {

class FastaFileReader {
   public:
    FastaFileReader() = delete;
    explicit FastaFileReader(const std::string& path);
    FastaFileReader(const FastaFileReader&) = delete;
    FastaFileReader& operator=(const FastaFileReader&) = delete;
    FastaFileReader(FastaFileReader&&) = delete;
    FastaFileReader& operator=(FastaFileReader&&) = delete;
    ~FastaFileReader();

    void parse(std::vector<FastaRecord>* fastaRecords);

   private:
    std::ifstream ifs_;
};

}  // namespace calq

#endif  // CALQ_FASTA_FILE_READER_H_

/**
 * @file fasta-file-reader.h
 */

#ifndef CALQ_FASTA_FILE_READER_H_
#define CALQ_FASTA_FILE_READER_H_

#include <string>
#include <vector>
#include "fasta-record.h"
#include "file-reader.h"

namespace calq {

class FastaFileReader : public FileReader {
   public:
    explicit FastaFileReader(const std::string &path);
    void parse(std::vector<FastaRecord> *fastaRecords);
};

}  // namespace calq

#endif  // CALQ_FASTA_FILE_READER_H_

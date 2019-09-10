#ifndef CALQ_FASTA_FILE_READER_H_
#define CALQ_FASTA_FILE_READER_H_

#include <string>
#include <vector>

#include "fasta-record.h"
#include "file-reader.h"

namespace calq {

class FastaFileReader : public FileReader {
   public:
    FastaFileReader(const std::string &path);

    ~FastaFileReader();

    void parse(std::vector<FastaRecord> *const fastaRecords);
};

}  // namespace calq

#endif  // CALQ_FASTA_FILE_READER_H_

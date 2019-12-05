#ifndef CALQ_FASTA_RECORD_H_
#define CALQ_FASTA_RECORD_H_

#include <string>

namespace calq {

struct FastaRecord {
   public:
    FastaRecord() = delete;
    FastaRecord(std::string head, std::string seq) : header(std::move(head)), sequence(std::move(seq)) {}
    FastaRecord(const FastaRecord &) = default;
    FastaRecord &operator=(const FastaRecord &) = default;
    FastaRecord(FastaRecord &&) = delete;
    FastaRecord &operator=(FastaRecord &&) = delete;
    ~FastaRecord() = default;

    std::string header;
    std::string sequence;
};

}  // namespace calq

#endif  // CALQ_FASTA_RECORD_H_

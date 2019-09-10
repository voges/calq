#ifndef CALQ_FASTA_RECORD_H_
#define CALQ_FASTA_RECORD_H_

#include <string>

namespace calq {

struct FastaRecord {
   public:
    FastaRecord(std::string head, std::string seq) : header(std::move(head)), sequence(std::move(seq)) {}

    std::string header;
    std::string sequence;
};

}  // namespace calq

#endif  // CALQ_FASTA_RECORD_H_

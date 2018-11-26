#ifndef CALQ_QUAL_DECODER_H_
#define CALQ_QUAL_DECODER_H_

#include <map>
#include <string>
#include <vector>

#include "calq/cq_file.h"
#include "calq/sam_record.h"

namespace calq {

class QualDecoder {
 public:
    QualDecoder();
    ~QualDecoder();

    void decodeMappedRecordFromBlock(const SAMRecord &samRecord, File *qualFile);
    void decodeUnmappedRecordFromBlock(const SAMRecord &samRecord, File *qualFile);
    size_t readBlock(CQFile *cqFile);

 private:
    uint32_t posOffset_;
    int qualityValueOffset_;

    std::string uqv_;
    std::string qvci_;
    std::vector< std::string > qvi_;

    size_t uqvIdx_;
    std::vector< size_t > qviIdx_;

    std::map<int, Quantizer> quantizers_;
};

}  // namespace calq

#endif  // CALQ_QUAL_DECODER_H_

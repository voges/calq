/** @file QualDecoder.h
 *  @brief This file contains the definition of the QualDecoder class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#ifndef CALQ_QUALCODEC_QUALDECODER_H_
#define CALQ_QUALCODEC_QUALDECODER_H_

#include <map>
#include <string>

#include "IO/CQ/CQFile.h"
#include "IO/SAM/SAMRecord.h"

namespace calq {

class QualDecoder {
public:
    explicit QualDecoder(const std::map<int, Quantizer> &quantizers);
    ~QualDecoder(void);

    size_t decodeMappedRecordFromBlock(const SAMRecord &samRecord, File *qualFile);
    size_t decodeUnmappedRecordFromBlock(const SAMRecord &samRecord, File *qualFile);
    size_t readBlock(CQFile *cqFile);

private:
    uint32_t posOffset_;
    int qualityValueOffset_;

    std::string unmappedQualityValues_;
    std::string mappedQuantizerIndices_;
    std::string mappedQualityValueIndices_;

    std::map<int, Quantizer> quantizers_;
};

} // namespace calq

#endif // CALQ_QUALCODEC_QUALDECODER_H_


/**
 * @file qual-decoder.h
 */

#ifndef CALQ_QUAL_DECODER_H_
#define CALQ_QUAL_DECODER_H_

#include <map>
#include <string>
#include <vector>
#include "data-structures.h"
#include "quantizer.h"

namespace calq {

struct DecodingRead {
    uint32_t posMin;
    std::string cigar;
};

class QualDecoder {
   public:
    QualDecoder() = delete;
    QualDecoder(const DecodingBlock &in, uint32_t positionOffset, uint8_t qualityOffset, EncodingBlock *out);
    QualDecoder(const QualDecoder &) = delete;
    QualDecoder &operator=(const QualDecoder &) = delete;
    QualDecoder(QualDecoder &&) = delete;
    QualDecoder &operator=(QualDecoder &&) = delete;
    ~QualDecoder() = default;

    void decodeMappedRecordFromBlock(const DecodingRead &samRecord);

   private:
    uint32_t posOffset_;
    int qualityValueOffset_;
    std::vector<size_t> qviIdx_;
    std::vector<Quantizer> quantizers_;
    EncodingBlock *out_;
    const DecodingBlock &in_;
};

}  // namespace calq

#endif  // CALQ_QUAL_DECODER_H_

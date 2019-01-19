#ifndef CALQ_QUAL_DECODER_H_
#define CALQ_QUAL_DECODER_H_

// -----------------------------------------------------------------------------

#include <map>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------

#include "calq/quantizer.h"

// -----------------------------------------------------------------------------

namespace calq {

// -----------------------------------------------------------------------------

struct DecodingRead;
struct EncodingBlock;
struct DecodingBlock;

// -----------------------------------------------------------------------------

class QualDecoder
{
 public:
    QualDecoder(const DecodingBlock& in,
                EncodingBlock *out
    );
    ~QualDecoder();

    void decodeMappedRecordFromBlock(const DecodingRead& samRecord);

 private:
    uint32_t posOffset_;
    int qualityValueOffset_;

    size_t uqvIdx_;
    std::vector<size_t> qviIdx_;

    std::map<int, Quantizer> quantizers_;

    EncodingBlock *out;
    const DecodingBlock& in;
};

// -----------------------------------------------------------------------------

}  // namespace calq

// -----------------------------------------------------------------------------

#endif  // CALQ_QUAL_DECODER_H_

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
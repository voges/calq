#include "calq/qual_decoder.h"
#include "calq/structs.h"

#include "calq/error_exception_reporter.h"

namespace calq {

QualDecoder::QualDecoder(const DecodingBlock& b, EncodingBlock* o)
        : posOffset_(0),
        qualityValueOffset_(0),
        uqvIdx_(0),
        qviIdx_(b.stepindices.size(), 0),
        quantizers_(),
        out(o),
        in(b){
    out->qvalues.clear();
}

QualDecoder::~QualDecoder() = default;

void QualDecoder::decodeMappedRecordFromBlock(const DecodingRead& samRecord){
    std::string qual;

    size_t cigarIdx = 0;
    size_t cigarLen = samRecord.cigar.length();
    size_t opLen = 0;
    size_t qvciPos = samRecord.posMin - posOffset_;

    for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++) {
        if (isdigit(samRecord.cigar[cigarIdx])) {
            opLen = opLen * 10 + (size_t) samRecord.cigar[cigarIdx] - (size_t) '0';
            continue;
        }

        switch (samRecord.cigar[cigarIdx]) {
            case 'M':
            case '=':
            case 'X':
                // Decode opLen quality value indices with computed quantizer indices
                for (size_t i = 0; i < opLen; i++) {
                    int quantizerIndex = in.quantizerIndices[qvciPos++] - '0';
                    int qualityValueIndex = in.stepindices.at(static_cast<size_t>(quantizerIndex))[qviIdx_[quantizerIndex]++] - '0';
                    int q = quantizers_.at(quantizerIndex).indexToReconstructionValue(qualityValueIndex);
                    qual += static_cast<char>(q + qualityValueOffset_);
                }
                break;
            case 'I':
            case 'S':
                // Decode opLen quality values with max quantizer index
                for (size_t i = 0; i < opLen; i++) {
                    int qualityValueIndex = in.stepindices.at(quantizers_.size() - 1)[qviIdx_[quantizers_.size() - 1]++] - '0';
                    int q = quantizers_.at(
                            static_cast<const int &>(quantizers_.size() - 1)).indexToReconstructionValue(qualityValueIndex);
                    qual += static_cast<char>(q + qualityValueOffset_);
                }
                break;
            case 'D':
            case 'N':
                qvciPos += opLen;
                break;  // do nothing as these bases are not present
            case 'H':
            case 'P':
                break;  // these have been clipped
            default:
                throwErrorException("Bad CIGAR string");
        }
        opLen = 0;
    }
    out->qvalues.push_back(std::move(qual));
}

}  // namespace calq

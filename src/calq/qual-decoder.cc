/**
 * @file qual-decoder.cc
 */

#include "qual-decoder.h"
#include "exceptions.h"

namespace calq {

QualDecoder::QualDecoder(const DecodingBlock& b, uint32_t positionOffset, uint8_t qualityOffset, EncodingBlock* o)
    : posOffset_(positionOffset),
      qualityValueOffset_(qualityOffset),
      qviIdx_(b.quantizerStepIndexes.size(), 0),
      quantizers_(0),
      out_(o),
      in_(b) {
    out_->qualityValues.clear();
    for (const auto& q : b.codebooks) {
        std::map<int, int> steps;
        for (size_t i = 0; i < q.size(); ++i) {
            steps[i] = q[i];
        }
        quantizers_.emplace_back(steps);
    }
}

QualDecoder::~QualDecoder() = default;

void QualDecoder::decodeMappedRecordFromBlock(const DecodingRead& samRecord) {
    std::string qual;

    size_t cigarIdx = 0;
    size_t cigarLen = samRecord.cigar.length();
    size_t opLen = 0;
    size_t qvciPos = samRecord.posMin - posOffset_;

    for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++) {
        if (isdigit(samRecord.cigar[cigarIdx])) {
            opLen = opLen * 10 + (size_t)samRecord.cigar[cigarIdx] - (size_t)'0';
            continue;
        }

        switch (samRecord.cigar[cigarIdx]) {
            case 'M':
            case '=':
            case 'X':
                // Decode opLen quality value indices with computed
                // quantizer indices
                for (size_t i = 0; i < opLen; i++) {
                    uint8_t quantizerIndex = in_.quantizerIndexes[qvciPos++] - '0';

                    uint8_t qualityValueIndex =
                        in_.quantizerStepIndexes.at(static_cast<size_t>(quantizerIndex))[qviIdx_[quantizerIndex]++] -
                        '0';

                    auto q = uint8_t(quantizers_.at(quantizerIndex).indexToReconstructedValue(qualityValueIndex));

                    qual += static_cast<char>(q + qualityValueOffset_);
                }
                break;
            case 'I':
            case 'S':
                // Decode opLen quality values with max quantizer index
                for (size_t i = 0; i < opLen; i++) {
                    int qualityValueIndex =
                        in_.quantizerStepIndexes.at(quantizers_.size() - 1)[qviIdx_[quantizers_.size() - 1]++] - '0';
                    int q = quantizers_.at(quantizers_.size() - 1).indexToReconstructedValue(qualityValueIndex);
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
    out_->qualityValues.push_back(std::move(qual));
}

}  // namespace calq

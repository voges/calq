/**
 * @file quantizer.cc
 */

#include "quantizer.h"
#include <utility>
#include "errors.h"

namespace calq {

Quantizer::Quantizer() : lut_(), inverseLut_() {}

Quantizer::Quantizer(std::map<int, int> inverseLut) : lut_(), inverseLut_(std::move(inverseLut)) {}

Quantizer::~Quantizer() = default;

int Quantizer::valueToIndex(const int value) const {
    if (lut_.find(value) == lut_.end()) {
        throwErrorException("Value out of range");
    }
    return lut_.at(value).first;
}

int Quantizer::indexToReconstructionValue(const int index) const {
    if (inverseLut_.find(index) == inverseLut_.end()) {
        throwErrorException("Index not found");
    }
    return inverseLut_.at(index);
}

const std::map<int, int>& Quantizer::inverseLut() const { return inverseLut_; }

}  // namespace calq

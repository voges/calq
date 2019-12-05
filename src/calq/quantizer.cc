#include "quantizer.h"
#include <utility>

namespace calq {

Quantizer::Quantizer() : lut_(), inverseLut_() {}

Quantizer::Quantizer(std::map<int, int> inverseLut) : lut_(), inverseLut_(std::move(inverseLut)) {}

int Quantizer::valueToIndex(const int value) const { return lut_.at(value).first; }

int Quantizer::indexToReconstructedValue(const int index) const { return inverseLut_.at(index); }

const std::map<int, int>& Quantizer::inverseLut() const { return inverseLut_; }

}  // namespace calq

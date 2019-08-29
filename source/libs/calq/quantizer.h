/**
 * @file quantizer.h
 */

#ifndef CALQ_QUANTIZER_H_
#define CALQ_QUANTIZER_H_

#include <map>
#include <utility>

namespace calq {

class Quantizer {
   public:
    Quantizer();
    explicit Quantizer(std::map<int, int> inverseLut);
    virtual ~Quantizer();
    int valueToIndex(int value) const;
    int indexToReconstructionValue(int index) const;
    const std::map<int, int>& inverseLut() const;

   protected:
    std::map<int, std::pair<int, int>> lut_;  // value->(index,reconstructionValue)
    std::map<int, int> inverseLut_;           // index->reconstructionValue
};

}  // namespace calq

#endif  // CALQ_QUANTIZER_H_

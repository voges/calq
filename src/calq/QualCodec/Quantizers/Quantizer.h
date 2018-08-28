/** @file Quantizer.h
 *  @brief This file contains the definition of the Quantizer base class.
 */

// Copyright 2015-2017 Leibniz Universitaet Hannover

#ifndef CALQ_QUALCODEC_QUANTIZERS_QUANTIZER_H_
#define CALQ_QUALCODEC_QUANTIZERS_QUANTIZER_H_

#include <map>
#include <utility>

namespace calq {

class Quantizer {
 public:
    Quantizer();
    explicit Quantizer(const std::map<int, int> &inverseLut);
    virtual ~Quantizer();

    int valueToIndex(const int &value) const;
    int indexToReconstructionValue(const int &index) const;
    int valueToReconstructionValue(const int &value) const;

    const std::map<int, int> &inverseLut() const;

    void print() const;

 protected:
    std::map<int, std::pair<int, int>> lut_;  // value->(index,reconstructionValue)
    std::map<int, int> inverseLut_;  // index->reconstructionValue
};

}  // namespace calq

#endif  // CALQ_QUALCODEC_QUANTIZERS_QUANTIZER_H_


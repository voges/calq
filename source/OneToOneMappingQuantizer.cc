/** @file OneToOneMappingQuantizer.cc
 *  @brief This file contains the implementation of the OneToOneMappingQuantizer class.
 */

// Copyright 2015-2017 Leibniz Universitaet Hannover

#include "OneToOneMappingQuantizer.h"

#include <math.h>

#include <queue>
#include <utility>

#include "ErrorExceptionReporter.h"

namespace calq {

OneToOneMappingQuantizer::OneToOneMappingQuantizer(int valueMin, int valueMax) : Quantizer() {
    if ((valueMin > valueMax)) {
        throwErrorException("Error in quantizer initialization");
    }

    // Fill the quantization table
    int index = 0;
    int reconstructionValue = 0;
    for (int value = valueMin; value <= valueMax; ++value) {
        std::pair<int, int> curr(index, reconstructionValue);
        lut_.insert(std::pair<int, std::pair<int, int>>(value, curr));
        inverseLut_.insert(curr);
        index++;
        reconstructionValue++;
    }
}

OneToOneMappingQuantizer::~OneToOneMappingQuantizer(void) {}

}  // namespace calq

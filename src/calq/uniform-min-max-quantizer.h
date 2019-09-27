/**
 * @file uniform-min-max-quantizer.h
 */

#ifndef CALQ_UNIFORM_MIN_MAX_QUANTIZER_H_
#define CALQ_UNIFORM_MIN_MAX_QUANTIZER_H_

#include "uniform-quantizer.h"

namespace calq {

class UniformMinMaxQuantizer : public UniformQuantizer {
   public:
    UniformMinMaxQuantizer(int minValue, int maxValue, int numSteps);
};

}  // namespace calq

#endif  // CALQ_UNIFORM_MIN_MAX_QUANTIZER_H_

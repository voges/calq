/**
 * @file uniform-min-max-quantizer.h
 */

#ifndef CALQ_UNIFORM_MIN_MAX_QUANTIZER_H_
#define CALQ_UNIFORM_MIN_MAX_QUANTIZER_H_

#include "uniform-quantizer.h"

namespace calq {

class UniformMinMaxQuantizer : public UniformQuantizer {
   public:
    UniformMinMaxQuantizer() = delete;
    UniformMinMaxQuantizer(int minValue, int maxValue, int numSteps);
    UniformMinMaxQuantizer(const UniformMinMaxQuantizer &) = delete;
    UniformMinMaxQuantizer &operator=(const UniformMinMaxQuantizer &) = delete;
    UniformMinMaxQuantizer(UniformMinMaxQuantizer &&) = delete;
    UniformMinMaxQuantizer &operator=(UniformMinMaxQuantizer &&) = delete;
    ~UniformMinMaxQuantizer() override = default;
};

}  // namespace calq

#endif  // CALQ_UNIFORM_MIN_MAX_QUANTIZER_H_

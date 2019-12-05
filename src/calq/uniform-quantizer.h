#ifndef CALQ_UNIFORM_QUANTIZER_H_
#define CALQ_UNIFORM_QUANTIZER_H_

#include "quantizer.h"

namespace calq {

class UniformQuantizer : public Quantizer {
   public:
    UniformQuantizer() = delete;
    UniformQuantizer(int minValue, int maxValue, int numSteps);
    UniformQuantizer(const UniformQuantizer &) = delete;
    UniformQuantizer &operator=(const UniformQuantizer &) = delete;
    UniformQuantizer(UniformQuantizer &&) = delete;
    UniformQuantizer &operator=(UniformQuantizer &&) = delete;
    ~UniformQuantizer() override = default;
};

}  // namespace calq

#endif  // CALQ_UNIFORM_QUANTIZER_H_

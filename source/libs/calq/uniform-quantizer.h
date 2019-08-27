#ifndef CALQ_UNIFORM_QUANTIZER_H_
#define CALQ_UNIFORM_QUANTIZER_H_

#include "quantizer.h"

namespace calq {

class UniformQuantizer : public Quantizer {
   public:
    UniformQuantizer(int minValue, int maxValue, int numSteps);
};

}  // namespace calq

#endif  // CALQ_UNIFORM_QUANTIZER_H_

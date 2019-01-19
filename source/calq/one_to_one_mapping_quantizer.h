#ifndef CALQ_ONE_TO_ONE_MAPPING_QUANTIZER_H_
#define CALQ_ONE_TO_ONE_MAPPING_QUANTIZER_H_

// -----------------------------------------------------------------------------

#include "calq/quantizer.h"

// -----------------------------------------------------------------------------

namespace calq {

// -----------------------------------------------------------------------------

class OneToOneMappingQuantizer : public Quantizer
{
 public:
    OneToOneMappingQuantizer(int valueMin,
                             int valueMax
    );
    ~OneToOneMappingQuantizer(void);
};

// -----------------------------------------------------------------------------

}  // namespace calq

// -----------------------------------------------------------------------------

#endif  // CALQ_ONE_TO_ONE_MAPPING_QUANTIZER_H_

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
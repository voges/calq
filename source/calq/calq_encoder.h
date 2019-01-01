#ifndef CALQ_CALQ_ENCODER_H_
#define CALQ_CALQ_ENCODER_H_

#include "calq/structs.h"

namespace calq {

void encode(const EncodingOptions& opt,
            const EncodingSideInformation& sideInformation,
            const EncodingBlock& input,
            DecodingBlock *output
);

}  // namespace calq

#endif  // CALQ_CALQ_ENCODER_H_
#ifndef CALQ_CALQ_DECODER_H_
#define CALQ_CALQ_DECODER_H_

#include "calq/structs.h"

namespace calq {

void decode(const DecodingSideInformation& sideInformation,
            const DecodingBlock& input,
            DecodingBlock* output
);

}  // namespace calq

#endif  // CALQ_CALQ_DECODER_H_

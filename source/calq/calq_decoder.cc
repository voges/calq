#include "calq/calq_decoder.h"

#include <chrono>

#include "calq/qual_decoder.h"
#include "calqapp/sam_file.h"
#include "calq/error_exception_reporter.h"
#include "calq/log.h"
#include "calq/structs.h"

namespace calq {

void decode(const DecodingSideInformation& sideInformation,
            const DecodingBlock& input,
            EncodingBlock *output
){

    // Decode the quality values
    QualDecoder qualDecoder(input, output);
    output->qvalues.clear();
    for (size_t i = 0; i < sideInformation.positions.size(); ++i)
    {
        DecodingRead r = {sideInformation.positions[i], sideInformation.cigars[i]};
        qualDecoder.decodeMappedRecordFromBlock(r);
    }
}

}  // namespace calq

#include "calq/calq_decoder.h"

#include <chrono>

#include "calq/qual_decoder.h"
#include "calq/sam_file.h"
#include "calq/error_exception_reporter.h"
#include "calq/log.h"
#include "calq/structs.h"

namespace calq {

void decode(const DecodingSideInformation& sideInformation,
            const DecodingBlock& input,
            EncodingBlock *output
){

    // Decode the quality values
    QualDecoder qualDecoder;
    for (size_t i = 0; i < sideInformation.positions.size(); ++i)
    {
        DecodingRead r;// = {};
        qualDecoder.decodeMappedRecordFromBlock(r);
    }
}

}  // namespace calq

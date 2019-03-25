#include "calq/calq_decoder.h"

#include <chrono>

#include "calq/exceptions.h"
#include "calq/log.h"
#include "calq/qual_decoder.h"

namespace calq {

CalqDecoder::CalqDecoder(const Options &options)
    : cqFile_(options.inputFileName, CQFile::MODE_READ),
      qualFile_(options.outputFileName, File::MODE_WRITE),
      sideInformationFile_(options.sideInformationFileName) {
    if (options.inputFileName.empty() == true) {
        throwErrorException("options.inputFileName is empty");
    }
    if (options.outputFileName.empty() == true) {
        throwErrorException("options.outputFileName is empty");
    }
    if (options.sideInformationFileName.empty() == true) {
        throwErrorException("options.sideInformationFileName is empty");
    }
}

CalqDecoder::~CalqDecoder(void) {}

void CalqDecoder::decode(void) {
    // Take time
    auto startTime = std::chrono::steady_clock::now();

    // Read CQ file header
    CALQ_LOG("Reading CQ file header");
    size_t blockSize = 0;
    cqFile_.readHeader(&blockSize);

    while (sideInformationFile_.readBlock(blockSize) != 0) {
//         CALQ_LOG("Decoding block %zu", sideInformationFile_.nrBlocksRead()-1);

        // Decode the quality values
        QualDecoder qualDecoder;
        qualDecoder.readBlock(&cqFile_);
        for (auto const &samRecord : sideInformationFile_.currentBlock.records) {
            if (samRecord.isMapped() == true) {
                qualDecoder.decodeMappedRecordFromBlock(samRecord, &qualFile_);
            } else {
                qualDecoder.decodeUnmappedRecordFromBlock(samRecord, &qualFile_);
            }
        }
    }

    auto stopTime = std::chrono::steady_clock::now();
    auto diffTime = stopTime - startTime;
    auto diffTimeMs = std::chrono::duration_cast<std::chrono::milliseconds>(diffTime).count();
    auto diffTimeS = std::chrono::duration_cast<std::chrono::seconds>(diffTime).count();
    auto diffTimeM = std::chrono::duration_cast<std::chrono::minutes>(diffTime).count();
    auto diffTimeH = std::chrono::duration_cast<std::chrono::hours>(diffTime).count();

    CALQ_LOG("DECOMPRESSION STATISTICS");
    CALQ_LOG("  Took %d ms ~= %d s ~= %d m ~= %d h",
        static_cast<int>(diffTimeMs),
        static_cast<int>(diffTimeS),
        static_cast<int>(diffTimeM),
        static_cast<int>(diffTimeH)
    );
    CALQ_LOG("  Speed (compressed size/time): %.2f MB/s",
        (static_cast<double>(cqFile_.nrReadBytes()) / static_cast<double>(MB)) / static_cast<double>(diffTimeS)
    );
    CALQ_LOG("  Decoded %zu block(s)", sideInformationFile_.nrBlocksRead());
}

}  // namespace calq

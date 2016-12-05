/** @file CalqDecoder.cc
 *  @brief This file contains the implementation of the CalqDecoder class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#include "CalqDecoder.h"

#include <chrono>

#include "Common/Exceptions.h"
#include "Common/log.h"
#include "QualCodec/QualDecoder.h"

namespace calq {

CalqDecoder::CalqDecoder(const Options &options)
    : cqFile_(options.inputFileName, CQFile::MODE_READ),
      qualFile_(options.outputFileName, File::MODE_WRITE),
      sideInformationFile_(options.sideInformationFileName)
{
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

void CalqDecoder::decode(void)
{
    // Take time
    auto startTime = std::chrono::steady_clock::now();

    // Read CQ file header
    CALQ_LOG("Reading CQ file header");
    size_t blockSize = 0;
    cqFile_.readHeader(&blockSize);

    while (sideInformationFile_.readBlock(blockSize) != 0) {
        CALQ_LOG("Decoding block %zu", sideInformationFile_.nrBlocksRead()-1);

        // Read the inverse quantization LUTs
        CALQ_LOG("Reading inverse quantization LUTs");
        std::map<int,Quantizer> quantizers;
        cqFile_.readQuantizers(&quantizers);

        // Decode the QVs
        CALQ_LOG("Decoding quality values");
        QualDecoder qualDecoder(quantizers);
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
    CALQ_LOG("  Took %d ms ~= %d s ~= %d m ~= %d h", (int)diffTimeMs, (int)diffTimeS, (int)diffTimeM, (int)diffTimeH);
    CALQ_LOG("  Speed (compressed size/time): %.2f MB/s", ((double)(cqFile_.nrReadBytes()/MB))/((double)diffTimeS));
    CALQ_LOG("  Decoded %zu block(s)", sideInformationFile_.nrBlocksRead());
}

} // namespace calq

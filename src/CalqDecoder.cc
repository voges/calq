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

    QualDecoder qualDecoder;
    while (sideInformationFile_.readBlock(blockSize) != 0) {
        CALQ_LOG("Decoding block %zu", sideInformationFile_.nrBlocksRead()-1);

        // Read the inverse quantization LUTs
        CALQ_LOG("Reading inverse quantization LUTs");
        std::map<int,Quantizer> quantizers;
        cqFile_.readQuantizers(quantizers);
        
        
        uint8_t nrQuantizers = 0;
        cqFile_.readUint8(&nrQuantizers);
        CALQ_LOG("Reading %u quantizers", nrQuantizers);
        for (int i = 0; i < nrQuantizers; ++i) {
            CALQ_LOG("Reading quantizer %d", i);
            uint8_t quantizerIdx = 0;
            cqFile_.readUint8(&quantizerIdx);
            CALQ_LOG("  Quantizer index: %d", quantizerIdx);
            uint8_t quantizerSteps = 0;
            cqFile_.readUint8(&quantizerSteps);
            CALQ_LOG("  Quantizer steps: %d", quantizerSteps);
            for (int s = 0; s < quantizerSteps; ++s) {
                uint8_t index = 0;
                cqFile_.readUint8(&index);
                uint8_t reconstructionValue = 0;
                cqFile_.readUint8(&reconstructionValue);
                CALQ_LOG("  %u -> %u", index, reconstructionValue);
            }
        }
        
        
        // Decode the QVs
        CALQ_LOG("Decoding quality values");
        qualDecoder.readBlock(cqFile_);
        for (auto const &samRecord : sideInformationFile_.currentBlock.records) {
            if (samRecord.isMapped() == true) {

            } else {
                
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
    CALQ_LOG("  Decoded %zu block(s)", sideInformationFile_.nrBlocksRead());
    CALQ_LOG("  Compressed size: %zu", cqFile_.nrReadBytes());
    CALQ_LOG("  Speed (compressed size/time): %.2f MB/s", ((double)(cqFile_.nrReadBytes()/MB))/(double)((double)diffTimeMs/1000));
}

} // namespace calq

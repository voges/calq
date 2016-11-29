/** @file CalqDecoder.cc
 *  @brief This file contains the implementation of the CalqDecoder class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "CalqDecoder.h"
#include "Common/Exceptions.h"
#include <chrono>

cq::CalqDecoder::CalqDecoder(const CLIOptions &cliOptions)
    : cqFile_(cliOptions.inputFileName, CQFile::MODE_READ)
    , qualFile_(cliOptions.outputFileName, File::MODE_WRITE)
    , sideInformationFile_(cliOptions.sideInformationFileName)
{
    if (cliOptions.inputFileName.empty() == true) {
        throwErrorException("cliOptions.inputFileName is empty");
    }
    if (cliOptions.outputFileName.empty() == true) {
        throwErrorException("cliOptions.outputFileName is empty");
    }
    if (cliOptions.sideInformationFileName.empty() == true) {
        throwErrorException("cliOptions.sideInformationFileName is empty");
    }
}

cq::CalqDecoder::~CalqDecoder(void)
{
    // empty
}

void cq::CalqDecoder::decode(void)
{
    auto startTime = std::chrono::steady_clock::now();
    size_t blockSize = 0;
    cqFile_.readHeader(&blockSize);

    while (sideInformationFile_.readBlock(blockSize) != 0) {
//         for (auto const &samRecord : sideInformationFile_.currentBlock.records) {
//             // empty
//         }
    }

    auto stopTime = std::chrono::steady_clock::now();
    auto diffTime = stopTime - startTime;
    auto diffTimeMs = std::chrono::duration_cast<std::chrono::milliseconds>(diffTime).count();
    auto diffTimeS = std::chrono::duration_cast<std::chrono::seconds>(diffTime).count();
    auto diffTimeM = std::chrono::duration_cast<std::chrono::minutes>(diffTime).count();
    auto diffTimeH = std::chrono::duration_cast<std::chrono::hours>(diffTime).count();

    CQ_LOG("DECOMPRESSION STATISTICS");
    CQ_LOG("  Took %d ms ~= %d s ~= %d m ~= %d h", (int)diffTimeMs, (int)diffTimeS, (int)diffTimeM, (int)diffTimeH);
    CQ_LOG("  Decoded %zu block(s)", sideInformationFile_.nrBlocksRead());
    CQ_LOG("  Compressed size: %zu", cqFile_.nrReadBytes());
    CQ_LOG("  Speed (compressed size/time): %.2f MB/s", ((double)(cqFile_.nrReadBytes()/MB))/(double)((double)diffTimeMs/1000));
}


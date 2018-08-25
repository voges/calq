/** @file CalqDecoder.cc
 *  @brief This file contains the implementation of the CalqDecoder class.
 */

// Copyright 2015-2017 Leibniz Universitaet Hannover

#include "CalqDecoder.h"

#include <chrono>

#include "QualCodec/QualDecoder.h"
#include "IO/SAM/SAMFile.h"

namespace calq {

CalqDecoder::CalqDecoder(const Options &options)
        : cqFile_(nullptr),
          qualFile_(nullptr),
          sideInformationFile_(nullptr) {
    if (options.inputFileName.empty()) {
        throwErrorException("options.inputFileName is empty");
    }
    if (options.outputFileName.empty()) {
        throwErrorException("options.outputFileName is empty");
    }
    if (options.sideInformationFileName.empty()) {
        throwErrorException("options.sideInformationFileName is empty");
    }

    try {
        cqFile_ = new CQFile(options.inputFileName, CQFile::MODE_READ);
        qualFile_ = new File(options.outputFileName, File::MODE_WRITE);
        sideInformationFile_ = new SAMFile(options.sideInformationFileName);
    } catch (const Exception &e) {
        if (cqFile_ != nullptr) {
            delete cqFile_;
            cqFile_ = nullptr;
        }
        if (qualFile_ != nullptr) {
            delete qualFile_;
            qualFile_ = nullptr;
        }
        if (sideInformationFile_ != nullptr) {
            delete sideInformationFile_;
            sideInformationFile_ = nullptr;
        }
        throw e;
    }
}

CalqDecoder::~CalqDecoder() {
    delete cqFile_;
    delete qualFile_;
    delete sideInformationFile_;
}

void CalqDecoder::decode() {
    // Take time
    auto startTime = std::chrono::steady_clock::now();

    // Read CQ file header
    CALQ_LOG("Reading CQ file header");
    size_t blockSize = 0;
    cqFile_->readHeader(&blockSize);

    while (sideInformationFile_->readBlock(blockSize) != 0) {
//         CALQ_LOG("Decoding block %zu", sideInformationFile_->nrBlocksRead()-1);

        // Decode the quality values
        QualDecoder qualDecoder;
        qualDecoder.readBlock(cqFile_);
        for (auto const &samRecord : sideInformationFile_->currentBlock.records) {
            if (samRecord.isMapped()) {
                qualDecoder.decodeMappedRecordFromBlock(samRecord, qualFile_);
            } else {
                qualDecoder.decodeUnmappedRecordFromBlock(samRecord, qualFile_);
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
    CALQ_LOG("  Took %d ms ~= %d s ~= %d m ~= %d h", (int) diffTimeMs, (int) diffTimeS, (int) diffTimeM, (int) diffTimeH);
    CALQ_LOG("  Speed (compressed size/time): %.2f MB/s",
             ((static_cast<double>(cqFile_->nrReadBytes()) / static_cast<double>(MB))) / (static_cast<double>(diffTimeS)));
    CALQ_LOG("  Decoded %zu block(s)", sideInformationFile_->nrBlocksRead());
}

}  // namespace calq


#include "calq/calq_decoder.h"

#include <chrono>

#include "calq/qual_decoder.h"
#include "calq/sam_file.h"
#include "calq/error_exception_reporter.h"
#include "calq/log.h"

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

    cqFile_ = calq::make_unique<CQFile>(options.inputFileName, CQFile::Mode::MODE_READ);
    qualFile_ = calq::make_unique<File>(options.outputFileName, File::Mode::MODE_WRITE);
    sideInformationFile_ = calq::make_unique<SAMFile>(options.sideInformationFileName);
}

CalqDecoder::~CalqDecoder() = default;

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
        qualDecoder.readBlock(cqFile_.get());
        for (auto const &samRecord : sideInformationFile_->currentBlock.records) {
            if (samRecord.isMapped()) {
                qualDecoder.decodeMappedRecordFromBlock(samRecord, qualFile_.get());
            } else {
                qualDecoder.decodeUnmappedRecordFromBlock(samRecord, qualFile_.get());
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
        (int) diffTimeMs, (int) diffTimeS, (int) diffTimeM, (int) diffTimeH
    );
    CALQ_LOG("  Speed (compressed size/time): %.2f MB/s",
             ((static_cast<double>(cqFile_->nrReadBytes()) / static_cast<double>(MB))) / (static_cast<double>(diffTimeS)));
    CALQ_LOG("  Decoded %zu block(s)", sideInformationFile_->nrBlocksRead());
}

}  // namespace calq

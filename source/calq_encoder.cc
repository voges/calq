/** @file CalqEncoder.cc
 *  @brief This file contains the implementation of the CalqEncoder class.
 */

// Copyright 2015-2017 Leibniz Universitaet Hannover

#include "calq_encoder.h"

#include <chrono>

#include "error_exception_reporter.h"
#include "log.h"
#include "helpers.h"
#include "fasta_file.h"
#include "qual_encoder.h"
#include "probability_distribution.h"
#include "uniform_min_max_quantizer.h"
#include "lloyd_max_quantizer.h"
#include "sam_file.h"


namespace calq {

CalqEncoder::CalqEncoder(const Options &opt) : cqFile_(nullptr),
                                               samFile_(nullptr),
                                               fastaFile_(nullptr),
                                               options(opt) {
    if (opt.blockSize < 1) {
        throwErrorException("blockSize must be greater than zero");
    }
    if (opt.inputFileName.empty()) {
        throwErrorException("inputFileName is empty");
    }
    if (opt.outputFileName.empty()) {
        throwErrorException("outputFileName is empty");
    }
    if (opt.polyploidy < 1) {
        throwErrorException("polyploidy must be greater than zero");
    }
    if (opt.qualityValueOffset < 1) {
        throwErrorException("qualityValueOffset must be greater than zero");
    }
    if (opt.referenceFileNames.empty()) {
        throwErrorException("referenceFileNames is empty");
    }

    cqFile_ = calq::make_unique<CQFile>(opt.outputFileName, CQFile::Mode::MODE_WRITE);
    samFile_ = calq::make_unique<SAMFile>(opt.inputFileName);
    fastaFile_ = calq::make_unique<FASTAFile>(opt.referenceFileNames);
}

CalqEncoder::~CalqEncoder() = default;

void CalqEncoder::encode() {
    size_t compressedMappedQualSize = 0;
    size_t compressedUnmappedQualSize = 0;
    size_t uncompressedMappedQualSize = 0;
    size_t uncompressedUnmappedQualSize = 0;

    // Take time
    auto startTime = std::chrono::steady_clock::now();

    // Write CQ file header
    CALQ_LOG("Writing CQ file header");
    cqFile_->writeHeader(options.blockSize);

    while (samFile_->readBlock(options.blockSize) != 0) {
//         CALQ_LOG("Processing block %zu", samFile_->nrBlocksRead()-1);

        ProbabilityDistribution pdf(options.qualityValueMin, options.qualityValueMax);

        // Check quality value range
        for (auto const &samRecord : samFile_->currentBlock.records) {
            if (samRecord.isMapped()) {
                for (auto const &q : samRecord.qual) {
                    if ((static_cast<int>(q) - options.qualityValueOffset) < options.qualityValueMin) {
                        throwErrorException("Quality value too small");
                    }
                    if ((static_cast<int>(q) - options.qualityValueOffset) > options.qualityValueMax) {
                        throwErrorException("Quality value too large");
                    }
                    pdf.addToPdf((static_cast<int>(q) - options.qualityValueOffset));
                }
            }
        }

        std::map<int, Quantizer> quantizers;

        for (auto i = static_cast<int>(options.quantizationMin); i <= static_cast<int>(options.quantizationMax); ++i) {
            if (options.quantizerType == Options::QuantizerType::UNIFORM) {
                UniformMinMaxQuantizer quantizer(static_cast<const int &>(options.qualityValueMin),
                                                 static_cast<const int &>(options.qualityValueMax), i);
                quantizers.insert(std::pair<int, Quantizer>(static_cast<const int &>(i - options.quantizationMin), quantizer));
            } else if (options.quantizerType == Options::QuantizerType::LLOYD_MAX) {
                LloydMaxQuantizer quantizer(static_cast<size_t>(i));
                quantizer.build(pdf);
                quantizers.insert(std::pair<int, Quantizer>(static_cast<const int &>(i - options.quantizationMin), quantizer));
            } else {
                throwErrorException("Quantization Type not supported");
            }
        }

        // Encode the quality values
        QualEncoder qualEncoder(options, quantizers);
        for (auto const &samRecord : samFile_->currentBlock.records) {
            if (samRecord.isMapped()) {
                qualEncoder.addMappedRecordToBlock(samRecord, *fastaFile_);
            } else {
                qualEncoder.addUnmappedRecordToBlock(samRecord);
            }
        }
        qualEncoder.finishBlock(*fastaFile_, samFile_->currentBlock.records.back().rname);
        qualEncoder.writeBlock(cqFile_.get());

        // Update statistics
        compressedMappedQualSize += qualEncoder.compressedMappedQualSize();
        compressedUnmappedQualSize += qualEncoder.compressedUnmappedQualSize();
        uncompressedMappedQualSize += qualEncoder.uncompressedMappedQualSize();
        uncompressedUnmappedQualSize += qualEncoder.uncompressedUnmappedQualSize();
    }

    auto stopTime = std::chrono::steady_clock::now();
    auto diffTime = stopTime - startTime;
    auto diffTimeMs = std::chrono::duration_cast<std::chrono::milliseconds>(diffTime).count();
    auto diffTimeS = std::chrono::duration_cast<std::chrono::seconds>(diffTime).count();
    auto diffTimeM = std::chrono::duration_cast<std::chrono::minutes>(diffTime).count();
    auto diffTimeH = std::chrono::duration_cast<std::chrono::hours>(diffTime).count();

    CALQ_LOG("COMPRESSION STATISTICS");
    CALQ_LOG("  Took %d ms ~= %d s ~= %d m ~= %d h", (int) diffTimeMs, (int) diffTimeS, (int) diffTimeM, (int) diffTimeH);
    CALQ_LOG("  Speed (uncompressed size/time): %.2f MB/s", ((static_cast<double>(uncompressedMappedQualSize + uncompressedUnmappedQualSize) /
                                                              static_cast<double>(MB))) / (static_cast<double>(diffTimeS)));
    CALQ_LOG("  Wrote %zu block(s)", samFile_->nrBlocksRead());
    CALQ_LOG("  Record(s):  %12zu", samFile_->nrRecordsRead());
    CALQ_LOG("    Mapped:   %12zu", samFile_->nrMappedRecordsRead());
    CALQ_LOG("    Unmapped: %12zu", samFile_->nrUnmappedRecordsRead());
    CALQ_LOG("  Uncompressed size: %12zu", uncompressedMappedQualSize + uncompressedUnmappedQualSize);
    CALQ_LOG("    Mapped:          %12zu", uncompressedMappedQualSize);
    CALQ_LOG("    Unmapped:        %12zu", uncompressedUnmappedQualSize);
    CALQ_LOG("  Compressed size: %12zu", cqFile_->nrWrittenBytes());
    CALQ_LOG("    File format:   %12zu", cqFile_->nrWrittenFileFormatBytes());
    CALQ_LOG("    Mapped:        %12zu", compressedMappedQualSize);
    CALQ_LOG("    Unmapped:      %12zu", compressedUnmappedQualSize);
    CALQ_LOG("  Compression ratio: %4.2f%%", (double) cqFile_->nrWrittenBytes() * 100 / (double) (uncompressedMappedQualSize + uncompressedUnmappedQualSize));
    CALQ_LOG("    Mapped:          %4.2f%%", (double) compressedMappedQualSize * 100 / (double) (uncompressedMappedQualSize));
    CALQ_LOG("    Unmapped:        %4.2f%%", (double) compressedUnmappedQualSize * 100 / (double) (uncompressedUnmappedQualSize));
    CALQ_LOG("  Compression factor: %4.2f", (double) (uncompressedMappedQualSize + uncompressedUnmappedQualSize) / (double) cqFile_->nrWrittenBytes());
    CALQ_LOG("    Mapped:           %4.2f", (double) (uncompressedMappedQualSize) / (double) compressedMappedQualSize);
    CALQ_LOG("    Unmapped:         %4.2f", (double) (uncompressedUnmappedQualSize) / (double) compressedUnmappedQualSize);
    CALQ_LOG("  Bits per quality value: %2.4f",
             (static_cast<double>(cqFile_->nrWrittenBytes()) * 8) / static_cast<double>(uncompressedMappedQualSize + uncompressedUnmappedQualSize));
    CALQ_LOG("    Mapped:               %2.4f", (static_cast<double>(compressedMappedQualSize) * 8) / static_cast<double>(uncompressedMappedQualSize));
    CALQ_LOG("    Unmapped:             %2.4f", (static_cast<double>(compressedUnmappedQualSize) * 8) / static_cast<double>(uncompressedUnmappedQualSize));
}

}  // namespace calq

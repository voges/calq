/** @file CalqEncoder.cc
 *  @brief This file contains the implementation of the CalqEncoder class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "CalqEncoder.h"
#include "Common/constants.h"
#include "Common/Exceptions.h"
#include "Common/helpers.h"
#include "IO/FASTA/FASTAFile.h"
#include "QualCodec/QualEncoder.h"
#include "QualCodec/UniformQuantizer.h"
#include <chrono>
#include <limits>

cq::CalqEncoder::CalqEncoder(const CLIOptions &cliOptions)
    : blockSize_(cliOptions.blockSize)
    , cqFile_(cliOptions.outputFileName, CQFile::MODE_WRITE)
    , polyploidy_(cliOptions.polyploidy)
    , referenceFileNames_(cliOptions.referenceFileNames)
    , samFile_(cliOptions.inputFileName)
{
    if (cliOptions.blockSize < 1) {
        throwErrorException("blockSize must be greater than zero");
    }
    if (cliOptions.inputFileName.empty() == true) {
        throwErrorException("inputFileName is empty");
    }
    if (cliOptions.outputFileName.empty() == true) {
        throwErrorException("outputFileName is empty");
    }
    if (cliOptions.polyploidy < 1) {
        throwErrorException("polyploidy must be greater than zero");
    }
    //if (referenceFileNames.empty() == true) {
    //    throwErrorException("referenceFileNames is empty");
    //}

    // Check and, in case they are provided, get reference sequences
    if (referenceFileNames_.empty() == true) {
        CQ_LOG("No reference file name(s) given - operating without reference sequence(s)");
    } else {
        CQ_LOG("Looking in %zu reference file(s) for reference sequence(s)", referenceFileNames_.size());
        for (auto const &referenceFileName : referenceFileNames_) {
            CQ_LOG("Parsing reference file: %s", referenceFileName.c_str());
            FASTAFile fastaFile(referenceFileName);
            CQ_LOG("Found %zu reference(s):", fastaFile.references.size());
            for (auto const &reference : fastaFile.references) {
                CQ_LOG("  %s (length: %zu)", reference.first.c_str(), reference.second.length());
            }
        }
    }
}

cq::CalqEncoder::~CalqEncoder(void)
{
    // empty
}

void cq::CalqEncoder::encode(void)
{
    auto startTime = std::chrono::steady_clock::now();;
    size_t uncompressedSize = 0;
    cqFile_.writeHeader(blockSize_);

    while (samFile_.readBlock(blockSize_) != 0) {
        CQ_LOG("Processing block %zu", samFile_.nrBlocksRead()-1);

        // Compute min and max QV for the current block
        int qMin = std::numeric_limits<int>::max();
        int qMax = std::numeric_limits<int>::min();
        for (auto const &samRecord : samFile_.currentBlock.records) {
            uncompressedSize += samRecord.qual.length();
            if (samRecord.isMapped() == true) {
                for (auto const &q : samRecord.qual) {
                    if ((int)q < qMin) qMin = q;
                    if ((int)q > qMax) qMax = q;
                }
            }
        }
        CQ_LOG("[qMin,qMax] = [%d,%d]", qMin, qMax);

        // Construct quantizers and store inverse quantization LUTs
        CQ_LOG("Constructing %u quantizers and storing inverse LUTs", QualEncoder::NR_QUANTIZERS);
        std::map<int,Quantizer> quantizers;
        unsigned int quantizerSteps = QualEncoder::QUANTIZER_STEPS_MIN;
        unsigned int quantizerIdx = QualEncoder::QUANTIZER_IDX_MIN;
        cqFile_.writeUint8((uint8_t)QualEncoder::NR_QUANTIZERS);
        for (unsigned int i = 0; i < QualEncoder::NR_QUANTIZERS; i++, quantizerIdx++, quantizerSteps++) {
            Quantizer quantizer = UniformQuantizer(qMin, qMax, quantizerSteps);
            quantizers.insert(std::pair<int,Quantizer>(quantizerIdx, quantizer));

            cqFile_.writeUint8(quantizerIdx);
            cqFile_.writeUint8(quantizerSteps);
            for (auto const &inverseLutEntry : quantizer.inverseLut()) {
                cqFile_.writeUint8((uint8_t)inverseLutEntry.first);
                cqFile_.writeUint8((uint8_t)inverseLutEntry.second);
            }
        }

        // Encode the QVs
        CQ_LOG("Encoding quality values");
        QualEncoder qualEncoder(polyploidy_, qMin, qMax, quantizers);
        qualEncoder.startBlock();
        for (auto const &samRecord : samFile_.currentBlock.records) {
            if (samRecord.isMapped() == true) {
                qualEncoder.addMappedRecordToBlock(samRecord);
            } else {
                qualEncoder.addUnmappedRecordToBlock(samRecord);
            }
        }
        qualEncoder.finishAndWriteBlock(cqFile_);
    }

    auto stopTime = std::chrono::steady_clock::now();
    auto diffTime = stopTime - startTime;
    auto diffTimeMs = std::chrono::duration_cast<std::chrono::milliseconds>(diffTime).count();
    auto diffTimeS = std::chrono::duration_cast<std::chrono::seconds>(diffTime).count();
    auto diffTimeM = std::chrono::duration_cast<std::chrono::minutes>(diffTime).count();
    auto diffTimeH = std::chrono::duration_cast<std::chrono::hours>(diffTime).count();

    CQ_LOG("COMPRESSION STATISTICS");
    CQ_LOG("  Took %d ms ~= %d s ~= %d m ~= %d h", (int)diffTimeMs, (int)diffTimeS, (int)diffTimeM, (int)diffTimeH);
    CQ_LOG("  Compressed %zu mapped + %zu unmapped = %zu record(s) in %zu block(s)", samFile_.nrMappedRecordsRead(), samFile_.nrUnmappedRecordsRead(), samFile_.nrRecordsRead(), samFile_.nrBlocksRead());
    CQ_LOG("  Uncompressed size: %zu", uncompressedSize);
    CQ_LOG("  Compressed size: %zu", cqFile_.nrWrittenBytes());
    CQ_LOG("  Compression ratio: %.2f%%", (double)cqFile_.nrWrittenBytes()*100/(double)uncompressedSize);
    CQ_LOG("  Compression factor: %.2f", (double)uncompressedSize/(double)cqFile_.nrWrittenBytes());
    CQ_LOG("  Bits per quality value: %.4f", ((double)cqFile_.nrWrittenBytes() * 8)/(double)uncompressedSize);
    CQ_LOG("  Speed (uncompressed size/time): %.2f MB/s", ((double)(uncompressedSize/MB))/(double)((double)diffTimeMs/1000));
}


/** @file CalqEncoder.cc
 *  @brief This file contains the implementation of the CalqEncoder class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#include "CalqEncoder.h"

#include <chrono>
#include <limits>

#include "Common/constants.h"
#include "Common/Exceptions.h"
#include "Common/log.h"
#include "IO/FASTA/FASTAFile.h"
#include "QualCodec/QualEncoder.h"
#include "QualCodec/UniformQuantizer.h"

namespace calq {

CalqEncoder::CalqEncoder(const Options &options)
    : blockSize_(options.blockSize),
      cqFile_(options.outputFileName, CQFile::MODE_WRITE),
      polyploidy_(options.polyploidy),
      qualityValueMax_(options.qualityValueMax),
      qualityValueMin_(options.qualityValueMin),
      qualityValueOffset_(options.qualityValueOffset),
      referenceFileNames_(options.referenceFileNames),
      samFile_(options.inputFileName)
{
    if (options.blockSize < 1) {
        throwErrorException("blockSize must be greater than zero");
    }
    if (options.inputFileName.empty() == true) {
        throwErrorException("inputFileName is empty");
    }
    if (options.outputFileName.empty() == true) {
        throwErrorException("outputFileName is empty");
    }
    if (options.polyploidy < 1) {
        throwErrorException("polyploidy must be greater than zero");
    }
    // TODO: check qualityValue*_
    //if (referenceFileNames.empty() == true) {
    //    throwErrorException("referenceFileNames is empty");
    //}

    // Check and, in case they are provided, get reference sequences
    if (referenceFileNames_.empty() == true) {
        CALQ_LOG("No reference file name(s) given - operating without reference sequence(s)");
    } else {
        CALQ_LOG("Looking in %zu reference file(s) for reference sequence(s)", referenceFileNames_.size());
        for (auto const &referenceFileName : referenceFileNames_) {
            CALQ_LOG("Parsing reference file: %s", referenceFileName.c_str());
            FASTAFile fastaFile(referenceFileName);
            CALQ_LOG("Found %zu reference(s):", fastaFile.references.size());
            for (auto const &reference : fastaFile.references) {
                CALQ_LOG("  %s (length: %zu)", reference.first.c_str(), reference.second.length());
            }
        }
    }
}

CalqEncoder::~CalqEncoder(void) {}

void CalqEncoder::encode(void)
{
    auto startTime = std::chrono::steady_clock::now();;
    size_t uncompressedSize = 0;
    cqFile_.writeHeader(blockSize_);

    while (samFile_.readBlock(blockSize_) != 0) {
        CALQ_LOG("Processing block %zu", samFile_.nrBlocksRead()-1);

        // Compute min and max QV for the current block
        for (auto const &samRecord : samFile_.currentBlock.records) {
            uncompressedSize += samRecord.qual.length();
            if (samRecord.isMapped() == true) {
                for (auto const &q : samRecord.qual) {
                    if (q-qualityValueOffset_ < qualityValueMin_) {
                        throwErrorException("Quality value too small");
                    }
                    if (q-qualityValueOffset_ > qualityValueMax_) {
                        throwErrorException("Quality value too large");
                    }
                }
            }
        }

        // Construct quantizers and store inverse quantization LUTs
        CALQ_LOG("Constructing %u quantizers and storing inverse LUTs", QualEncoder::NR_QUANTIZERS);
        std::map<int,Quantizer> quantizers;
        unsigned int quantizerSteps = QualEncoder::QUANTIZER_STEPS_MIN;
        unsigned int quantizerIdx = QualEncoder::QUANTIZER_IDX_MIN;
        cqFile_.writeUint8((uint8_t)QualEncoder::NR_QUANTIZERS);
        for (unsigned int i = 0; i < QualEncoder::NR_QUANTIZERS; i++, quantizerIdx++, quantizerSteps++) {
            Quantizer quantizer = UniformQuantizer(qualityValueMin_, qualityValueMax_, quantizerSteps);
            quantizers.insert(std::pair<int,Quantizer>(quantizerIdx, quantizer));

            cqFile_.writeUint8(quantizerIdx);
            cqFile_.writeUint8(quantizerSteps);
            for (auto const &inverseLutEntry : quantizer.inverseLut()) {
                cqFile_.writeUint8((uint8_t)inverseLutEntry.first);
                cqFile_.writeUint8((uint8_t)inverseLutEntry.second);
            }
        }

        // Encode the QVs
        CALQ_LOG("Encoding quality values");
        QualEncoder qualEncoder(polyploidy_, qualityValueMax_, qualityValueMin_, qualityValueOffset_, quantizers);
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

    CALQ_LOG("COMPRESSION STATISTICS");
    CALQ_LOG("  Took %d ms ~= %d s ~= %d m ~= %d h", (int)diffTimeMs, (int)diffTimeS, (int)diffTimeM, (int)diffTimeH);
    CALQ_LOG("  Compressed %zu mapped + %zu unmapped = %zu record(s) in %zu block(s)", samFile_.nrMappedRecordsRead(), samFile_.nrUnmappedRecordsRead(), samFile_.nrRecordsRead(), samFile_.nrBlocksRead());
    CALQ_LOG("  Uncompressed size: %zu", uncompressedSize);
    CALQ_LOG("  Compressed size: %zu", cqFile_.nrWrittenBytes());
    CALQ_LOG("  Compression ratio: %.2f%%", (double)cqFile_.nrWrittenBytes()*100/(double)uncompressedSize);
    CALQ_LOG("  Compression factor: %.2f", (double)uncompressedSize/(double)cqFile_.nrWrittenBytes());
    CALQ_LOG("  Bits per quality value: %.4f", ((double)cqFile_.nrWrittenBytes() * 8)/(double)uncompressedSize);
    CALQ_LOG("  Speed (uncompressed size/time): %.2f MB/s", ((double)(uncompressedSize/MB))/(double)((double)diffTimeMs/1000));
}

} // namespace calq


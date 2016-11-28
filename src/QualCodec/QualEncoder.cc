/** @file QualEncoder.cc
 *  @brief This file contains the implementation of the QualEncoder class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "QualCodec/QualEncoder.h"
#include "Common/constants.h"
#include "Common/Exceptions.h"
#include "Compressors/range/range.h"
#include "Compressors/rle/rle.h"
#include <math.h>
//#include <chrono>
//#include <limits>
//#include <stdlib.h>
//#include <string.h>

static size_t rangeEncodeAndWrite(cq::CQFile &cqFile, std::string str)
{
    if (str.empty() == true) {
        throwErrorException("str is empty");
    }

    unsigned char *buffer = (unsigned char *)str.c_str();
    size_t bufferSize = str.length();

    size_t numBlocks = (size_t)ceil((double)bufferSize / (double)1*MB);
    size_t ret = cqFile.writeUint64(numBlocks);

    size_t encodedBytes = 0;
    while (encodedBytes < bufferSize) {
        unsigned int bytesToEncode = 0;
        if ((bufferSize - encodedBytes) > 1*MB) {
            bytesToEncode = 1*MB;
        } else {
            bytesToEncode = bufferSize - encodedBytes;
        }

        unsigned int compressedSize = 0;
        unsigned char *compressed = range_compress_o1(buffer+encodedBytes, (unsigned int)bytesToEncode, &compressedSize);

        if (compressedSize >= bytesToEncode) {
            ret += cqFile.writeUint8(0);
            ret += cqFile.writeUint32(bytesToEncode);
            ret += cqFile.write(buffer+encodedBytes, bytesToEncode);
        } else {
            ret += cqFile.writeUint8(1);
            ret += cqFile.writeUint32(compressedSize);
            ret += cqFile.write(compressed, compressedSize);
        }

        encodedBytes += bytesToEncode;
        free(compressed);
    }

    return ret;
}

cq::QualEncoder::QualEncoder(const unsigned int &polyploidy,
                             const int &qMin,
                             const int &qMax)
    : compressedMappedQualSize_(0)
    , compressedUnmappedQualSize_(0)
    , numMappedRecords_(0)
    , numUnmappedRecords_(0)
    , uncompressedMappedQualSize_(0)
    , uncompressedUnmappedQualSize_(0)
    , unmappedQual_("")
    , mappedQuantizerIndices_()
    , mappedQuantizerIndicesPosMin_(std::numeric_limits<uint32_t>::max())
    , mappedQuantizerIndicesPosMax_(std::numeric_limits<uint32_t>::min())
    , mappedQualIndices_()
    , pileupQueue_()
    , genotyper_(polyploidy, QUANTIZER_IDX_MIN, QUANTIZER_IDX_MAX, qMin, qMax)
    //, quantizers_()
    , samRecordQueue_()
{
    // Check arguments
    if (polyploidy == 0) {
       throwErrorException("Polyploidy must be greater than zero");
    }

//    // Init uniform quantizers
//    LOG("Initializing %u uniform quantizers", QUANTIZER_NUM);
//    for (unsigned int i = 0; i < QUANTIZER_NUM; i++) {
//        unsigned int numberOfSteps = QUANTIZER_STEP_MIN;
//        UniformQuantizer uniformQuantizer(qvMin, qvMax, numberOfSteps);
//        //uniformQuantizer.print();
//        uniformQuantizers.insert(std::pair<int,UniformQuantizer>(i, uniformQuantizer));
//    }
}

cq::QualEncoder::~QualEncoder(void)
{
    // empty
}

void cq::QualEncoder::startBlock(void)
{
    CQ_LOG("Starting block");

    blockStartTime_ = std::chrono::steady_clock::now();
    blockStopTime_ = std::chrono::steady_clock::now();

    compressedMappedQualSize_ = 0;
    compressedUnmappedQualSize_ = 0;
    numMappedRecords_ = 0;
    numUnmappedRecords_ = 0;
    uncompressedMappedQualSize_ = 0;
    uncompressedUnmappedQualSize_ = 0;

    unmappedQual_ = "";
    mappedQuantizerIndices_.clear();
    mappedQuantizerIndicesPosMin_ = std::numeric_limits<uint32_t>::max();
    mappedQuantizerIndicesPosMax_ = std::numeric_limits<uint32_t>::min();
    mappedQualIndices_.clear();

    pileupQueue_.clear();

    // nothing to be resetted for genotyper_

    // nothing to be resetted for quantizers_

    samRecordQueue_.clear();
}

void cq::QualEncoder::addUnmappedRecordToBlock(const SAMRecord &samRecord)
{
    uncompressedUnmappedQualSize_ += samRecord.qual.length();
    encodeUnmappedQual(samRecord.qual);
    numUnmappedRecords_++;
}

void cq::QualEncoder::addMappedRecordToBlock(const SAMRecord &samRecord)
{
    uncompressedMappedQualSize_ += samRecord.qual.length();

    if (numMappedRecords() == 0) {
        pileupQueue_.setPosMin(samRecord.posMin);
    }

    if (samRecord.posMax > pileupQueue_.posMax()) {
        pileupQueue_.setPosMax(samRecord.posMax);
    }

    samRecord.addToPileupQueue(pileupQueue_);
    samRecordQueue_.push_back(samRecord);

    while (pileupQueue_.posMin() < samRecord.posMin) {
        int k = genotyper_.computeQuantizerIndex(pileupQueue_.front().seq, pileupQueue_.front().qual);
        mappedQuantizerIndices_.push_back(k);
        pileupQueue_.pop_front();
    }

    while (samRecordQueue_.front().posMax < pileupQueue_.posMin()) {
        encodeMappedQual(samRecordQueue_.front());
        samRecordQueue_.pop_front();
    }

    numMappedRecords_++;
}

size_t cq::QualEncoder::finishAndWriteBlock(CQFile &cqFile)
{
    CQ_LOG("Finishing and writing block");

    // Compute all remaining quantizers
    while (pileupQueue_.posMin() <= pileupQueue_.posMax()) {
        int k = genotyper_.computeQuantizerIndex(pileupQueue_.front().seq, pileupQueue_.front().qual);
        mappedQuantizerIndices_.push_back(k);
        pileupQueue_.pop_front();
    }

    // Process all remaining records from queue
    while (samRecordQueue_.empty() == false) {
        encodeMappedQual(samRecordQueue_.front());
    }

//    // Write block number and numbers of records to output file
//    compressedMappedSizeOfBlock += cqFile.writeUint64(numBlocks);
//    compressedMappedSizeOfBlock += cqFile.writeUint64(numRecordsInBlock);
    compressedMappedQualSize_ += cqFile.writeUint64(numMappedRecords_);
//    compressedMappedSizeOfBlock += cqFile.writeUint64(numUnmappedRecordsInBlock);
//

//
//    // RLE and range encoding for quality value indices
//    //std::cout << "Quality value indices: " << qvi << std::endl;
//    unsigned char *qviBuffer = (unsigned char *)qvi.c_str();
//    size_t qviBufferSize = qvi.length();
//    if (qviBufferSize > 0) {
//        size_t qviRLESize = 0;
//        unsigned char *qviRLE = rle_encode(qviBuffer, qviBufferSize, &qviRLESize, QUANTIZER_STEP_MAX, (unsigned char)'0');
//
//        size_t numQviBlocks = (size_t)ceil((double)qviRLESize / (double)MB);
//        compressedMappedSizeOfBlock += cqFile.writeUint64(numQviBlocks);
//
//        size_t encodedQviBytes = 0;
//        while (encodedQviBytes < qviRLESize) {
//            unsigned int bytesToEncode = 0;
//            if ((qviRLESize - encodedQviBytes) > MB) {
//                bytesToEncode = MB;
//            } else {
//                bytesToEncode = qviRLESize - encodedQviBytes;
//            }
//
//            unsigned int qviRangeSize = 0;
//            unsigned char *qviRange = range_compress_o1(qviRLE+encodedQviBytes, (unsigned int)bytesToEncode, &qviRangeSize);
//
//            if (qviRangeSize >= bytesToEncode) {
//                compressedMappedSizeOfBlock += cqFile.writeUint8(0);
//                compressedMappedSizeOfBlock += cqFile.writeUint32(bytesToEncode);
//                compressedMappedSizeOfBlock += cqFile.write(qviRLE+encodedQviBytes, bytesToEncode);
//            } else {
//                compressedMappedSizeOfBlock += cqFile.writeUint8(1);
//                compressedMappedSizeOfBlock += cqFile.writeUint32(qviRangeSize);
//                compressedMappedSizeOfBlock += cqFile.write(qviRange, qviRangeSize);
//            }
//
//            encodedQviBytes += bytesToEncode;
//            free(qviRange);
//        }
//        free(qviRLE);
//    } else {
//        compressedMappedSizeOfBlock += cqFile.writeUint64(0);
//    }
//
//    // Range encoding for unmapped quality values
//    //std::cout << "Unmapped quality values: " << uqv << std::endl;
//    unsigned char *uqvBuffer = (unsigned char *)uqv.c_str();
//    size_t uqvBufferSize = uqv.length();
//    if (uqvBufferSize > 0) {
//        size_t numUqvBlocks = (size_t)ceil((double)uqvBufferSize / (double)MB);
//        compressedUnmappedSizeOfBlock += cqFile.writeUint64(numUqvBlocks);
//
//        size_t encodedUqvBytes = 0;
//        while (encodedUqvBytes < uqvBufferSize) {
//            unsigned int bytesToEncode = 0;
//            if ((uqvBufferSize - encodedUqvBytes) > MB) {
//                bytesToEncode = MB;
//            } else {
//                bytesToEncode = uqvBufferSize - encodedUqvBytes;
//            }
//
//            unsigned int uqvRangeSize = 0;
//            unsigned char *uqvRange = range_compress_o1(uqvBuffer+encodedUqvBytes, (unsigned int)bytesToEncode, &uqvRangeSize);
//
//            if (uqvRangeSize >= bytesToEncode) {
//                compressedUnmappedSizeOfBlock += cqFile.writeUint8(0);
//                compressedUnmappedSizeOfBlock += cqFile.writeUint32(bytesToEncode);
//                compressedUnmappedSizeOfBlock += cqFile.write(uqvBuffer+encodedUqvBytes, bytesToEncode);
//            } else {
//                compressedUnmappedSizeOfBlock += cqFile.writeUint8(1);
//                compressedUnmappedSizeOfBlock += cqFile.writeUint32(uqvRangeSize);
//                compressedUnmappedSizeOfBlock += cqFile.write(uqvRange, uqvRangeSize);
//            }
// 
//            encodedUqvBytes += bytesToEncode;
//            free(uqvRange);
//        }
//    } else {
//        compressedUnmappedSizeOfBlock += cqFile.writeUint64(0);
//    }
//
//    // Update sizes & counters
//    compressedSizeOfBlock = compressedMappedSizeOfBlock + compressedUnmappedSizeOfBlock;
//    compressedSize += compressedMappedSizeOfBlock + compressedUnmappedSizeOfBlock;
//    compressedMappedSize += compressedMappedSizeOfBlock;
//    compressedUnmappedSize += compressedUnmappedSizeOfBlock;
//    numBlocks++;
//
//    // Take time
    blockStopTime_ = std::chrono::steady_clock::now();
    auto diffTime = blockStopTime_ - blockStartTime_;
    auto diffTimeMs = std::chrono::duration_cast<std::chrono::milliseconds>(diffTime).count();
    auto diffTimeS = std::chrono::duration_cast<std::chrono::seconds>(diffTime).count();

//
//    LOG("Finished block %zu (contains %zu record(s))", (numBlocks-1), numRecordsInBlock);
    CQ_LOG("Took %ld ms ~= %ld s", diffTimeMs, diffTimeS);
//    LOG("Speed (uncompressed size/time): %.2f MB/s", ((double)(uncompressedSizeOfBlock/MB))/(double)((double)diffTimeMs/1000));
//    LOG("Bits per quality value: %.4f", ((double)compressedSizeOfBlock * 8)/(double)uncompressedSizeOfBlock);
//    LOG("  Mapped: %.4f", ((double)compressedMappedSizeOfBlock * 8)/(double)uncompressedMappedSizeOfBlock);
//    LOG("  Unmapped: %.4f", ((double)compressedUnmappedSizeOfBlock * 8)/(double)uncompressedUnmappedSizeOfBlock);
//
    return compressedQualSize();
}

void cq::QualEncoder::printBlockStatistics(void) const
{
    CQ_LOG("Block statistics:");
    CQ_LOG("  Uncompressed size:  %9zu", uncompressedQualSize());
    CQ_LOG("    Mapped:           %9zu", uncompressedMappedQualSize());
    CQ_LOG("    Unmapped:         %9zu", uncompressedUnmappedQualSize());
    CQ_LOG("  Compressed size:    %9zu", compressedQualSize());
    CQ_LOG("    Mapped:           %9zu", compressedMappedQualSize());
    CQ_LOG("    Unmapped:         %9zu", compressedUnmappedQualSize());
    CQ_LOG("  Records:            %9zu", numRecords());
    CQ_LOG("    Mapped:           %9zu", numMappedRecords());
    CQ_LOG("    Unmapped:         %9zu", numUnmappedRecords());
    CQ_LOG("  Compression ratio:  %12.2f%%", (double)compressedQualSize()*100/(double)uncompressedQualSize());
    CQ_LOG("    Mapped:           %12.2f%%", (double)compressedMappedQualSize()*100/(double)uncompressedMappedQualSize());
    CQ_LOG("    Unmapped:         %12.2f%%", (double)compressedUnmappedQualSize()*100/(double)uncompressedUnmappedQualSize());
    CQ_LOG("  Compression factor: %12.2f", (double)uncompressedQualSize()/(double)compressedQualSize());
    CQ_LOG("    Mapped:           %12.2f", (double)uncompressedMappedQualSize()/(double)compressedMappedQualSize());
    CQ_LOG("    Unmapped:         %12.2f", (double)uncompressedUnmappedQualSize()/(double)compressedUnmappedQualSize());
    CQ_LOG("  Bits per quality value: %.4f", ((double)compressedQualSize() * 8)/(double)uncompressedQualSize());
    CQ_LOG("    Mapped: %.4f", ((double)compressedMappedQualSize() * 8)/(double)uncompressedMappedQualSize());
    CQ_LOG("    Unmapped: %.4f", ((double)compressedUnmappedQualSize() * 8)/(double)uncompressedUnmappedQualSize());
}

size_t cq::QualEncoder::compressedMappedQualSize(void) const { return compressedMappedQualSize_; }
size_t cq::QualEncoder::compressedUnmappedQualSize(void) const { return compressedUnmappedQualSize_; }
size_t cq::QualEncoder::compressedQualSize(void) const { return compressedMappedQualSize_ + compressedUnmappedQualSize_; }
size_t cq::QualEncoder::numMappedRecords(void) const { return numMappedRecords_; }
size_t cq::QualEncoder::numUnmappedRecords(void) const { return numUnmappedRecords_; }
size_t cq::QualEncoder::numRecords(void) const { return numMappedRecords_ + numUnmappedRecords_; }
size_t cq::QualEncoder::uncompressedMappedQualSize(void) const { return uncompressedMappedQualSize_; }
size_t cq::QualEncoder::uncompressedUnmappedQualSize(void) const { return uncompressedUnmappedQualSize_; }
size_t cq::QualEncoder::uncompressedQualSize(void) const { return uncompressedMappedQualSize_ + uncompressedUnmappedQualSize_; }

void cq::QualEncoder::encodeMappedQual(const SAMRecord &samRecord)
{
//    // Iterators
//    uint32_t mrIdx = 0; // iterator for quality values in the mapped record
//    uint32_t qiIdx = 0; // iterator for quantizer indices
//
//    // Iterate through CIGAR string and code the quality values
    size_t cigarIdx = 0;
    size_t cigarLen = samRecord.cigar.length();
    size_t opLen = 0; // length of current CIGAR operation
//
    for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++) {
       if (isdigit(samRecord.cigar[cigarIdx])) {
           opLen = opLen*10 + (size_t)samRecord.cigar[cigarIdx] - (size_t)'0';
           continue;
       }
//
//        switch (mappedRecord.cigar[cigarIdx]) {
//        case 'M':
//        case '=':
//        case 'X':
//            // Encode opLen quality values with computed quantizer indices
//            for (size_t i = 0; i < opLen; i++) {
//                int qualityValue = (int)mappedRecord.qualityValues[mrIdx++];
//                int quantizerIndex = quantizerIndices[qiIdx++];
//                int qualityValueIndex = uniformQuantizers.at(quantizerIndex).valueToIndex(qualityValue);
//                qi += std::to_string(quantizerIndex);
//                qvi += std::to_string(qualityValueIndex);
//                if (quantizedPrintout == true) {
//                    int qualityValueQuantized = uniformQuantizers.at(quantizerIndex).valueToReconstructionValue(qualityValue);
//                    qualMappedFile << (char)qualityValueQuantized;
//                }
//            }
//            break;
//        case 'I':
//        case 'S':
//            // Encode opLen quality values with max quantizer index
//            for (size_t i = 0; i < opLen; i++) {
//                int qualityValue = (int)mappedRecord.qualityValues[mrIdx++];
//                int qualityValueIndex = uniformQuantizers.at(QUANTIZER_IDX_MAX).valueToIndex(qualityValue);
//                qi += std::to_string(QUANTIZER_IDX_MAX);
//                qvi += std::to_string(qualityValueIndex);
//                if (quantizedPrintout == true) {
//                    int qualityValueQuantized = uniformQuantizers.at(QUANTIZER_IDX_MAX).valueToReconstructionValue(qualityValue);
//                    qualMappedFile << (char)qualityValueQuantized;
//                }
//            }
//            break;
//        case 'D':
//        case 'N':
//            qiIdx += opLen;
//            break; // do nothing as these bases are not present
//        case 'H':
//        case 'P':
//            break; // these have been clipped
//        default:
//            LOG("CIGAR string: %s", mappedRecord.cigar.c_str());
//            throwErrorException("Bad CIGAR string");
//        }
//        opLen = 0;
    }
//
//    if (quantizedPrintout == true) {
//        qualMappedFile << std::endl;
//    }
}

void cq::QualEncoder::encodeUnmappedQual(const std::string &qual)
{
    unmappedQual_ += qual;
}


/** @file QualCodec.cc
 *  @brief This file contains the implementations of the QualEncoder and
 *         QualDecoder classes.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "QualCodec/QualCodec.h"
#include "Common/Exceptions.h"
//#include "Compressors/range/range.h"
//#include "Compressors/rle/rle.h"
//#include <chrono>
//#include <limits>
//#include <stdlib.h>
//#include <string.h>

cq::QualEncoder::QualEncoder(const unsigned int &polyploidy,
                             const int &qMin,
                             const int &qMax)
    : m_compressedMappedQualSize(0)
    , m_compressedUnmappedQualSize(0)
    , m_numMappedRecords(0)
    , m_numUnmappedRecords(0)
    , m_uncompressedMappedQualSize(0)
    , m_uncompressedUnmappedQualSize(0)
    , m_unmappedQual("")
    , m_mappedQuantizerIndices()
    , m_mappedQuantizerIndicesPosMin(std::numeric_limits<uint32_t>::max())
    , m_mappedQuantizerIndicesPosMax(std::numeric_limits<uint32_t>::min())
    , m_mappedQualIndices()
    , m_seqPileup()
    , m_qualPileup()
    , m_pileupPosMin(std::numeric_limits<uint32_t>::max())
    , m_pileupPosMax(std::numeric_limits<uint32_t>::min())
    , m_genotyper(polyploidy, QUANTIZER_IDX_MIN, QUANTIZER_IDX_MAX, qMin, qMax)
    //, m_quantizers()
    , m_samRecordQueue()
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

    m_blockStartTime = std::chrono::steady_clock::now();
    m_blockStopTime = std::chrono::steady_clock::now();

    m_compressedMappedQualSize = 0;
    m_compressedUnmappedQualSize = 0;
    m_numMappedRecords = 0;
    m_numUnmappedRecords = 0;
    m_uncompressedMappedQualSize = 0;
    m_uncompressedUnmappedQualSize = 0;

    m_unmappedQual = "";
    m_mappedQuantizerIndices.clear();
    m_mappedQuantizerIndicesPosMin = std::numeric_limits<uint32_t>::max();
    m_mappedQuantizerIndicesPosMax = std::numeric_limits<uint32_t>::min();
    m_mappedQualIndices.clear();

    m_seqPileup.clear();
    m_qualPileup.clear();
    m_pileupPosMin = std::numeric_limits<uint32_t>::max();
    m_pileupPosMax = std::numeric_limits<uint32_t>::min();

    // nothing to be resetted for m_genotyper

    // nothing to be resetted for m_quantizers

    m_samRecordQueue.clear();
}

void cq::QualEncoder::addUnmappedRecordToBlock(const SAMRecord &samRecord)
{
    m_uncompressedUnmappedQualSize += samRecord.qual.length();
    encodeUnmappedQual(samRecord.qual);
    m_numUnmappedRecords++;
}

void cq::QualEncoder::addMappedRecordToBlock(const SAMRecord &samRecord)
{
    m_uncompressedMappedQualSize += samRecord.qual.length();

    if (numMappedRecords() == 0) {
        m_pileupPosMin = samRecord.posMin;
        m_mappedQuantizerIndicesPosMin = m_pileupPosMin;
    }

    if (samRecord.posMax > m_pileupPosMax) {
        m_pileupPosMax = samRecord.posMax;
        size_t pileupLength = m_pileupPosMax - m_pileupPosMin + 1;
        m_seqPileup.resize(pileupLength);
        m_qualPileup.resize(pileupLength);
    }

    samRecord.extractPileup(m_pileupPosMin, m_pileupPosMax, m_seqPileup, m_qualPileup);
    m_samRecordQueue.push_back(samRecord);

    while (m_pileupPosMin < samRecord.posMin) {
        int k = m_genotyper.computeQuantizerIndex(m_seqPileup[0], m_qualPileup[0]);

        m_mappedQuantizerIndices.push_back(k);
        m_mappedQuantizerIndicesPosMax = m_pileupPosMin;
        m_seqPileup.pop_front();
        m_qualPileup.pop_front();
        m_pileupPosMin++;
    }

    while (m_samRecordQueue.front().posMax < m_pileupPosMin) {
        encodeMappedQual(m_samRecordQueue.front());
        m_samRecordQueue.pop_front();
    }

    m_numMappedRecords++;
}

size_t cq::QualEncoder::finishAndWriteBlock(CQFile &cqFile)
{
    CQ_LOG("Finishing and writing block");

    // Compute all remaining quantizers
    while (m_pileupPosMin <= m_pileupPosMax) {
        int k = m_genotyper.computeQuantizerIndex(m_seqPileup[0], m_qualPileup[0]);

        m_mappedQuantizerIndices.push_back(k);
        m_mappedQuantizerIndicesPosMax = m_pileupPosMin;
        m_seqPileup.pop_front();
        m_qualPileup.pop_front();
        m_pileupPosMin++;
    }

    // Process all remaining records from queue
    while (m_samRecordQueue.empty() == false) {
        encodeMappedQual(m_samRecordQueue.front());
        m_samRecordQueue.pop_front();
    }

//    // Write block number and numbers of records to output file
//    compressedMappedSizeOfBlock += cqFile.writeUint64(numBlocks);
//    compressedMappedSizeOfBlock += cqFile.writeUint64(numRecordsInBlock);
    m_compressedMappedQualSize += cqFile.writeUint64(m_numMappedRecords);
//    compressedMappedSizeOfBlock += cqFile.writeUint64(numUnmappedRecordsInBlock);
//
//    // RLE and range encoding for quantizer indices
//    //std::cout << "Quantizer indices: " << qi << std::endl;
//    unsigned char *qiBuffer = (unsigned char *)qi.c_str();
//    size_t qiBufferSize = qi.length();
//    if (qiBufferSize > 0) {
//        size_t qiRLESize = 0;
//        unsigned char *qiRLE = rle_encode(qiBuffer, qiBufferSize, &qiRLESize, QUANTIZER_NUM, (unsigned char)'0');
//
//        size_t numQiBlocks = (size_t)ceil((double)qiRLESize / (double)MB);
//        compressedMappedSizeOfBlock += cqFile.writeUint64(numQiBlocks);
//
//        size_t encodedQiBytes = 0;
//        while (encodedQiBytes < qiRLESize) {
//            unsigned int bytesToEncode = 0;
//            if ((qiRLESize - encodedQiBytes) > MB) {
//                bytesToEncode = MB;
//            } else {
//                bytesToEncode = qiRLESize - encodedQiBytes;
//            }
//
//            unsigned int qiRangeSize = 0;
//            unsigned char *qiRange = range_compress_o1(qiRLE+encodedQiBytes, (unsigned int)bytesToEncode, &qiRangeSize);
//
//            if (qiRangeSize >= bytesToEncode) {
//                compressedMappedSizeOfBlock += cqFile.writeUint8(0);
//                compressedMappedSizeOfBlock += cqFile.writeUint32(bytesToEncode);
//                compressedMappedSizeOfBlock += cqFile.write(qiRLE+encodedQiBytes, bytesToEncode);
//            } else {
//                compressedMappedSizeOfBlock += cqFile.writeUint8(1);
//                compressedMappedSizeOfBlock += cqFile.writeUint32(qiRangeSize);
//                compressedMappedSizeOfBlock += cqFile.write(qiRange, qiRangeSize);
//            }
//
//            encodedQiBytes += bytesToEncode;
//            free(qiRange);
//        }
//        free(qiRLE);
//    } else {
//        compressedMappedSizeOfBlock += cqFile.writeUint64(0);
//    }
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
    m_blockStopTime = std::chrono::steady_clock::now();
    auto diffTime = m_blockStopTime - m_blockStartTime;
    auto diffTimeMs = std::chrono::duration_cast<std::chrono::milliseconds>(diffTime).count();
    auto diffTimeS = std::chrono::duration_cast<std::chrono::seconds>(diffTime).count();
//
//    // Print statistics for this block
//    if (verbose == true) {
//        LOG("Block statistics:");
//        LOG("  Block number:       %9zu", (numBlocks-1));
//        LOG("  Uncompressed size:  %9zu", uncompressedSizeOfBlock);
//        LOG("    Mapped:           %9zu", uncompressedMappedSizeOfBlock);
//        LOG("    Unmapped:         %9zu", uncompressedUnmappedSizeOfBlock);
//        LOG("  Compressed size:    %9zu", compressedSizeOfBlock);
//        LOG("    Mapped:           %9zu", compressedMappedSizeOfBlock);
//        LOG("    Unmapped:         %9zu", compressedUnmappedSizeOfBlock);
//        LOG("  Records:            %9zu", numRecordsInBlock);
//        LOG("    Mapped:           %9zu", numMappedRecordsInBlock);
//        LOG("    Unmapped:         %9zu", numUnmappedRecordsInBlock);
//        LOG("  Compression ratio:  %12.2f%%", (double)compressedSizeOfBlock*100/(double)uncompressedSizeOfBlock);
//        LOG("    Mapped:           %12.2f%%", (double)compressedMappedSizeOfBlock*100/(double)uncompressedMappedSizeOfBlock);
//        LOG("    Unmapped:         %12.2f%%", (double)compressedUnmappedSizeOfBlock*100/(double)uncompressedUnmappedSizeOfBlock);
//        LOG("  Compression factor: %12.2f", (double)uncompressedSizeOfBlock/(double)compressedSizeOfBlock);
//        LOG("    Mapped:           %12.2f", (double)uncompressedMappedSizeOfBlock/(double)compressedMappedSizeOfBlock);
//        LOG("    Unmapped:         %12.2f", (double)uncompressedUnmappedSizeOfBlock/(double)compressedUnmappedSizeOfBlock);
//    }
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

size_t cq::QualEncoder::compressedMappedQualSize(void) const { return m_compressedMappedQualSize; }
size_t cq::QualEncoder::compressedUnmappedQualSize(void) const { return m_compressedUnmappedQualSize; }
size_t cq::QualEncoder::compressedQualSize(void) const { return m_compressedMappedQualSize + m_compressedUnmappedQualSize; }
size_t cq::QualEncoder::numMappedRecords(void) const { return m_numMappedRecords; }
size_t cq::QualEncoder::numUnmappedRecords(void) const { return m_numUnmappedRecords; }
size_t cq::QualEncoder::numRecords(void) const { return m_numMappedRecords + m_numUnmappedRecords; }
size_t cq::QualEncoder::uncompressedMappedQualSize(void) const { return m_uncompressedMappedQualSize; }
size_t cq::QualEncoder::uncompressedUnmappedQualSize(void) const { return m_uncompressedUnmappedQualSize; }
size_t cq::QualEncoder::uncompressedQualSize(void) const { return m_uncompressedMappedQualSize + m_uncompressedUnmappedQualSize; }

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
    m_unmappedQual += qual;
}

cq::QualDecoder::QualDecoder(void)
{
   // empty
}

cq::QualDecoder::~QualDecoder(void)
{
   //empty
}

//void QualDecoder::printStats(void) const
//{
//    if (verbose == true) {
//        LOG("QualDecoder statistics:");
//        LOG("  Blocks:             %9zu", numBlocks);
//        LOG("  Uncompressed size:  %9zu", uncompressedSize);
//        LOG("  Compressed size:    %9zu", compressedSize);
//        LOG("  Records:            %9zu", numRecords);
//        LOG("    Mapped:           %9zu", numMappedRecords);
//        LOG("    Unmapped:         %9zu", numUnmappedRecords);
//        LOG("  Compression ratio:  %12.2f%%", (double)compressedSize*100/(double)uncompressedSize);
//        LOG("  Compression factor: %12.2f", (double)uncompressedSize/(double)compressedSize);
//    }
//
//    LOG("Decoded %zu record(s) in %zu block(s)", numRecords, numBlocks);
//    LOG("Bits per quality value: %.4f", ((double)compressedSize * 8)/(double)uncompressedSize);
//}
//
//size_t QualDecoder::decodeBlock(void)
//{
//    size_t ret = 0;
//
//    LOG("Decoding block %zu", numBlocks);
//    auto startTime = std::chrono::steady_clock::now();
//
//    // Reset all variables in block scope
//    uncompressedSizeOfBlock = 0;
//    compressedSizeOfBlock = 0;
//    numRecordsInBlock = 0;
//    numMappedRecordsInBlock = 0;
//    numUnmappedRecordsInBlock = 0;
//
//    // Read and check block number
//    uint64_t blockNumber = 0;
//    ret += cqFile.readUint64(&blockNumber);
//    if (numBlocks != blockNumber) {
//        LOG("Block number: %lu", blockNumber);
//        throwErrorException("Corrupted block number in CQ file");
//    }
//
//    // Read record numbers for this block
//    ret += cqFile.readUint64((uint64_t *)&numRecordsInBlock);
//    ret += cqFile.readUint64((uint64_t *)&numMappedRecordsInBlock);
//    ret += cqFile.readUint64((uint64_t *)&numUnmappedRecordsInBlock);
//
//    std::string qiRLEString("");
//    uint64_t numQiBlocks = 0;
//    ret += cqFile.readUint64(&numQiBlocks);
//    for (size_t i = 0; i < numQiBlocks; i++) {
//        uint8_t qiRangeFlag = 0;
//        ret += cqFile.readUint8(&qiRangeFlag);
//        if (qiRangeFlag == 1) {
//            uint32_t qiRangeSize = 0;
//            ret += cqFile.readUint32(&qiRangeSize);
//            unsigned char *qiRange = (unsigned char *)malloc(qiRangeSize);
//            ret += cqFile.read(qiRange, qiRangeSize);
//            unsigned int qiRLESize = 0;
//            unsigned char *qiRLE= range_decompress_o1(qiRange, &qiRLESize);
//            free(qiRange);
//            std::string qiRLEStringTmp((char *)qiRLE, qiRLESize);
//            free(qiRLE);
//            qiRLEString += qiRLEStringTmp;
//        } else { // qiRangeFlag == 0
//            unsigned int qiRLESize = 0;
//            ret += cqFile.readUint32(&qiRLESize);
//            unsigned char *qiRLE = (unsigned char *)malloc(qiRLESize);
//            ret += cqFile.read(qiRLE, qiRLESize);
//            std::string qiRLEStringTmp((char *)qiRLE, qiRLESize);
//            free(qiRLE);
//            qiRLEString += qiRLEStringTmp;
//        }
//
//        size_t qiSize = 0;
//        unsigned char *qi = rle_decode((unsigned char *)qiRLEString.c_str(), qiRLEString.length(), &qiSize, QUANTIZER_NUM, (unsigned int)'0');
//        //std::string qiString((char *)qi, qiSize);
//        //std::cout << "Decoded quantizer indices: " << qiString << std::endl;
//        free(qi);
//    }
//
//    std::string qviRLEString("");
//    uint64_t numQviBlocks = 0;
//    ret += cqFile.readUint64(&numQviBlocks);
//    for (size_t i = 0; i < numQviBlocks; i++) {
//        uint8_t qviRangeFlag = 0;
//        ret += cqFile.readUint8(&qviRangeFlag);
//        if (qviRangeFlag == 1) {
//            uint32_t qviRangeSize = 0;
//            ret += cqFile.readUint32(&qviRangeSize);
//            unsigned char *qviRange = (unsigned char *)malloc(qviRangeSize);
//            ret += cqFile.read(qviRange, qviRangeSize);
//            unsigned int qviRLESize = 0;
//            unsigned char *qviRLE= range_decompress_o1(qviRange, &qviRLESize);
//            free(qviRange);
//            std::string qviRLEStringTmp((char *)qviRLE, qviRLESize);
//            free(qviRLE);
//            qviRLEString += qviRLEStringTmp;
//        } else { // qviRangeFlag == 0
//            unsigned int qviRLESize = 0;
//            ret += cqFile.readUint32(&qviRLESize);
//            unsigned char *qviRLE = (unsigned char *)malloc(qviRLESize);
//            ret += cqFile.read(qviRLE, qviRLESize);
//            std::string qviRLEStringTmp((char *)qviRLE, qviRLESize);
//            free(qviRLE);
//            qviRLEString += qviRLEStringTmp;
//        }
//
//        size_t qviSize = 0;
//        unsigned char *qvi = rle_decode((unsigned char *)qviRLEString.c_str(), qviRLEString.length(), &qviSize, QUANTIZER_STEP_MAX, (unsigned int)'0');
//        //std::string qviString((char *)qvi, qviSize);
//        //std::cout << "Decoded quality value indices: " << qviString << std::endl;
//        free(qvi);
//    }
//
//    //std::string uqvString("");
//    uint64_t numUqvBlocks = 0;
//    ret += cqFile.readUint64(&numUqvBlocks);
//    for (size_t i = 0; i < numUqvBlocks; i++) {
//        uint8_t uqvRangeFlag = 0;
//        ret += cqFile.readUint8(&uqvRangeFlag);
//        if (uqvRangeFlag == 1) {
//            uint32_t uqvRangeSize = 0;
//            ret += cqFile.readUint32(&uqvRangeSize);
//            unsigned char *uqvRange = (unsigned char *)malloc(uqvRangeSize);
//            ret += cqFile.read(uqvRange, uqvRangeSize);
//            unsigned int uqvSize = 0;
//            unsigned char *uqv= range_decompress_o1(uqvRange, &uqvSize);
//            free(uqvRange);
//            //std::string uqvStringTmp((char *)uqv, uqvSize);
//            free(uqv);
//            //uqvString += uqvStringTmp;
//        } else { // uqvRangeFlag == 0
//            unsigned int uqvSize = 0;
//            ret += cqFile.readUint32(&uqvSize);
//            unsigned char *uqv = (unsigned char *)malloc(uqvSize);
//            ret += cqFile.read(uqv, uqvSize);
//            //std::string uqvStringTmp((char *)uqv, uqvSize);
//            free(uqv);
//            //uqvString += uqvStringTmp;
//        }
//        //std::cout << "Decoded unmapped quality values: " << uqvString << std::endl;
//    }
//
//    // Update sizes & counters
//    compressedSizeOfBlock += ret;
//    compressedSize += compressedSizeOfBlock;
//    numBlocks++;
//    numRecords += numRecordsInBlock;
//    numMappedRecords += numMappedRecordsInBlock;
//    numUnmappedRecords += numUnmappedRecordsInBlock;
//
//    // Take time
//    auto stopTime = std::chrono::steady_clock::now();
//    auto diffTime = stopTime - startTime;
//    auto diffTimeMs = std::chrono::duration_cast<std::chrono::milliseconds>(diffTime).count();
//    auto diffTimeS = std::chrono::duration_cast<std::chrono::seconds>(diffTime).count();
//
//    // Print statistics for this block
//    if (verbose == true) {
//        LOG("Block statistics:");
//        LOG("  Block number:       %9zu", blockNumber);
//        LOG("  Uncompressed size:  %9zu", uncompressedSizeOfBlock);
//        LOG("  Compressed size:    %9zu", compressedSizeOfBlock);
//        LOG("  Records:            %9zu", numRecordsInBlock);
//        LOG("    Mapped:           %9zu", numMappedRecordsInBlock);
//        LOG("    Unmapped:         %9zu", numUnmappedRecordsInBlock);
//        LOG("  Compression ratio:  %12.2f%%", (double)compressedSizeOfBlock*100/(double)uncompressedSizeOfBlock);
//        LOG("  Compression factor: %12.2f", (double)uncompressedSizeOfBlock/(double)compressedSizeOfBlock);
//    }
//
//    LOG("Finished block %zu (contains %zu records)", blockNumber, numRecordsInBlock);
//    LOG("Took %ld ms ~= %ld s", diffTimeMs, diffTimeS);
//    LOG("Speed (uncompressed size/time): %.2f MB/s", ((double)(uncompressedSizeOfBlock/MB))/(double)((double)diffTimeMs/1000));
//    LOG("Bits per quality value: %.4f", ((double)compressedSizeOfBlock * 8)/(double)uncompressedSizeOfBlock);
//
//    return ret;
//}


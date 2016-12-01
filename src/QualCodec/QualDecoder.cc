/** @file QualDecoder.cc
 *  @brief This file contains the implementation of the QualDecoder class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */
#include "QualCodec/QualDecoder.h"

#include "Common/Exceptions.h"
#include "Common/log.h"

namespace calq {

QualDecoder::QualDecoder(void)
    : unmappedQualityValues_("")
    , mappedQuantizerIndices_("")
    , mappedQualityValueIndices_("") {}

QualDecoder::~QualDecoder(void) {}

size_t QualDecoder::readBlock(CQFile &cqFile)
{
    CALQ_LOG("Reading block");

    size_t ret = 0;

    CALQ_LOG("Reading unmapped quality values");
    uint8_t flag = 0;
    cqFile.readUint8(&flag);
    if (flag)
        ret += cqFile.readBuffer(&unmappedQualityValues_);
    std::cout << unmappedQualityValues_ << std::endl;

    CALQ_LOG("Reading mapped quantizer indices");
    ret += cqFile.readBuffer(&mappedQuantizerIndices_);
    std::cout << mappedQuantizerIndices_ << std::endl;

    CALQ_LOG("Reading mapped quality value indices");
    ret += cqFile.readBuffer(&mappedQualityValueIndices_);

    return ret;
}

void QualDecoder::reset(void)
{
    unmappedQualityValues_ = "";
    mappedQuantizerIndices_ = "";
    mappedQualityValueIndices_ = "";
}




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

} // namespace calq


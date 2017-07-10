/** @file QualEncoder.cc
 *  @brief This file contains the implementation of the QualEncoder class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#include "QualCodec/QualEncoder.h"

#include <math.h>

#include <fstream>
#include <map>
#include <string>

#include "Common/constants.h"
#include "Common/Exceptions.h"
#include "Common/log.h"
#include "config.h"

namespace calq {

QualEncoder::QualEncoder(const int &polyploidy,
                         const int &qualityValueMax,
                         const int &qualityValueMin,
                         const int &qualityValueOffset,
                         const std::map<int, Quantizer> &quantizers)
    : startTime_(std::chrono::steady_clock::now()),
      stopTime_(std::chrono::steady_clock::now()),

      compressedMappedQualSize_(0),
      compressedUnmappedQualSize_(0),
      nrMappedRecords_(0),
      nrUnmappedRecords_(0),
      uncompressedMappedQualSize_(0),
      uncompressedUnmappedQualSize_(0),
#if MPEG_CE5_DESCRIPTOR_STREAMS_COMPRESSION_EXTENSION
      compressedQviSize_(0),
      compressedQvciSize_(0),
      uncompressedQvciSize_(0),
      uncompressedQviSize_(0),
#endif

      qualityValueOffset_(qualityValueOffset),
      posOffset_(0),

      unmappedQualityValues_(""),
      mappedQuantizerIndices_(),
      mappedQualityValueIndices_(),
#if MPEG_CE5_DESCRIPTOR_STREAMS_COMPRESSION_EXTENSION
      qvi_(),
#endif

      samPileupDeque_(),

      genotyper_(polyploidy, qualityValueMin, qualityValueMax, qualityValueOffset, quantizers.size()),

      quantizers_(quantizers),

      samRecordDeque_()
{
    if (polyploidy < 1) {
        throwErrorException("polyploidy must be greater than zero");
    }
    if (qualityValueMax < 0) {
        throwErrorException("qualityValueMax must be zero or greater");
    }
    if (qualityValueMin < 0) {
        throwErrorException("qualityValueMin must be zero or greater");
    }
    if (qualityValueOffset < 1) {
        throwErrorException("qualityValueOffset must be greater than zero");
    }
    if (quantizers.empty() == true) {
        throwErrorException("quantizers is empty");
    }

#if MPEG_CE5_DESCRIPTOR_STREAMS_COMPRESSION_EXTENSION
    for (int i = QUANTIZER_IDX_MIN; i <= QUANTIZER_IDX_MAX; ++i) {
        qvi_.push_back(std::deque<int>());
    }
//     CALQ_LOG("Created a vector of %d double-ended queues for temporary QVI storage", qvi_.size());
#endif
}

QualEncoder::~QualEncoder(void)
{
    // empty
}

void QualEncoder::addUnmappedRecordToBlock(const SAMRecord &samRecord)
{
    encodeUnmappedQual(samRecord.qual);

    uncompressedUnmappedQualSize_ += samRecord.qual.length();
    nrUnmappedRecords_++;
}

void QualEncoder::addMappedRecordToBlock(const SAMRecord &samRecord)
{
    if (nrMappedRecords() == 0) {
        posOffset_ = samRecord.posMin;
        samPileupDeque_.setPosMin(samRecord.posMin);
        samPileupDeque_.setPosMax(samRecord.posMax);
    }

    if (samRecord.posMax > samPileupDeque_.posMax()) {
        samPileupDeque_.setPosMax(samRecord.posMax);
    }

    samRecord.addToPileupQueue(&samPileupDeque_);
    samRecordDeque_.push_back(samRecord);

    while (samPileupDeque_.posMin() < samRecord.posMin) {
        int k = genotyper_.computeQuantizerIndex(samPileupDeque_.front().seq, samPileupDeque_.front().qual);
        mappedQuantizerIndices_.push_back(k);
        samPileupDeque_.pop_front();
    }

    while (samRecordDeque_.front().posMax < samPileupDeque_.posMin()) {
        encodeMappedQual(samRecordDeque_.front());
        samRecordDeque_.pop_front();
    }

    uncompressedMappedQualSize_ += samRecord.qual.length();
    nrMappedRecords_++;
}

void QualEncoder::finishBlock(void)
{
    // Compute all remaining quantizers
    while (samPileupDeque_.empty() == false) {
        int k = genotyper_.computeQuantizerIndex(samPileupDeque_.front().seq, samPileupDeque_.front().qual);
        mappedQuantizerIndices_.push_back(k);
        samPileupDeque_.pop_front();
    }

    // Process all remaining records from queue
    while (samRecordDeque_.empty() == false) {
        encodeMappedQual(samRecordDeque_.front());
        samRecordDeque_.pop_front();
    }
}

size_t QualEncoder::writeBlock(CQFile *cqFile)
{
//     CALQ_LOG("Writing block");

    compressedMappedQualSize_ = 0;
    compressedUnmappedQualSize_ = 0;

//     CALQ_LOG("  Writing pos offset and quality value offset");
    compressedMappedQualSize_ += cqFile->writeUint32(posOffset_);
    compressedMappedQualSize_ += cqFile->writeUint32((uint32_t)qualityValueOffset_);

//     CALQ_LOG("  Writing unmapped quality values");
//     std::cout << "uqv: " << unmappedQualityValues_ << std::endl;
    unsigned char *uqv = (unsigned char *)unmappedQualityValues_.c_str();
    size_t uqvSize = unmappedQualityValues_.length();
    uint8_t uqvFlags = 0;
    if (uqvSize > 0) {
        uqvFlags = 1;
        compressedUnmappedQualSize_ += cqFile->writeUint8(uqvFlags);
        compressedUnmappedQualSize_ += cqFile->writeQualBlock(uqv, uqvSize);
    } else {
        uqvFlags = 0;
        compressedUnmappedQualSize_ += cqFile->writeUint8(uqvFlags);
//         CALQ_LOG("  No unmapped quality values in this block");
    }

//     CALQ_LOG("  Writing mapped quantizer indices");
    std::string tmp("");
    for (auto const &mappedQuantizerIndex : mappedQuantizerIndices_) {
        tmp += std::to_string(mappedQuantizerIndex+1);
    }

//     std::cout << "mqi: " << tmp << std::endl;
    unsigned char *mqi = (unsigned char *)tmp.c_str();
    size_t mqiSize = tmp.length();
    uint8_t mqiFlags = 0;
    if (mqiSize > 0) {
        mqiFlags = 1;
        compressedMappedQualSize_ += cqFile->writeUint8(mqiFlags);
        compressedMappedQualSize_ += cqFile->writeQualBlock(mqi, mqiSize);
    } else {
        mqiFlags = 0;
        compressedMappedQualSize_ += cqFile->writeUint8(mqiFlags);
//         CALQ_LOG("  No mapped quantizer indices in this block");
    }

//     CALQ_LOG("  Writing mapped quality value indices");
    tmp = "";
    for (auto const &mappedQualityValueIndex : mappedQualityValueIndices_) {
        tmp += std::to_string(mappedQualityValueIndex);
    }
//     std::cout << "mqvi: " << tmp << std::endl;
    unsigned char *mqvi = (unsigned char *)tmp.c_str();
    size_t mqviSize = tmp.length();
    uint8_t mqviFlags = 0;
    if (mqviSize > 0) {
        mqviFlags = 1;
        compressedMappedQualSize_ += cqFile->writeUint8(mqviFlags);
        compressedMappedQualSize_ += cqFile->writeQualBlock(mqvi, mqviSize);
    } else {
        mqviFlags = 0;
        compressedMappedQualSize_ += cqFile->writeUint8(mqviFlags);
//         CALQ_LOG("  No mapped quality value in this block");
    }

    return compressedQualSize();
}

#if MPEG_CE5_DESCRIPTOR_STREAMS_OUTPUT
void QualEncoder::writeMPEGBlock(File *mpegQVCIFile, File *mpegQVIFile)
{
    std::string tmp("");
    for (auto const &mappedQuantizerIndex : mappedQuantizerIndices_) {
        tmp += std::to_string(mappedQuantizerIndex + 1);
    }
//          std::cout << "QVCodebookIdentifier: " << tmp << std::endl;

    unsigned char *QVCodebookIdentifier = (unsigned char *)tmp.c_str();
    size_t QVCodebookIdentifierBufferSize = tmp.length();
    uint32_t PositionOffset = posOffset_;

    mpegQVCIFile->writeUint32(PositionOffset);
    mpegQVCIFile->writeUint64(QVCodebookIdentifierBufferSize);
    mpegQVCIFile->write(QVCodebookIdentifier, QVCodebookIdentifierBufferSize);

    tmp = "";
    for (auto const &mappedQualityValueIndex : mappedQualityValueIndices_) {
        tmp += std::to_string(mappedQualityValueIndex);
    }
//          std::cout << "QVIndex: " << tmp << std::endl;

    unsigned char *QVIndex = (unsigned char *)tmp.c_str();
    size_t QVIndexBufferSize = tmp.length();

    mpegQVIFile->writeUint64(QVIndexBufferSize);
    mpegQVIFile->write(QVIndex, QVIndexBufferSize);

}
#endif

#if MPEG_CE5_DESCRIPTOR_STREAMS_COMPRESSION_EXTENSION
void QualEncoder::writeContextAdaptiveMPEGBlock(CQFile *qvciFile, CQFile *qviFile)
{
    compressedQvciSize_+= qvciFile->writeUint32(posOffset_);

    std::string tmp("");
    for (auto const &mappedQuantizerIndex : mappedQuantizerIndices_) {
        tmp += std::to_string(mappedQuantizerIndex+1);
    }

    unsigned char *qvci = (unsigned char *)tmp.c_str();
    size_t qvciSize = tmp.length();
    uint8_t qvciFlags = 0;
    if (qvciSize > 0) {
        qvciFlags = 1;
        uncompressedQvciSize_ += qvciSize;
        compressedQvciSize_ += qvciFile->writeUint8(qvciFlags);
        compressedQvciSize_ += qvciFile->writeQualBlock(qvci, qvciSize);
    } else {
        qvciFlags = 0;
        compressedQvciSize_ += qvciFile->writeUint8(qvciFlags);
    }

    tmp = "";
    int idx = 0;
    for (std::deque<int> const &qviStream : qvi_) {
        tmp = "";
        for (auto const &i : qviStream) {
            tmp += std::to_string(i);
        }
        unsigned char *qvi = (unsigned char *)tmp.c_str();
        size_t qviSize = tmp.length();
        uint8_t qviFlags = 0;
        if (qviSize > 0) {
            qviFlags = 1;
            compressedQviSize_ += qviFile->writeUint8(qviFlags);
            compressedQviSize_ += qviFile->writeQualBlock(qvi, qviSize);
        } else {
            qviFlags = 0;
            compressedQviSize_ += qviFile->writeUint8(qviFlags);
        }
    }
}
#endif

size_t QualEncoder::compressedMappedQualSize(void) const { return compressedMappedQualSize_; }
size_t QualEncoder::compressedUnmappedQualSize(void) const { return compressedUnmappedQualSize_; }
size_t QualEncoder::compressedQualSize(void) const { return compressedMappedQualSize_ + compressedUnmappedQualSize_; }
size_t QualEncoder::nrMappedRecords(void) const { return nrMappedRecords_; }
size_t QualEncoder::nrUnmappedRecords(void) const { return nrUnmappedRecords_; }
size_t QualEncoder::nrRecords(void) const { return nrMappedRecords_ + nrUnmappedRecords_; }
size_t QualEncoder::uncompressedMappedQualSize(void) const { return uncompressedMappedQualSize_; }
size_t QualEncoder::uncompressedUnmappedQualSize(void) const { return uncompressedUnmappedQualSize_; }
size_t QualEncoder::uncompressedQualSize(void) const { return uncompressedMappedQualSize_ + uncompressedUnmappedQualSize_; }
#if MPEG_CE5_DESCRIPTOR_STREAMS_COMPRESSION_EXTENSION
size_t QualEncoder::compressedQviSize(void) const { return compressedQviSize_; }
size_t QualEncoder::compressedQvciSize(void) const { return compressedQvciSize_; }
size_t QualEncoder::uncompressedQvciSize(void) const { return uncompressedQvciSize_; }
size_t QualEncoder::uncompressedQviSize(void) const { return uncompressedQviSize_; }
#endif

void QualEncoder::encodeMappedQual(const SAMRecord &samRecord)
{
//     printf("Encoding SAM record");
//     samRecord.printShort();

    size_t cigarIdx = 0;
    size_t cigarLen = samRecord.cigar.length();
    size_t opLen = 0; // length of current CIGAR operation
    size_t qualIdx = 0;
    size_t quantizerIndicesIdx = samRecord.posMin - posOffset_;

    for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++) {
       if (isdigit(samRecord.cigar[cigarIdx])) {
           opLen = opLen*10 + (size_t)samRecord.cigar[cigarIdx] - (size_t)'0';
           continue;
       }

       switch (samRecord.cigar[cigarIdx]) {
       case 'M':
       case '=':
       case 'X':
           // Encode opLen quality values with computed quantizer indices
           for (size_t i = 0; i < opLen; i++) {
               int q = (int)samRecord.qual[qualIdx++] - qualityValueOffset_;
               int quantizerIndex = mappedQuantizerIndices_[quantizerIndicesIdx++];
               int qualityValueIndex = quantizers_.at(quantizerIndex).valueToIndex(q);
               mappedQualityValueIndices_.push_back(qualityValueIndex);
#if MPEG_CE5_DESCRIPTOR_STREAMS_COMPRESSION_EXTENSION
               int qvci = quantizerIndex;
               int qvi = qualityValueIndex;
               qvi_.at(qvci).push_back(qvi);
               uncompressedQviSize_++;
#endif
//                printf("Encoding %c with k=%d -> %d\n", (char)(q+qualityValueOffset_), quantizerIndex, qualityValueIndex);
           }
           break;
       case 'I':
       case 'S':
           // Encode opLen quality values with max quantizer index
           for (size_t i = 0; i < opLen; i++) {
               int q = (int)samRecord.qual[qualIdx++] - qualityValueOffset_;
               int quantizerIndex = quantizers_.size() - 1;
               int qualityValueIndex = quantizers_.at(quantizerIndex).valueToIndex(q);
               mappedQualityValueIndices_.push_back(qualityValueIndex);
#if MPEG_CE5_DESCRIPTOR_STREAMS_COMPRESSION_EXTENSION
               int qvci = quantizerIndex;
               int qvi = qualityValueIndex;
               qvi_.at(qvci).push_back(qvi);
               uncompressedQviSize_++;
#endif
//                printf("Encoding %c with kmax=%d -> %d\n", (char)(q+qualityValueOffset_), quantizerIndex, qualityValueIndex);
           }
           break;
       case 'D':
       case 'N':
           quantizerIndicesIdx += opLen;
           break; // do nothing as these bases are not present
       case 'H':
       case 'P':
           break; // these have been clipped
       default:
           throwErrorException("Bad CIGAR string");
       }
       opLen = 0;
    }
}

void QualEncoder::encodeUnmappedQual(const std::string &qual)
{
    unmappedQualityValues_ += qual;
}

} // namespace calq


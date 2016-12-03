/** @file QualEncoder.cc
 *  @brief This file contains the implementation of the QualEncoder class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#include "QualCodec/QualEncoder.h"

#include <math.h>

#include "Common/constants.h"
#include "Common/Exceptions.h"
#include "Common/log.h"
#include "Compressors/rle/rle.h"

namespace calq {

QualEncoder::QualEncoder(const int &polyploidy,
                         const int &qualityValueMax,
                         const int &qualityValueMin,
                         const int &qualityValueOffset,
                         const std::map<int,Quantizer> &quantizers)
    : startTime_(std::chrono::steady_clock::now()),
      stopTime_(std::chrono::steady_clock::now()),

      compressedMappedQualSize_(0),
      compressedUnmappedQualSize_(0),
      nrMappedRecords_(0),
      nrUnmappedRecords_(0),
      uncompressedMappedQualSize_(0),
      uncompressedUnmappedQualSize_(0),

      qualityValueOffset_(qualityValueOffset),
      posOffset_(0),

      unmappedQual_(""),
      mappedQuantizerIndices_(),
      mappedQualIndices_(),

      samPileupDeque_(),

      genotyper_(polyploidy, qualityValueMin, qualityValueMax, quantizers.size()),

      quantizers_(quantizers),

      samRecordDeque_()
{
    if (polyploidy == 0) {
       throwErrorException("Polyploidy must be greater than zero");
    }
    if (qualityValueMin > qualityValueMax) {
        throwErrorException("qMin is greater than qMax");
    }
    if (quantizers.empty() == true) {
        throwErrorException("quantizers is empty");
    }
}

QualEncoder::~QualEncoder(void)
{
    // empty
}

void QualEncoder::addUnmappedRecordToBlock(const SAMRecord &samRecord)
{
    uncompressedUnmappedQualSize_ += samRecord.qual.length();
    encodeUnmappedQual(samRecord.qual);
    nrUnmappedRecords_++;
}

void QualEncoder::addMappedRecordToBlock(const SAMRecord &samRecord)
{
    uncompressedMappedQualSize_ += samRecord.qual.length();

    if (nrMappedRecords() == 0) {
        posOffset_ = samRecord.posMin;
        samPileupDeque_.setPosMin(samRecord.posMin);
        samPileupDeque_.setPosMax(samRecord.posMax);
    }

    if (samRecord.posMax > samPileupDeque_.posMax()) {
        samPileupDeque_.setPosMax(samRecord.posMax);
    }

    samRecord.addToPileupQueue(samPileupDeque_);
    samRecordDeque_.push_back(samRecord);

    while (samPileupDeque_.posMin() < samRecord.posMin) {
        int k = genotyper_.computeQuantizerIndex(samPileupDeque_.front().seq, samPileupDeque_.front().qual);
        mappedQuantizerIndices_.push_back(k);
//         std::cout << samPileupDeque_.front().pos << " ";
//         std::cout << samPileupDeque_.front().seq << " ";
//         std::cout << samPileupDeque_.front().qual << " > " << k << std::endl;
        samPileupDeque_.pop_front();
    }

    while (samRecordDeque_.front().posMax < samPileupDeque_.posMin()) {
        encodeMappedQual(samRecordDeque_.front());
        samRecordDeque_.pop_front();
    }

    nrMappedRecords_++;
}

void QualEncoder::finishBlock(void)
{
    // Compute all remaining quantizers
    while (samPileupDeque_.empty() == false) {
        int k = genotyper_.computeQuantizerIndex(samPileupDeque_.front().seq, samPileupDeque_.front().qual);
        mappedQuantizerIndices_.push_back(k);
        if (k < -1) throwErrorException("here i am");
//         std::cout << samPileupDeque_.front().pos << " ";
//         std::cout << samPileupDeque_.front().seq << " ";
//         std::cout << samPileupDeque_.front().qual << " > " << k << std::endl;
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
    CALQ_LOG("Writing block");

    size_t ret = 0;

    CALQ_LOG("Writing unmapped quality values");
    //std::cout << "uqv: " << unmappedQual_ << std::endl;
    unsigned char *uqv = (unsigned char *)unmappedQual_.c_str();
    size_t uqvSize = unmappedQual_.length();
    uint8_t uqvFlag = 0;
    if (uqvSize > 0) {
        uqvFlag = 1;
        cqFile->writeUint8(uqvFlag);
        ret += cqFile->writeQualBlock(uqv, uqvSize);
    } else {
        uqvFlag = 0;
        cqFile->writeUint8(uqvFlag);
        CALQ_LOG("No unmapped quality values in this block");
    }

    CALQ_LOG("Writing mapped quantizer indices");
    std::string tmp("");
    for (auto const &mappedQuantizerIndex : mappedQuantizerIndices_) {
        tmp += std::to_string(mappedQuantizerIndex+1);
    }
    //std::cout << "mqi: " << tmp << std::endl;
    unsigned char *mqi = (unsigned char *)tmp.c_str();
    size_t mqiSize = tmp.length();
    uint8_t mqiFlag = 0;
    if (mqiSize > 0) {
        mqiFlag = 1;
        cqFile->writeUint8(mqiFlag);
        ret += cqFile->writeQualBlock(mqi, mqiSize);
    } else {
        mqiFlag = 0;
        cqFile->writeUint8(mqiFlag);
        CALQ_LOG("No mapped quantizer indices in this block");
    }

    CALQ_LOG("Writing mapped quality value indices");
    tmp = "";
    for (auto const &mappedQualIndex : mappedQualIndices_) {
        tmp += std::to_string(mappedQualIndex);
    }
    //std::cout << "mqvi: " << tmp << std::endl;
    unsigned char *mqvi = (unsigned char *)tmp.c_str();
    size_t mqviSize = tmp.length();
    uint8_t mqviFlag = 0;
    if (mqviSize > 0) {
        mqviFlag = 1;
        cqFile->writeUint8(mqviFlag);
        ret += cqFile->writeQualBlock(mqvi, mqviSize);
    } else {
        mqviFlag = 0;
        cqFile->writeUint8(mqviFlag);
        CALQ_LOG("No mapped quality value in this block");
    }

    compressedMappedQualSize_ += ret;
    return ret;
}

size_t QualEncoder::compressedMappedQualSize(void) const { return compressedMappedQualSize_; }
size_t QualEncoder::compressedUnmappedQualSize(void) const { return compressedUnmappedQualSize_; }
size_t QualEncoder::compressedQualSize(void) const { return compressedMappedQualSize_ + compressedUnmappedQualSize_; }
size_t QualEncoder::nrMappedRecords(void) const { return nrMappedRecords_; }
size_t QualEncoder::nrUnmappedRecords(void) const { return nrUnmappedRecords_; }
size_t QualEncoder::nrRecords(void) const { return nrMappedRecords_ + nrUnmappedRecords_; }
size_t QualEncoder::uncompressedMappedQualSize(void) const { return uncompressedMappedQualSize_; }
size_t QualEncoder::uncompressedUnmappedQualSize(void) const { return uncompressedUnmappedQualSize_; }
size_t QualEncoder::uncompressedQualSize(void) const { return uncompressedMappedQualSize_ + uncompressedUnmappedQualSize_; }

void QualEncoder::encodeMappedQual(const SAMRecord &samRecord)
{
    size_t cigarIdx = 0;
    size_t cigarLen = samRecord.cigar.length();
    size_t opLen = 0; // length of current CIGAR operation
    size_t qualIdx = 0;
    size_t mqiIdx = samRecord.posMin - posOffset_;

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
               int mappedQuantizerIndex = mappedQuantizerIndices_[mqiIdx++];
//                std::cout << "mappedQuantizerIndices_[" <<mqiIdx-1<< "]="<< mappedQuantizerIndex << std::endl;
               int mappedQualIndex = quantizers_.at(mappedQuantizerIndex).valueToIndex(q);
               mappedQualIndices_.push_back(mappedQualIndex);
           }
           break;
       case 'I':
       case 'S':
           // Encode opLen quality values with max quantizer index
           for (size_t i = 0; i < opLen; i++) {
               int q = (int)samRecord.qual[qualIdx++] - qualityValueOffset_;
               int mappedQualIndex = quantizers_.at(quantizers_.size()-1).valueToIndex(q);
               mappedQualIndices_.push_back(mappedQualIndex);
           }
           break;
       case 'D':
       case 'N':
           mqiIdx += opLen;
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
    unmappedQual_ += qual;
}

} // namespace calq


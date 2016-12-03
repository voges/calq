/** @file QualDecoder.cc
 *  @brief This file contains the implementation of the QualDecoder class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */
#include "QualCodec/QualDecoder.h"

#include <map>
#include <string>

#include "Common/Exceptions.h"
#include "Common/log.h"

namespace calq {

QualDecoder::QualDecoder(const std::map<int, Quantizer> &quantizers)
    : posOffset_(0),
      qualityValueOffset_(0),
      unmappedQualityValues_(""),
      mappedQuantizerIndices_(""),
      mappedQualityValueIndices_(""),
      unmappedQualityValuesPosition_(0),
      mappedQualityValueIndicesPosition_(0),
      quantizers_(quantizers) {}

QualDecoder::~QualDecoder(void) {}

void QualDecoder::decodeMappedRecordFromBlock(const SAMRecord &samRecord, File *qualFile)
{
    size_t qualityValueIndicesLen = 0;

    size_t cigarIdx = 0;
    size_t cigarLen = samRecord.cigar.length();
    size_t opLen = 0; // length of current CIGAR operation

    for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++) {
       if (isdigit(samRecord.cigar[cigarIdx])) {
           opLen = opLen*10 + (size_t)samRecord.cigar[cigarIdx] - (size_t)'0';
           continue;
       }
       qualityValueIndicesLen += opLen;
    }

    std::string qualityValueIndices = mappedQualityValueIndices_.substr(mappedQualityValueIndicesPosition_, qualityValueIndicesLen);
    mappedQualityValueIndicesPosition_ += qualityValueIndicesLen;

    std::string qual("");

    cigarIdx = 0;
    cigarLen = samRecord.cigar.length();
    opLen = 0; // length of current CIGAR operation
//     size_t qualIdx = 0;
//     size_t mqiIdx = samRecord.posMin - posOffset_;

    for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++) {
       if (isdigit(samRecord.cigar[cigarIdx])) {
           opLen = opLen*10 + (size_t)samRecord.cigar[cigarIdx] - (size_t)'0';
           continue;
       }

       switch (samRecord.cigar[cigarIdx]) {
       case 'M':
       case '=':
       case 'X':
           // Decode opLen quality value indices with computed quantizer indices
           for (size_t i = 0; i < opLen; i++) {
//                int q = (int)samRecord.qual[qualIdx++] - qualityValueOffset_;
//                int mappedQuantizerIndex = mappedQuantizerIndices_[mqiIdx++];
//                int mappedQualityValueIndex = quantizers_.at(mappedQuantizerIndex).valueToIndex(q);
//                mappedQualityValueIndices_.push_back(mappedQualityValueIndex);
           }
           break;
       case 'I':
       case 'S':
           // Decode opLen quality values with max quantizer index
           for (size_t i = 0; i < opLen; i++) {
//                int q = (int)samRecord.qual[qualIdx++] - qualityValueOffset_;
//                int mappedQualityValueIndex = quantizers_.at(quantizers_.size()-1).valueToIndex(q);
//                mappedQualityValueIndices_.push_back(mappedQualityValueIndex);
           }
           break;
       case 'D':
       case 'N':
//            mqiIdx += opLen;
           break; // do nothing as these bases are not present
       case 'H':
       case 'P':
           break; // these have been clipped
       default:
           throwErrorException("Bad CIGAR string");
       }
       opLen = 0;
    }

    qualFile->write((unsigned char *)qual.c_str(), qual.length());
    qualFile->writeByte('\n');
}

void QualDecoder::decodeUnmappedRecordFromBlock(const SAMRecord &samRecord, File *qualFile)
{
    size_t qualLen = samRecord.seq.length();
    std::string qual = unmappedQualityValues_.substr(unmappedQualityValuesPosition_, qualLen);
    unmappedQualityValuesPosition_ += qualLen;
    qualFile->write((unsigned char *)qual.c_str(), qual.length());
    qualFile->writeByte('\n');
}

size_t QualDecoder::readBlock(CQFile *cqFile)
{
    CALQ_LOG("Reading block");

    size_t ret = 0;

//     CALQ_LOG("  Reading pos offset and quality value offset");
    ret += cqFile->readUint32(&posOffset_);
    ret += cqFile->readUint32((uint32_t *)&qualityValueOffset_);

//     CALQ_LOG("  Reading unmapped quality values");
    uint8_t uqvFlags = 0;
    ret += cqFile->readUint8(&uqvFlags);
    if (uqvFlags & 0x1) {
        ret += cqFile->readQualBlock(&unmappedQualityValues_);
    } else {
//         CALQ_LOG("  No unmapped quality values in this block");
    }
//     std::cout << "uqv: " << unmappedQualityValues_ << std::endl;

//     CALQ_LOG("  Reading mapped quantizer indices");
    uint8_t mqiFlags = 0;
    ret += cqFile->readUint8(&mqiFlags);
    if (mqiFlags & 0x1) {
        ret += cqFile->readQualBlock(&mappedQuantizerIndices_);
    } else {
        CALQ_LOG("  No mapped quantizer indices in this block");
    }
//     std::cout << "mqi: " << mappedQuantizerIndices_ << std::endl;

//     CALQ_LOG("  Reading mapped quality value indices");
    uint8_t mqviFlags = 0;
    ret += cqFile->readUint8(&mqviFlags);
    if (mqviFlags & 0x1) {
        ret += cqFile->readQualBlock(&mappedQualityValueIndices_);
    } else {
//         CALQ_LOG("  No mapped quality value indices in this block");
    }
//     std::cout << "mqvi: " << mappedQualityValueIndices_ << std::endl;

    return ret;
}

} // namespace calq


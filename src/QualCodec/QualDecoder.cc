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

static size_t readLength(const std::string &cigar)
{
    size_t readLen = 0;
    size_t cigarIdx = 0;
    size_t cigarLen = cigar.length();
    size_t opLen = 0; // length of current CIGAR operation

    for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++) {
       if (isdigit(cigar[cigarIdx])) {
           opLen = opLen*10 + (size_t)cigar[cigarIdx] - (size_t)'0';
           continue;
       }

       switch (cigar[cigarIdx]) {
       case 'M':
       case '=':
       case 'X':
           readLen += opLen;
           break;
       case 'I':
       case 'S':
           readLen += opLen;
           break;
       case 'D':
       case 'N':
           break; // do nothing as these bases are not present
       case 'H':
       case 'P':
           break; // these have been clipped
       default:
           throwErrorException("Bad CIGAR string");
       }

       opLen = 0;
    }

    return readLen;
}

QualDecoder::QualDecoder(const std::map<int, Quantizer> &quantizers)
    : posOffset_(0),
      qualityValueOffset_(0),
      unmappedQualityValues_(""),
      mappedQuantizerIndices_(""),
      mappedQualityValueIndices_(""),
      unmappedQualityValuesPosition_(0),
      mappedQualityValueIndicesPosition_(0),
      quantizers_(quantizers)
{
    CALQ_LOG("Initialized with %zu quantizers", quantizers_.size());
}

QualDecoder::~QualDecoder(void) {}

void QualDecoder::decodeMappedRecordFromBlock(const SAMRecord &samRecord, File *qualFile)
{
//     printf("Decoding SAM record\n");
//     samRecord.printShort();

    // Compute the read length from the CIGAR string
    size_t qualityValueIndicesLen = readLength(samRecord.cigar);

    // Get the quality value indices
    std::string qualityValueIndices = mappedQualityValueIndices_.substr(mappedQualityValueIndicesPosition_, qualityValueIndicesLen);
    if (qualityValueIndices.empty() == true) {
        throwErrorException("Decoding quality value indices failed");
    }
    mappedQualityValueIndicesPosition_ += qualityValueIndicesLen;

    // Reconstruct the quality values
    std::string qual("");

    size_t cigarIdx = 0;
    size_t cigarLen = samRecord.cigar.length();
    size_t opLen = 0;

    size_t quantizerIndicesIdx = samRecord.posMin - posOffset_;
    size_t qualityValueIndicesIdx = 0;

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
               int quantizerIndex = mappedQuantizerIndices_[quantizerIndicesIdx++] - 1 - '0';
               int qualityValueIndex = qualityValueIndices[qualityValueIndicesIdx++] - '0';
               int q = quantizers_.at(quantizerIndex).indexToReconstructionValue(qualityValueIndex);
//                printf("Decoded %c with k=%d -> %c\n", (char)(qualityValueIndex+'0'), quantizerIndex, (char)(q + qualityValueOffset_));
               qual += q + qualityValueOffset_;
           }
           break;
       case 'I':
       case 'S':
           // Decode opLen quality values with max quantizer index
           for (size_t i = 0; i < opLen; i++) {
               int quantizerIndex = quantizers_.size() - 1;
               int qualityValueIndex = qualityValueIndices[qualityValueIndicesIdx++] - '0';
               int q = quantizers_.at(quantizerIndex).indexToReconstructionValue(qualityValueIndex);
//                printf("Decoded %c with kmax=%d -> %c\n", (char)(qualityValueIndex+'0'), quantizerIndex, (char)(q + qualityValueOffset_));
               qual += q + qualityValueOffset_;
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

    qualFile->write((unsigned char *)qual.c_str(), qual.length());
    qualFile->writeByte('\n');
}

void QualDecoder::decodeUnmappedRecordFromBlock(const SAMRecord &samRecord, File *qualFile)
{
    // Compute the read length from the CIGAR string
    size_t qualLen = readLength(samRecord.cigar);

    // Get the quality values
    std::string qual = unmappedQualityValues_.substr(unmappedQualityValuesPosition_, qualLen);
    if (qual.empty() == true) {
        throwErrorException("Decoding quality values failed");
    }
    unmappedQualityValuesPosition_ += qualLen;

    // Write the quality values
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


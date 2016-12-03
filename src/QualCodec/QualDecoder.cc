/** @file QualDecoder.cc
 *  @brief This file contains the implementation of the QualDecoder class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */
#include "QualCodec/QualDecoder.h"

#include "Common/Exceptions.h"
#include "Common/log.h"
#include "Compressors/rle/rle.h"

namespace calq {

QualDecoder::QualDecoder(std::map<int,Quantizer> &quantizers)
    : unmappedQualityValues_(""),
      mappedQuantizerIndices_(""),
      mappedQualityValueIndices_(""),
      quantizers_(quantizers) {}

QualDecoder::~QualDecoder(void) {}

size_t QualDecoder::decodeMappedRecordFromBlock(const SAMRecord &samRecord, File *qualFile)
{
    size_t ret = 0;

    return ret;
}

size_t QualDecoder::decodeUnmappedRecordFromBlock(const SAMRecord &samRecord, File *qualFile)
{
    size_t ret = 0;

    return ret;
}

size_t QualDecoder::readBlock(CQFile *cqFile)
{
    CALQ_LOG("Reading block");

    size_t ret = 0;

    CALQ_LOG("Reading unmapped quality values");
    uint8_t uqvFlag = 0;
    cqFile->readUint8(&uqvFlag);
    if (uqvFlag & 0x1) {
        ret += cqFile->readQualBlock(&unmappedQualityValues_);
    } else {
        CALQ_LOG("No unmapped quality values in this block");
    }
    std::cout << "uqv: " << unmappedQualityValues_ << std::endl;

    CALQ_LOG("Reading mapped quantizer indices");
    uint8_t mqiFlag = 0;
    cqFile->readUint8(&mqiFlag);
    if (mqiFlag & 0x1) {
        ret += cqFile->readQualBlock(&mappedQuantizerIndices_);
    } else {
        CALQ_LOG("No mapped quantizer indices in this block");
    }
    std::cout << "mqi: " << mappedQuantizerIndices_ << std::endl;

    CALQ_LOG("Reading mapped quality value indices");
    uint8_t mqviFlag = 0;
    cqFile->readUint8(&mqviFlag);
    if (mqviFlag & 0x1) {
        ret += cqFile->readQualBlock(&mappedQualityValueIndices_);
    } else {
        CALQ_LOG("No mapped quality value indices in this block");
    }
    std::cout << "mqvi: " << mappedQualityValueIndices_ << std::endl;

    return ret;
}

} // namespace calq


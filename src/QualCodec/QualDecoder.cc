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
    : unmappedQualityValues_(""),
      mappedQuantizerIndices_(""),
      mappedQualityValueIndices_("") {}

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
    uint8_t flag = 0;
    cqFile->readUint8(&flag);
    if (flag)
        ret += cqFile->readQualBlock(&unmappedQualityValues_);
    std::cout << unmappedQualityValues_ << std::endl;

    CALQ_LOG("Reading mapped quantizer indices");
    ret += cqFile->readQualBlock(&mappedQuantizerIndices_);
    std::cout << mappedQuantizerIndices_ << std::endl;

    CALQ_LOG("Reading mapped quality value indices");
    ret += cqFile->readQualBlock(&mappedQualityValueIndices_);

    return ret;
}

} // namespace calq


/** @file QualDecoder.cc
 *  @brief This file contains the implementation of the QualDecoder class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */
#include "QualCodec/QualDecoder.h"

#include <map>

#include "Common/Exceptions.h"
#include "Common/log.h"

namespace calq {

QualDecoder::QualDecoder(const std::map<int, Quantizer> &quantizers)
    : posOffset_(0),
      qualityValueOffset_(0),
      unmappedQualityValues_(""),
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


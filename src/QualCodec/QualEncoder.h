/** @file QualEncoder.h
 *  @brief This file contains the definition of the QualEncoder class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#ifndef CALQ_QUALCODEC_QUALENCODER_H_
#define CALQ_QUALCODEC_QUALENCODER_H_

#include <chrono>
#include <deque>
#include <map>
#include <string>

#include "IO/CQ/CQFile.h"
#include "IO/SAM/SAMPileupDeque.h"
#include "IO/SAM/SAMRecord.h"
#include "QualCodec/Genotyper.h"
#include "QualCodec/Quantizer.h"

namespace calq {

class QualEncoder {
public:
    static const int QUANTIZER_STEPS_MIN = 2;
    static const int QUANTIZER_STEPS_MAX = 5;
    static const int NR_QUANTIZERS = QUANTIZER_STEPS_MAX-QUANTIZER_STEPS_MIN+1;
    static const int QUANTIZER_IDX_MIN = 0;
    static const int QUANTIZER_IDX_MAX = NR_QUANTIZERS-1;

    explicit QualEncoder(const int &polyploidy,
                         const int &qualityValueMax,
                         const int &qualityValueMin,
                         const int &qualityValueOffset,
                         const std::map<int, Quantizer> &quantizers);
    ~QualEncoder(void);

    void addUnmappedRecordToBlock(const SAMRecord &samRecord);
    void addMappedRecordToBlock(const SAMRecord &samRecord);
    void finishBlock(void);
    size_t writeBlock(CQFile *cqFile);

    size_t compressedMappedQualSize(void) const;
    size_t compressedUnmappedQualSize(void) const;
    size_t compressedQualSize(void) const;
    size_t nrMappedRecords(void) const;
    size_t nrUnmappedRecords(void) const;
    size_t nrRecords(void) const;
    size_t uncompressedMappedQualSize(void) const;
    size_t uncompressedUnmappedQualSize(void) const;
    size_t uncompressedQualSize(void) const;

private:
    void encodeMappedQual(const SAMRecord &samRecord);
    void encodeUnmappedQual(const std::string &qual);

    // Timekeeping
    std::chrono::time_point<std::chrono::steady_clock> startTime_;
    std::chrono::time_point<std::chrono::steady_clock> stopTime_;

    // Sizes & counters
    size_t compressedMappedQualSize_;
    size_t compressedUnmappedQualSize_;
    size_t nrMappedRecords_;
    size_t nrUnmappedRecords_;
    size_t uncompressedMappedQualSize_;
    size_t uncompressedUnmappedQualSize_;

    int qualityValueOffset_;

    // 0-based position offset of this block
    uint32_t posOffset_;

    // Buffers
    std::string unmappedQualityValues_;
    std::deque<int> mappedQuantizerIndices_;
    std::deque<int> mappedQualityValueIndices_;

    // Pileup
    SAMPileupDeque samPileupDeque_;

    // Genotyper
    Genotyper genotyper_;

    // Quantizers
    std::map<int, Quantizer> quantizers_;

    // Double-ended queue holding the SAM records; records get popped when they
    // are finally encoded
    std::deque<SAMRecord> samRecordDeque_;
};

} // namespace calq

#endif // CALQ_QUALCODEC_QUALENCODER_H_


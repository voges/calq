/** @file QualEncoder.h
 *  @brief This file contains the definition of the QualEncoder class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef CQ_QUALENCODER_H
#define CQ_QUALENCODER_H

#include "IO/CQ/CQFile.h"
#include "IO/SAM/SAMPileupDeque.h"
#include "IO/SAM/SAMRecord.h"
#include "QualCodec/Genotyper.h"
#include <chrono>
#include <deque>

namespace cq {

class QualEncoder {
public:
    const unsigned int QUANTIZER_STEPS_MIN = 2;
    const unsigned int QUANTIZER_STEPS_MAX = 8;
    const unsigned int NUM_QUANTIZERS = QUANTIZER_STEPS_MAX-QUANTIZER_STEPS_MIN+1;
    const unsigned int QUANTIZER_IDX_MIN = 0;
    const unsigned int QUANTIZER_IDX_MAX = NUM_QUANTIZERS-1;

public:
    explicit QualEncoder(const unsigned int &polyploidy,
                         const int &qMin,
                         const int &qMax);
    ~QualEncoder(void);

    void startBlock(void);
    void addUnmappedRecordToBlock(const SAMRecord &samRecord);
    void addMappedRecordToBlock(const SAMRecord &samRecord);
    size_t finishAndWriteBlock(CQFile &cqFile);
    void printBlockStatistics(void) const;

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

private:
    // Timekeeping
    std::chrono::time_point<std::chrono::steady_clock> blockStartTime_;
    std::chrono::time_point<std::chrono::steady_clock> blockStopTime_;

    // Sizes & counters
    size_t compressedMappedQualSize_;
    size_t compressedUnmappedQualSize_;
    size_t nrMappedRecords_;
    size_t nrUnmappedRecords_;
    size_t uncompressedMappedQualSize_;
    size_t uncompressedUnmappedQualSize_;

    // 0-based position offset of this block
    uint32_t posOff_;

    // Buffers
    std::string unmappedQual_;
    std::string mappedQuantizerIndices_;
    std::string mappedQualIndices_;

    // Pileup
    SAMPileupDeque samPileupDeque_;

    // Genotyper
    Genotyper genotyper_;

    // Quantizers
    //std::map<int,Quantizer> quantizers_;

    // Double-ended queue holding the SAM records; records get popped when they
    // are finally encoded
    std::deque<SAMRecord> samRecordDeque_;
};

}

#endif // CQ_QUALENCODER_H


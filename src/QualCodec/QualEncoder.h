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
#include "IO/SAM/SAMRecord.h"
#include "QualCodec/Genotyper.h"
#include "QualCodec/Pileup.h"
//#include "QualCodec/UniformQuantizer.h"
#include <chrono>
#include <deque>

namespace cq {

class QualEncoder {
public:
    static const unsigned int QUANTIZER_STEPS_MIN = 2;
    static const unsigned int QUANTIZER_STEPS_MAX = 8;
    static const unsigned int NUM_QUANTIZERS = QUANTIZER_STEPS_MAX-QUANTIZER_STEPS_MIN+1;
    static const unsigned int QUANTIZER_IDX_MIN = 0;
    static const unsigned int QUANTIZER_IDX_MAX = NUM_QUANTIZERS-1;

public:
    explicit QualEncoder(const unsigned int &polyploidy,
                         const int &qMin,
                         const int &qMax
    );
    ~QualEncoder(void);

    void startBlock(void);
    void addUnmappedRecordToBlock(const SAMRecord &samRecord);
    void addMappedRecordToBlock(const SAMRecord &samRecord);
    size_t finishAndWriteBlock(CQFile &cqFile);
    void printBlockStatistics(void) const;

    size_t compressedMappedQualSize(void) const;
    size_t compressedUnmappedQualSize(void) const;
    size_t compressedQualSize(void) const;
    size_t numMappedRecords(void) const;
    size_t numUnmappedRecords(void) const;
    size_t numRecords(void) const;
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
    size_t numMappedRecords_;
    size_t numUnmappedRecords_;
    size_t uncompressedMappedQualSize_;
    size_t uncompressedUnmappedQualSize_;

    // Buffers
    std::string m_unmappedQual_
    std::deque<int> mappedQuantizerIndices_;
    uint32_t mappedQuantizerIndicesPosMin_;
    uint32_t mappedQuantizerIndicesPosMax_;
    std::deque<int> mappedQualIndices_;

    // Pileup
    PileupQueue pileupQueue_;

    // Genotyper
    Genotyper genotyper_;

    // Quantizers
    //std::map<int,Quantizer> quantizers_;

    // Double-ended queue holding the SAM records; records get popped when they
    // are finally encoded
    std::deque<SAMRecord> samRecordQueue_;
};

}

#endif // CQ_QUALENCODER_H


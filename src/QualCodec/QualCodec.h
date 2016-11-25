/** @file QualCodec.h
 *  @brief This file contains the definitions of the QualEncoder and
 *         QualDecoder classes.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef CQ_QUALCODEC_H
#define CQ_QUALCODEC_H

#include "IO/CQ/CQFile.h"
#include "IO/SAM/SAMRecord.h"
#include "QualCodec/Genotyper.h"
//#include "QualCodec/UniformQuantizer.h"
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
    std::chrono::time_point<std::chrono::steady_clock> m_blockStartTime;
    std::chrono::time_point<std::chrono::steady_clock> m_blockStopTime;

    // Sizes & counters
    size_t m_compressedMappedQualSize;
    size_t m_compressedUnmappedQualSize;
    size_t m_numMappedRecords;
    size_t m_numUnmappedRecords;
    size_t m_uncompressedMappedQualSize;
    size_t m_uncompressedUnmappedQualSize;

    // Buffers
    std::string m_unmappedQual;
    std::deque<int> m_mappedQuantizerIndices;
    uint32_t m_mappedQuantizerIndicesPosMin;
    uint32_t m_mappedQuantizerIndicesPosMax;
    std::deque<int> m_mappedQualIndices;

    // Pileup
    std::deque<std::string> m_seqPileup;
    std::deque<std::string> m_qualPileup;
    uint32_t m_pileupPosMin;
    uint32_t m_pileupPosMax;

    // Genotyper
    Genotyper m_genotyper;

    // Quantizers
    //std::map<int,Quantizer> m_quantizers;

    // Double-ended queue holding the SAM records; records get popped when they
    // are finally encoded
    std::deque<SAMRecord> m_samRecordQueue;
};

class QualDecoder {
public:
   explicit QualDecoder(void);
   ~QualDecoder(void);
};

}

#endif // CQ_QUALCODEC_H


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

#include "IO/File.h"
//#include "Genotyper/Genotyper.h"
//#include "Quantizers/UniformQuantizer.h"
#include "IO/SAM/SAMRecord.h"
#include <chrono>
//#include <fstream>
//#include <map>
//#include <math.h>
//#include <queue>
//#include <vector>
//

namespace cq {

class QualEncoder {
public:
    QualEncoder(File &cqFile, const unsigned int &polyploidy);
    ~QualEncoder(void);

    void startBlock(void);
    void addUnmappedRecordToBlock(const SAMRecord &samRecord);
    //void addMappedRecordToBlock(const SAMRecord &samRecord);
    size_t finishBlock(void);

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
    void encodeUnmappedQual(const std::string &qual);

private:
    std::chrono::time_point<std::chrono::steady_clock> m_blockStartTime;
    std::chrono::time_point<std::chrono::steady_clock> m_blockStopTime;

    size_t m_compressedMappedQualSize;
    size_t m_compressedUnmappedQualSize;
    size_t m_numMappedRecords;
    size_t m_numUnmappedRecords;
    size_t m_uncompressedMappedQualSize;
    size_t m_uncompressedUnmappedQualSize;
    
    std::string m_unmappedQual;
    std::string m_mappedQuantizerIndices;
    std::string m_mappedQualIndices;

//    void printStats(void) const;
//
    
    
//
//
//private:
//    void encodeMappedQualityValues(const MappedRecord &mappedRecord);
//    

private:
//    ////////////////////////////////////////////////////////////////////////////
//    // Class scope:
//    // These member variables are used throughout coding multiple blocks
//    ////////////////////////////////////////////////////////////////////////////
//
//    // Generic
//    bool quantizedPrintout;
//    bool stats;
//    bool verbose;
//    int qvMin;
//    int qvMax;
//
//    std::ofstream qualMappedFile;
//    std::ofstream qualUnmappedFile;
//    std::ofstream statsFile;
//
//    File &cqFile;
//    Genotyper genotyper;
//    std::map<int,UniformQuantizer> uniformQuantizers;
//
//    // Sizes & counters
//    size_t uncompressedSize;
//    size_t uncompressedMappedSize;
//    size_t uncompressedUnmappedSize;
//    size_t compressedSize;
//    size_t compressedMappedSize;
//    size_t compressedUnmappedSize;
//    size_t numRecords;
//    size_t numMappedRecords;
//    size_t numUnmappedRecords;
//
//    ////////////////////////////////////////////////////////////////////////////
//    // Block scope:
//    // The following member variables are used per block; all of them are reset
//    // by startBlock.
//    ////////////////////////////////////////////////////////////////////////////
//
//    // Start and stop times for a block
//    std::chrono::time_point<std::chrono::steady_clock> blockStartTime;
//    std::chrono::time_point<std::chrono::steady_clock> blockStopTime;
//
//    // Sizes & counters
//    size_t uncompressedSizeOfBlock;
//    size_t uncompressedMappedSizeOfBlock;
//    size_t uncompressedUnmappedSizeOfBlock;
//    size_t compressedSizeOfBlock;
//    size_t compressedMappedSizeOfBlock;
//    size_t compressedUnmappedSizeOfBlock;
//    size_t numRecordsInBlock;
//    size_t numMappedRecordsInBlock;
//    size_t numUnmappedRecordsInBlock;
//
//    // Strings to buffer the data to be transmitted
//    std::string qi; // quantizer indices
//    std::string qvi; // quality value indices
//    std::string uqv; // unmapped quality values
//
//    // Reference loaded for the current block
//    std::string referenceName;
//    std::string reference;
//    uint32_t referencePosMin;
//    uint32_t referencePosMax;
//
//    // Queue holding the mapped records; records get popped when they are
//    // finally encoded
//    std::queue<MappedRecord> mappedRecordQueue;
//
//    // Current observed nucleotides and quality values; when a quantizer index
//    // is computed for a certain position, the vectors are shrinked
//    std::deque<std::string> observedNucleotides;
//    std::deque<std::string> observedQualityValues;
//    uint32_t observedPosMin;
//    uint32_t observedPosMax;
//
//    // Computed quantizer indices; when the indices are not needed anymore; the
//    // vector is shrinked
//    std::vector<int> quantizerIndices;
//    uint32_t quantizerIndicesPosMin;
//    uint32_t quantizerIndicesPosMax;
//
};
//
//class QualDecoder {
//public:
//    QualDecoder(File &cqFile, File &qualFile, const CLIOptions &cliOptions);
//    ~QualDecoder(void);
//
//    void printStats(void) const;
//
//    size_t decodeBlock(void);
//
//private:
//    ////////////////////////////////////////////////////////////////////////////
//    // Class scope:
//    // These member variables are used throughout decoding multiple blocks
//    ////////////////////////////////////////////////////////////////////////////
//
//    // Generic
//    bool verbose;
//    File &cqFile;
//    File &qualFile;
//
//    // Sizes & counters
//    size_t uncompressedSize;
//    size_t compressedSize;
//    size_t numBlocks;
//    size_t numRecords;
//    size_t numMappedRecords;
//    size_t numUnmappedRecords;
//
//    ////////////////////////////////////////////////////////////////////////////
//    // Block scope:
//    // The following member variables are used per block; all of them are reset
//    // by decodeBlock.
//    ////////////////////////////////////////////////////////////////////////////
//
//    // Sizes & counters
//    size_t uncompressedSizeOfBlock;
//    size_t compressedSizeOfBlock;
//    size_t numRecordsInBlock;
//    size_t numMappedRecordsInBlock;
//    size_t numUnmappedRecordsInBlock;
//};

}

#endif // CQ_QUALCODEC_H


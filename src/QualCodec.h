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

#ifndef QUALCODEC_H
#define QUALCODEC_H

#include "Common/File.h"
#include "Genotyper/Genotyper.h"
#include "Parsers/FASTAReference.h"
#include "Quantizers/UniformQuantizer.h"
#include "Records/SAMRecord.h"
#include "Records/MappedRecord.h"
#include <math.h>
#include <fstream>
#include <map>
#include <queue>
#include <vector>

class QualEncoder {
public:
    QualEncoder(File &cqFile, 
                const unsigned int &polyploidy,
                const int &qvMin,
                const int &qvMax);
    ~QualEncoder(void);

    void startBlock(void);
    void addUnmappedRecordToBlock(const SAMRecord &samRecord);
    void addMappedRecordToBlock(const SAMRecord &samRecord);
    size_t finishBlock(void);

public:
    std::vector<FASTAReference> fastaReferences;

private:
    ////////////////////////////////////////////////////////////////////////////
    // Class scope:
    // These member variables are used throughout coding multiple blocks
    ////////////////////////////////////////////////////////////////////////////

    // Generic
    bool encoderStats;
    bool quantizedPrintout;
    bool verbose;
    int qvMin;
    int qvMax;

    File &cqFile;
    Genotyper genotyper;
    std::map<int,UniformQuantizer> uniformQuantizers;

    // Sizes & counters
    size_t uncompressedSize;
    size_t compressedSize;
    size_t numBlocks;
    size_t numRecords;
    size_t numMappedRecords;
    size_t numUnmappedRecords;

    ////////////////////////////////////////////////////////////////////////////
    // Block scope:
    // The following member variables are used per block; all of them are reset
    // by startBlock.
    ////////////////////////////////////////////////////////////////////////////

    // Sizes & counters
    size_t uncompressedSizeOfBlock;
    size_t compressedSizeOfBlock;
    size_t numRecordsInBlock;
    size_t numMappedRecordsInBlock;
    size_t numUnmappedRecordsInBlock;

    // Strings to buffer the data to be transmitted
    std::string qi; // quantizer indices
    std::string qvi; // quality value indices
    std::string uqv; // unmapped quality values

    // Reference loaded for the current block
    std::string referenceName;
    std::string reference;
    uint32_t referencePosMin;
    uint32_t referencePosMax;

    // Queue holding the mapped records; records get popped when they are
    // finally encoded
    std::queue<MappedRecord> mappedRecordQueue;

    // Current observed nucleotides and quality values; when a quantizer index
    // is computed for a certain position, the vectors are shrinked
    std::vector<std::string> observedNucleotides;
    std::vector<std::string> observedQualityValues;
    uint32_t observedPosMin;
    uint32_t observedPosMax;

    // Computed quantizer indices; when the indices are not needed anymore; the
    // vector is shrinked
    std::vector<int> quantizerIndices;
    uint32_t quantizerIndicesPosMin;
    uint32_t quantizerIndicesPosMax;

private:
    void loadFastaReference(const std::string &rname);
    void encodeMappedQualityValues(const MappedRecord &mappedRecord);
    void encodeUnmappedQualityValues(const std::string &qualityValues);
};

class QualDecoder {
public:
    QualDecoder(File &cqFile, File &qualFile);
    ~QualDecoder(void);

    size_t decodeBlock(void);

private:
    ////////////////////////////////////////////////////////////////////////////
    // Class scope:
    // These member variables are used throughout decoding multiple blocks
    ////////////////////////////////////////////////////////////////////////////

    // Generic
    bool verbose;
    File &cqFile;
    File &qualFile;

    // Sizes & counters
    size_t uncompressedSize;
    size_t compressedSize;
    size_t numBlocks;
    size_t numRecords;
    size_t numMappedRecords;
    size_t numUnmappedRecords;

    ////////////////////////////////////////////////////////////////////////////
    // Block scope:
    // The following member variables are used per block; all of them are reset
    // by decodeBlock.
    ////////////////////////////////////////////////////////////////////////////

    // Sizes & counters
    size_t uncompressedSizeOfBlock;
    size_t compressedSizeOfBlock;
    size_t numRecordsInBlock;
    size_t numMappedRecordsInBlock;
    size_t numUnmappedRecordsInBlock;
};

#endif // QUALCODEC_H


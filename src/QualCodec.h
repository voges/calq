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

/** @brief Class: QualEncoder
 *
 *  The QualEncoder class encodes the quality scores.
 */
class QualEncoder {
public:
    QualEncoder(File &cqFile, 
                const std::vector<FASTAReference> &fastaReferences,
                const unsigned int &polyploidy);
    ~QualEncoder(void);

    void startBlock(void);
    void addUnmappedRecordToBlock(const SAMRecord &samRecord);
    void addMappedRecordToBlock(const SAMRecord &samRecord);
    size_t finishBlock(void);

private:
    ////////////////////////////////////////////////////////////////////////////
    // These member variables are used throughout coding multiple blocks
    ////////////////////////////////////////////////////////////////////////////

    std::vector<FASTAReference> fastaReferences;
    size_t numBlocks;
    size_t numMappedRecords;
    size_t numUnmappedRecords;
    File &cqFile;

    ////////////////////////////////////////////////////////////////////////////
    // The following member variables are used per block; all of them are reset
    // by startBlock.
    ////////////////////////////////////////////////////////////////////////////

    // Counters
    size_t numMappedRecordsInBlock;
    size_t numUnmappedRecordsInBlock;

    // Strings (i.e. byte buffers) holding the data to be transmitted
    std::string qi; // quantizer indices
    std::string qvi; // quality value indices
    std::string uqv; // unmapped quality values

    // Reference loaded for the current block
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

    // Class to perform the magic ;)
    Genotyper genotyper;

    // Quantizers
    std::map<int,UniformQuantizer> uniformQuantizers;

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

/** @brief Class: QualDecoder
 *
 *  The QualDecoder class decodes the encoded bitstream.
 */
class QualDecoder {
public:
    QualDecoder(File &cqFile,
                File &qualFile,
                const std::vector<FASTAReference> &fastaReferences);
    ~QualDecoder(void);

    void decodeBlock(void);

private:
    File &cqFile;
    std::vector<FASTAReference> fastaReferences;
    File &qualFile;
};

#endif // QUALCODEC_H


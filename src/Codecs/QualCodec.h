/** @file QualCodec.h
 *  @brief This file contains the definitions of the QualEncoder and
 *         QualDecoder classes.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: what (who)
 */

#ifndef QUALCODEC_H
#define QUALCODEC_H

#include "bitstream.h"
#include "Codecs/MappedRecord.h"
#include "Parsers/FASTAReference.h"
#include "Parsers/SAMRecord.h"
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
    QualEncoder(ofbitstream &ofbs, const std::vector<FASTAReference> &fastaReferences);
    ~QualEncoder(void);

    void startBlock(void);
    void addUnmappedRecordToBlock(const SAMRecord &samRecord);
    void addMappedRecordToBlock(const SAMRecord &samRecord);
    size_t finishBlock(void);

private:
    // these member variables are used throughout coding multiple blocks
    std::vector<FASTAReference> fastaReferences;
    size_t numBlocks;
    size_t numMappedRecords;
    size_t numUnmappedRecords;
    ofbitstream &ofbs;

    // these member variables are used per block
    std::string reference;
    std::queue<MappedRecord> mappedRecordQueue;
    uint32_t maxPosition;
    uint32_t minPosition;
    std::vector<std::string> observedNucleotides;
    std::vector<std::string> observedQualityValues;
    size_t observedSize;
    std::vector<int> quantizerIndices;

private:
    void loadFastaReference(const std::string &rname);
    void encodeMappedRecord(const MappedRecord &mappedRecord);
    void encodeUnmappedRecord(const std::string &qual);
};

/** @brief Class: QualDecoder
 *
 *  The QualDecoder class decodes the encoded bitstream.
 */
class QualDecoder {
public:
    QualDecoder(ifbitstream &ifbs, std::ofstream &ofs, const std::vector<FASTAReference> &fastaReferences);
    ~QualDecoder(void);

private:
    std::vector<FASTAReference> fastaReferences;
    ifbitstream &ifbs;
    std::ofstream &ofs;
};

#endif // QUALCODEC_H


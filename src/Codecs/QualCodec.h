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

#include "Compressors/bitstream.h"
#include "Compressors/caac/modelA.h"
#include "Compressors/caac/compressor.h"
#include "Compressors/caac/decompressor.h"
#include "Genotyper/Genotyper.h"
#include "Genotyper/Genotyper2.h"
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
    QualEncoder(ofbitstream &ofbs, 
                const std::vector<FASTAReference> &fastaReferences,
                const int &polyploidy);
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
    
    //ari member variables
    modelA<int, 16, 14> cmodel;
    Compressor<modelA<int, 16, 14>, ofbitstream> caac;


    // these member variables are used per block

    // reference loaded for the current block
    std::string reference;
    uint32_t referencePosMin;
    uint32_t referencePosMax;

    // queue holding the mapped records; records get popped when they are
    // finally encoded
    std::queue<MappedRecord> mappedRecordQueue;

    // current observed nucleotides and quality values; when a quantizer index
    // is computed for a certain position, the vectors are shrinked
    std::vector<std::string> observedNucleotides;
    std::vector<std::string> observedQualityValues;
    uint32_t observedPosMin;
    uint32_t observedPosMax;

    // class to perform the magic
    Genotyper genotyper;
    Genotyper2 genotyper2;

    // quantizers
    std::map<int,UniformQuantizer> uniformQuantizers;

    // computed quantizer indices; when the indices are not needed anymore; the
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
    QualDecoder(ifbitstream &ifbs, std::ofstream &ofs, const std::vector<FASTAReference> &fastaReferences);
    ~QualDecoder(void);

    void decodeBlock(void);

private:
    std::vector<FASTAReference> fastaReferences;
    ifbitstream &ifbs;
    std::ofstream &ofs;
    
    //ari member variables
    modelA<int, 16, 14> cmodel;
    Decompressor<modelA<int, 16, 14>, ifbitstream> caad;

};

#endif // QUALCODEC_H


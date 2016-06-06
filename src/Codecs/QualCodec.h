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
#include "Codecs/Predictor.h"
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
    QualEncoder::QualEncoder(ofbitstream &ofbs, const std::vector<FASTAReference> &fastaReferences);
    ~QualEncoder(void);

    void startBlock(const std::string &rname);
    void addRecordToBlock(const SAMRecord &samRecord);
    size_t finishBlock(void);

private:
    // these member variables are used throughout coding multiple blocks
    std::vector<FASTAReference> fastaReferences;
    size_t numEncodedRecords;
    ofbitstream &ofbs;

    // these member variables are used per block
    Predictor predictor;
    std::string rnamePrev;
    FASTAReference currFastaReference;
    std::map<size_t, uint8_t> quantizers;
    std::queue<SAMRecord> recordBuffer;

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


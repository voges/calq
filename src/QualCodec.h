/** @file QualCodec.h
 *  @brief This file contains the definitions of the QualEncoder and
 *         QualDecoder classes, respectively.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#ifndef QUALCODEC_H
#define QUALCODEC_H

#include "bitstream.h"
#include "Predictor.h"
#include "SAMRecord.h"
#include <vector>

/** @brief Class: QualEncoder
 *
 *  The QualEncoder class provides two methods as interface:
 *  - encodeRecord: This function is used to encode the quality scores in a
 *                  given SAM record.
 *  - finishBlock:  Upon having processed some number of SAM records, this
 *                  function can be triggered to finish a block. Some metadata
 *                  is written to the given stream and the QualEncoder instance
 *                  is reset. A new block can be triggered with another call
 *                  to encodeRecord.
 */
class QualEncoder {
public:
    QualEncoder(ofbitstream &ofbs);
    ~QualEncoder(void);
    void encodeRecord(const SAMRecord &samRecord);
    size_t finishBlock(void);
    void createCSV(void);

private:
    ofbitstream &ofbs;
    size_t recordCnt; ///< number of records processed in the current block
    Predictor predictor;
};

/** @brief Class: QualDecoder
 *
 *  The QualDecoder class provides one method:
 *  - decodeBlock: Read and decode a block of quality score from ifbs; the
 *                 decoded quality score are stored in qual.
 */
class QualDecoder {
public:
    QualDecoder(ifbitstream &ifbs);
    ~QualDecoder(void);
    void decodeBlock(std::vector<std::string> &qual);

private:
    ifbitstream &ifbs;
};

#endif // QUALCODEC_H


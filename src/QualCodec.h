/** @file QualCodec.h
 *  @brief
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#ifndef QUALCODEC_H
#define QUALCODEC_H

#include "bitstream.h"
#include "SAMRecord.h"
#include <vector>

/** @brief Class: QualEncoder
 *
 *  The QualEncoder class provides two methods:
 *  - encodeRecord: This function is used to encode the quality scores in a
 *                  given SAM record.
 *  - finishBlock:  Upon having processed some number of SAM records, this
 *                  function can be triggered to finish a block. Some metadata
 *                  is written to the given file and the QualEncoder instance
 *                  is reset. A new block can be triggered with another call
 *                  to encodeRecord.
 */
class QualEncoder {
private:
    size_t recordCnt; ///< number of records processed in the current block
    size_t posMin;    ///< smallest mapping position encountered in block
    size_t posMax;    ///< largest mapping position encountered in block
    std::vector<size_t> depths;
    ofbitstream &ofbs;

public:
    QualEncoder(ofbitstream &ofbs);
    ~QualEncoder(void);

    void encodeRecord(const SAMRecord &samRecord);
    size_t finishBlock(void);
};

/** @brief Class: QualDecoder
 *
 *  The QualDecoder class provides one methods:
 *  - decodeBlock: ...
 */
class QualDecoder {
private:

public:
    QualDecoder(ifbitstream &ifbs);
    ~QualDecoder(void);

    void decodeBlock(std::vector<std::string> &qual);
};

#endif // QUALCODEC_H


/** @file CalqCodec.h
 *  @brief This file contains the definitions of the CalqEncoder,
 *         CalqDecoder, and CalqInfoTool classes.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#ifndef CALQCODEC_H
#define CALQCODEC_H

#include "bitstream.h"
//#include "FASTAParser.h"
#include "QualCodec.h"
#include "SAMParser.h"
#include <fstream>
#include <string>

/** @brief Class: CalqEncoder
 *
 *  The CalqEncoder provides only one interface to the outside which is the
 *  member function encode; it encodes the SAM file with the name infileName
 *  with the specified blockSize and writes the encoded bitstream to the CQ
 *  file with the name outfileName.
 */
class CalqEncoder {
public:
    CalqEncoder(const std::string &inFileName,
                const std::string &outFileName,
                const std::string &referenceFileName);
    ~CalqEncoder(void);

    void encode(void);

private:
    std::vector<std::streampos> blockPositionList;
    const size_t blockSize;
    //FASTAParser fastaParser;
    const std::string inFileName;
    ofbitstream ofbs;
    const std::string outFileName;
    QualEncoder qualEncoder;
    const std::string referenceFileName;
    SAMParser samParser;
};

/** @brief Class: CalqDecoder
 *
 *  The CalqDecoder provides only one interface to the outside which is the
 *  member function decode; it decodes the CQ file with the name infileName
 *  and writes the decoded quality scores to the file with the name outfileName.
 */
class CalqDecoder {
public:
    CalqDecoder(const std::string &inFileName,
                const std::string &outFileName,
                const std::string &referenceFileName);
    ~CalqDecoder(void);

    void decode(void);

private:
    //FASTAParser fastaParser;
    ifbitstream ifbs;
    const std::string inFileName;
    std::ofstream ofs;
    const std::string outFileName;
    QualDecoder qualDecoder;
    const std::string referenceFileName;
};

#endif // CALQCODEC_H


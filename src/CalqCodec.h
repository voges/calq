/** @file CalqCodec.h
 *  @brief This file contains the definitions of the CalqEncoder,
 *         CalqDecoder, and CalqInfoTool classes.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#ifndef CALQCODEC_H
#define CALQCODEC_H

#include "bitstream.h"
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
    CalqEncoder(const std::string &infileName, const std::string &outfileName, const size_t &blockSize);
    ~CalqEncoder(void);
    void encode(void);

private:
    std::string infileName;
    std::string outfileName;
    SAMParser samParser;
    ofbitstream ofbs;
    size_t blockSize;
    QualEncoder qualEncoder;
    std::vector<std::streampos> blockPositionList;
};

/** @brief Class: CalqDecoder
 *
 *  The CalqDecoder provides only one interface to the outside which is the
 *  member function decode; it decodes the CQ file with the name infileName
 *  and writes the decoded quality scores to the file with the name outfileName.
 */
class CalqDecoder {
public:
    CalqDecoder(const std::string &infileName, const std::string &outfileName);
    ~CalqDecoder(void);
    void decode(void);

private:
    std::string infileName;
    std::string outfileName;
    ifbitstream ifbs;
    std::ofstream ofs;
    QualDecoder qualDecoder;
};

/** @brief Class: CalqInfoTool
 *
 *  The CalqInfoTool provides only one interface to the outside which is the
 *  member function extractInfo; it reads meta information from a CQ file and
 *  writes it to stdout.
 */
class CalqInfoTool {
public:
    CalqInfoTool(const std::string &infileName);
    ~CalqInfoTool(void);
    void extractInfo(void);

private:
    std::string infileName;
    ifbitstream ifbs;
};

#endif // CALQCODEC_H


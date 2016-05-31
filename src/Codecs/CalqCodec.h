/** @file CalqCodec.h
 *  @brief This file contains the definitions of the CalqEncoder and
 *         CalqDecoder classes.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#ifndef CALQCODEC_H
#define CALQCODEC_H

#include "bitstream.h"
#include "Parsers/FASTAParser.h"
#include "Parsers/SAMParser.h"
#include "QualCodec.h"
#include <fstream>
#include <string>
#include <vector>

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
                const std::vector<std::string> &referenceFileNames);
    ~CalqEncoder(void);

    void encode(void);

private:
    FASTAParser fastaParser;
    const std::string samFileName;
    ofbitstream ofbs;
    const std::string cqFileName;
    QualEncoder qualEncoder;
    const std::vector<std::string> fastaFileNames;
    std::vector<FASTAReference> fastaReferences;
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
                const std::vector<std::string> &referenceFileNames);
    ~CalqDecoder(void);

    void decode(void);

private:
    //FASTAParser fastaParser;
    ifbitstream ifbs;
    const std::string inFileName;
    std::ofstream ofs;
    const std::string outFileName;
    QualDecoder qualDecoder;
    const std::vector<std::string> referenceFileNames;
};

#endif // CALQCODEC_H


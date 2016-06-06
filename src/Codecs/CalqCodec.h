/** @file CalqCodec.h
 *  @brief This file contains the definitions of the CalqEncoder and
 *         CalqDecoder classes.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: what (who)
 */

#ifndef CALQCODEC_H
#define CALQCODEC_H

#include "bitstream.h"
#include "Codecs/QualCodec.h"
#include "Parsers/FASTAParser.h"
#include "Parsers/SAMParser.h"
#include <fstream>
#include <string>
#include <vector>

/** @brief Class: CalqCodec
 *
 *  The CalqEncoder and the CalqDecoder classes inherit common methods from
 *  this class.
 */
class CalqCodec {
public:
    /** @brief Constructor: CalqCodec
     *
     *  Initializes a new CalqCodec instance and reads the reference
     *  sequence(s) from all provided FASTA files. The reference
     *  sequence(s) are stored in 'fastaReferences'.
     *
     *  @param inFileName File name of the input file to be en- or decoded
     *  @param outfileName File name of the output file
     *  @param fastaFileNames Vector containing the file names of the FASTA
     *         which contain the reference sequence(s) used during en- or
     *         decoding
     */
    CalqCodec(const std::string &inFileName,
              const std::string &outFileName,
              const std::vector<std::string> &fastaFileNames);
              
    /** @brief Destructor: CalqCodec
     *
     *  Destructs a CalqCodec instance.
     */
    virtual ~CalqCodec(void);

protected:
    const std::vector<std::string> fastaFileNames;
    std::vector<FASTAReference> fastaReferences;
    const std::string inFileName;
    const std::string outFileName;

private:
    FASTAParser fastaParser;
};

/** @brief Class: CalqEncoder
*
*  The CalqEncoder provides only one interface to the outside which is the
*  member function encode; it encodes the SAM file with the name infileName
*  with the specified blockSize and writes the encoded bitstream to the CQ
*  file with the name outfileName.
*/
class CalqEncoder: public CalqCodec {
public:

    /** @brief Constructor: CalqEncoder
     *
     *  Initializes a new CalqEncoder instance with a SAM file to be encoded,
     *  a CQ file to write the bitstream to and a set of FASTA files to read
     *  the reference sequence(s) from.
     *
     *  @param samFileName The name of the SAM file to be encoded. The file has
     *         to exist and  has to be readable.
     *  @param cqFileName The name of the CQ file to write the encoded
     *         bitstream to. The file may not exist already.
     *  @param fastaFileNames A vector containing the names of all FASTA
     *         files which contents shall be used as reference sequence(s).
     *         The FASTA files have to exists and have to be readable.
     *  @return Void
     */
    CalqEncoder(const std::string &samFileName,
                const std::string &cqFileName,
                const std::vector<std::string> &fastaFileNames);
    ~CalqEncoder(void);

    void encode(void);

private:
    ofbitstream ofbs;
    QualEncoder qualEncoder;
    SAMParser samParser;
};

/** @brief Class: CalqDecoder
 *
 *  The CalqDecoder provides only one interface to the outside which is the
 *  member function decode; it decodes the CQ file with the name infileName
 *  and writes the decoded quality scores to the file with the name outfileName.
 */
class CalqDecoder: public CalqCodec {
public:
    CalqDecoder(const std::string &cqFileName,
                const std::string &qualFileName,
                const std::vector<std::string> &fastaFileNames);
    ~CalqDecoder(void);

    void decode(void);
 
private:
    ifbitstream ifbs;
    std::ofstream ofs;
    QualDecoder qualDecoder;
};

#endif // CALQCODEC_H


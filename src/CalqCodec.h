/** @file CalqCodec.h
 *  @brief This file contains the definitions of the CalqEncoder and
 *         CalqDecoder classes.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef CALQCODEC_H
#define CALQCODEC_H

#include "Common/File.h"
#include "Parsers/FASTAParser.h"
#include "Parsers/SAMParser.h"
#include "QualCodec.h"
#include <fstream>
#include <string>
#include <vector>

class CalqCodec {
public:
    CalqCodec(const std::string &inFileName,
              const std::string &outFileName);
    virtual ~CalqCodec(void);

protected:
    const std::string inFileName;
    const std::string outFileName;
};

class CalqEncoder: public CalqCodec {
public:
    CalqEncoder(const std::string &samFileName,
                const std::string &cqFileName,
                const std::vector<std::string> &fastaFileNames,
                const unsigned int &blockSize,
                const unsigned int &polyploidy,
                const bool &quantizedPrintout,
                const int &qvMin,
                const int &qvMax);
    ~CalqEncoder(void);

    void encode(void);

private:
    unsigned int blockSize;
    File cqFile;
    std::vector<FASTAReference> fastaReferences;
    unsigned int polyploidy;
    QualEncoder qualEncoder;
    SAMParser samParser;

    size_t writeFileHeader(void);
};

class CalqDecoder: public CalqCodec {
public:
    CalqDecoder(const std::string &cqFileName,
                const std::string &qualFileName,
                const std::string &samFileName);
    ~CalqDecoder(void);

    void decode(void);

private:
    File cqFile;
    File qualFile;
    QualDecoder qualDecoder;
    SAMParser samParser;

    size_t readFileHeader(void);
};

#endif // CALQCODEC_H


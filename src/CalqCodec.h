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

#include "Common/CLIOptions.h"


//#include "QualCodec.h"

#include "IO/File.h"
#include "IO/SAMFile.h"

#include <fstream>
#include <string>
#include <vector>

class CalqEncoder {
public:
    explicit CalqEncoder(const CLIOptions &cliOptions);
    ~CalqEncoder(void);

    void encode(void);

private:
    // From CLIOptions
    bool force;
    std::string inputFileName;
    std::string outputFileName;
    int blockSize;
    int polyploidy;
    int qualityValueMax;
    int qualityValueMin;
    std::vector<std::string> referenceFileNames;

    // Internal
    SAMFile samFile;
    File cqFile;
    //std:vector<FASTAReference> fastaReferences;

    //QualEncoder qualEncoder;


    size_t writeFileHeader(void);
};

class CalqDecoder {
public:
    explicit CalqDecoder(const CLIOptions &cliOptions);
    ~CalqDecoder(void);

    void decode(void);

private:
    File cqFile;
    File qualFile;
    //QualDecoder qualDecoder;
    //SAMParser samParser;

    size_t readFileHeader(void);
};

#endif // CALQCODEC_H


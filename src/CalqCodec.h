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
#include "IO/CQFile.h"
#include "IO/SAMFile.h"
//#include "QualCodec.h"

class CalqEncoder {
public:
    explicit CalqEncoder(const CLIOptions &cliOptions);
    ~CalqEncoder(void);

    void encode(void);

private:
    size_t blockSize;
    CQFile cqFile;
    bool force;
    unsigned int polyploidy;
    //QualEncoder qualEncoder;
    std::vector<std::string> referenceFileNames;
    SAMFile samFile;
    int qualityValueMax;
    int qualityValueMin;
};

class CalqDecoder {
public:
    explicit CalqDecoder(const CLIOptions &cliOptions);
    ~CalqDecoder(void);

    void decode(void);

private:
    CQFile cqFile;
    File qualFile;
    //QualDecoder qualDecoder;
    SAMFile sideInformationFile;
};

#endif // CALQCODEC_H


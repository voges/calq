#ifndef CALQCODEC_H
#define CALQCODEC_H

#include "bitstream.h"
#include "QualCodec.h"
#include "SAMParser.h"
#include <string>

class CalqEncoder {
public:
    CalqEncoder(std::string &infileName, std::string &outfileName, size_t &blockSize);
    ~CalqEncoder(void);

    void encode(void);

private:
    std::string infileName;
    std::string outfileName;

    SAMParser samParser;
    ofbitstream ofbs;

    size_t blockSize;

    QualEncoder qualEncoder;
};

class CalqDecoder {
public:
    CalqDecoder(std::string &infileName, std::string &outfileName);
    ~CalqDecoder(void);

    void decode(void);

private:
    std::string infileName;
    std::string outfileName;

    //QualDecoder qualDecoder;
};

class CalqInfoTool {
public:
    CalqInfoTool(std::string &infileName);
    ~CalqInfoTool(void);

    void extractInfo(void);

private:
    std::string infileName;
};

#endif // CALQCODEC_H


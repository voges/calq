#ifndef SAMPARSER_H
#define SAMPARSER_H

#include "SAMRecord.h"
#include <fstream>
#include <string>

class SAMParser {
public:
    SAMParser(const std::string &filename);
    ~SAMParser(void);

    std::string header; ///< SAM header
    SAMRecord curr;     ///< current SAM record

    void next(void);
    bool hasNext(void);

private:
    std::ifstream ifs;

    void parse(void);
    void readSAMHeader(void);
};

#endif // SAMPARSER_H


/** @file SAMParser.h
 *  @brief This file contains the definition of the SAMParser class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: what (who)
 */

#ifndef SAMPARSER_H
#define SAMPARSER_H

#include "Parsers/SAMRecord.h"
#include <fstream>
#include <string>

/** @brief Class: SAMParser
 *
 *  This class parses a SAM file. The constructor reads and stores the SAM
 *  header. The member function next the reads the next record and stores it in
 *  curr where it is accessible to the outside world.
 */
class SAMParser {
public:
    SAMParser(const std::string &filename);
    ~SAMParser(void);

    std::string header; // SAM header
    SAMRecord curr;     // current SAM record

    bool next(void);
private:
    std::ifstream ifs;

    void parse(void);
    void readSAMHeader(void);
};

#endif // SAMPARSER_H


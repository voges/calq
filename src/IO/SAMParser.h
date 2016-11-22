/** @file SAMParser.h
 *  @brief This file contains the definition of the SAMParser class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef SAMPARSER_H
#define SAMPARSER_H

#include "IO/SAMRecord.h"
#include <fstream>
#include <string>

/** @brief Class: SAMParser
 *
 *  This class parses a SAM file. The constructor reads and stores the SAM
 *  header and stores the first record in 'curr'. The member function 'next'
 *  reads the next record and stores it 'curr'.
 */
class SAMParser {
public:
    /** @brief Constructor: SAMParser
     *
     *  Initializes a new SAMParser instance with a SAM file with the name
     *  'fileName'. Checks if 'fileName' has the correct file name extension
     *  and checks if the file is accessible. Reads the SAM header by calling
     *  'readSAMheader'. The SAM header is stored in 'header'. Reads the first
     *  record and stores it in 'curr'.
     */
    explicit SAMParser(const std::string &fileName);

    /** @brief Destructor: SAMParser
     *
     *  Destructs a SAMParser instance.
     */
    ~SAMParser(void);

    /** @brief Member function: next
     *
     *  Reads the next record and stores it in 'curr'.
     *
     *  @return Returns true if a next record has been found and read and false
     *          otherwise.
     */
    bool next(void);

    std::string header; /// SAM header
    SAMRecord curr;     /// current SAM record
private:
    std::ifstream ifs;

    void parse(void);
    void readSAMHeader(void);
};

#endif // SAMPARSER_H


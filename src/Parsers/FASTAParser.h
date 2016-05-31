/** @file FASTAParser.h
 *  @brief This file contains the definition of the FASTAParser class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#ifndef FASTAPARSER_H
#define FASTAPARSER_H

#include "constants.h"
#include "Parsers/FASTAReference.h"
#include <string>
#include <vector>

/** @brief Class: FASTAParser
 *
 *  This class can parse FASTA files with the public member function
 * 'parseFile'.
 */
class FASTAParser {
public:
    FASTAParser(void);
    ~FASTAParser(void);

    void parseFile(const std::string &filename, std::vector<FASTAReference> &fastaReferences);

private:
    char line[4 * KB]; // usually, lines are limited to 80 chars, so 4 KB
                       // should be enough space
};

#endif // FASTAPARSER_H


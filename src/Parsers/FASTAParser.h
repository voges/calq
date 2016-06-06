/** @file FASTAParser.h
 *  @brief This file contains the definition of the FASTAParser class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: what (who)
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
    /** @brief Constructor: FASTAParser
     *
     *  Initializes a new FASTAParser instance.
     */
    FASTAParser(void);

    /** @brief Destructor: FASTAParser
     *
     *  Destructs a FASTAParser instance.
     */
    ~FASTAParser(void);

    /** @brief Member function: parseFile
     *
     *  Parses a FASTA file with name 'filename' and appends the found
     *  sequences and the associated headers to the vector 'fastaReferences'.
     *  Checks if 'filename' has the correct filename extension and also
     *  checks if the file is accessible.
     *
     *  @param filename FASTA file name
     *  @param fastaReferences A vector containing elements of type
     *         FASTAReference.
     *  @return Void.
     */
    void parseFile(const std::string &filename, std::vector<FASTAReference> &fastaReferences);
private:
    char line[4 * KB]; // usually, lines are limited to 80 chars, so 4 KB should be enough
};

#endif // FASTAPARSER_H


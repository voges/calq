/** @file FASTAParser.h
 *  @brief This file contains the definition of the FASTAParser class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef FASTAPARSER_H
#define FASTAPARSER_H

#include "Common/constants.h"
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

    /** @brief Copy constructor: FASTAParser
     *
     *  Initializes a new FASTAParser instance as a copy.
     */
    FASTAParser(const FASTAParser &fastaParser);

    /** @brief Destructor: FASTAParser
     *
     *  Destructs a FASTAParser instance.
     */
    ~FASTAParser(void);

    /** @brief Member function: parseFile
     *
     *  Parses a FASTA file with name 'fileName' and appends the found
     *  sequences and the associated headers to the vector 'fastaReferences'.
     *
     *  @param fileName FASTA file name
     *  @param fastaReferences A vector containing elements of type
     *         FASTAReference.
     *  @return Void.
     */
    void parseFile(const std::string &fileName, std::vector<FASTAReference> &fastaReferences);
private:
    char *line;
    size_t lineSize;
};

#endif // FASTAPARSER_H


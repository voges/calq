/** @file FASTAParser.cc
 *  @brief This file contains the implementation of the FASTAParser class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: what (who)
 */

#include "Parsers/FASTAParser.h"
#include "Common/Exceptions.h"
#include "Common/fileSystemHelpers.h"
#include <fstream>
#include <iostream>

FASTAParser::FASTAParser(void)
{
    // empty
}

FASTAParser::~FASTAParser(void) 
{
    // empty
}

void FASTAParser::parseFile(const std::string &fileName, std::vector<FASTAReference> &fastaReferences)
{
    if (   fileNameExtension(fileName) != std::string("fa")
        && fileNameExtension(fileName) != std::string("fasta")) {
        throwErrorException("Reference file extension must be 'fa' or 'fasta'");
    }

    if (!fileExists(fileName)) {
        throwErrorException("Cannot access reference file");
    }

    std::ifstream ifs;
    ifs.open(fileName.c_str(), std::ios::in);

    FASTAReference fastaReference;
    bool first = true;

    while (ifs.getline(line, sizeof(line))) {
        if (line[0] == '>') {
            if (first) {
                fastaReference.header = line+1; // do not take the '>'
                first = false;
            } else {
                fastaReferences.push_back(fastaReference);
                fastaReference.header = line+1; // do not take the '>'
                fastaReference.sequence = "";
            }
            // remove everything after the first space from header line
            std::size_t foundSpace = fastaReference.header.find_first_of(" ");
            if (foundSpace != std::string::npos) {
                fastaReference.header = fastaReference.header.substr(0, foundSpace);
            }
        } else {
            fastaReference.sequence += line;
        }
    }

    fastaReferences.push_back(fastaReference);

    ifs.close();
}


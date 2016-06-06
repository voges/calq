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
#include "common.h"
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

void FASTAParser::parseFile(const std::string &filename, 
                            std::vector<FASTAReference> &fastaReferences)
{
    if (   filenameExtension(filename) != std::string("fa")
        && filenameExtension(filename) != std::string("fasta")) {
        throwErrorException("Reference file extension must be 'fa' or 'fasta'");
    }

    if (!fileExists(filename)) {
        throwErrorException("Cannot access reference file");
    }

    std::ifstream ifs;
    ifs.open(filename.c_str(), std::ios::in);

    FASTAReference fastaReference;
    bool first = true;

    while (ifs.getline(line, sizeof(line))) {
        if (line[0] == '>') {
            if (first) {
                fastaReference.header = line;
                first = false;
            } else {
                fastaReferences.push_back(fastaReference);
                fastaReference.header = line;
                fastaReference.sequence = "";
            }
        } else {
            fastaReference.sequence += line;
        }
    }

    fastaReferences.push_back(fastaReference);

    ifs.close();
}


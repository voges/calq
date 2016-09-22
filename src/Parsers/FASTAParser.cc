/** @file FASTAParser.cc
 *  @brief This file contains the implementation of the FASTAParser class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "FASTAParser.h"
#include "Common/Exceptions.h"
#include "Common/helpers.h"
#include <fstream>
#include <iostream>

FASTAParser::FASTAParser(void)
    : line(NULL)
    , lineSize(sizeof(char) * (4*KB))
{
    line = (char *)malloc(lineSize); // usually lines are limited to 80 chars, so 4 KB should be enough
}

FASTAParser::FASTAParser(const FASTAParser &fastaParser)
    : line(fastaParser.line)
    , lineSize(fastaParser.lineSize)
{
    // empty
}

FASTAParser::~FASTAParser(void) 
{
    free(line);
}

void FASTAParser::parseFile(const std::string &fileName, std::vector<FASTAReference> &fastaReferences)
{
    if (fileName.empty() == true) {
        throwErrorException("No file name given");
    }
    if (   fileNameExtension(fileName) != std::string("fa")
        && fileNameExtension(fileName) != std::string("fasta")) {
        throwErrorException("Reference file extension must be 'fa' or 'fasta'");
    }
    if (fileExists(fileName) == false) {
        throwErrorException("Cannot access reference file");
    }
    if (fastaReferences.size() != 0) {
        throwErrorException("Vector 'fastaReferences' is not empty");
    }

    std::ifstream ifs;
    ifs.open(fileName.c_str(), std::ios::in);

    FASTAReference fastaReference;
    bool first = true;

    while (ifs.getline(line, lineSize)) {
        if (line[0] == '>') {
            if (first) {
                fastaReference.header = line+1; // do not take the '>'
                first = false;
            } else {
                fastaReferences.push_back(fastaReference);
                fastaReference.header = line+1; // do not take the '>'
                fastaReference.sequence = "";
            }
            // Remove everything after the first space from header line
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


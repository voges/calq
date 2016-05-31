/** @file FASTAParser.cc
 *  @brief This file contains the implementation of the FASTAParser class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#include "FASTAParser.h"
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

void FASTAParser::parseFile(const std::string &filename, std::vector<FASTAReference> &fastaReferences)
{
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


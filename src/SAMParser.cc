/** @file SAMParser.cc
 *  @brief This file contains the implementation of the SAMParser class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#include "SAMParser.h"
#include "Exceptions.h"
#include <iostream>
#include <stdlib.h>
#include <string.h>

SAMParser::SAMParser(const std::string &filename)
    : header()
    , ifs()
{
    ifs.open(filename.c_str(), std::ios::in);
    readSAMHeader();
}

SAMParser::~SAMParser (void) 
{
    ifs.close();
}

void SAMParser::readSAMHeader(void)
{
    bool foundHeader = false;

    while (ifs.getline(curr.line, sizeof(curr.line))) {
        if (*(curr.line) == '@') {
            header += curr.line;
            header += "\n";
            foundHeader = true;
        } else {
            if (foundHeader == false) {
                throwErrorException("SAM header is missing");
            }
            size_t offset = -strlen(curr.line);
            ifs.seekg(offset, std::ios::cur);
            break;
        }
    }
}

void SAMParser::parse(void)
{
    size_t l = strlen(curr.line) - 1;

    while (l && (curr.line[l] == '\r' || curr.line[l] == '\n')) {
        curr.line[l--] = '\0';
    }

    char *c = curr.qname = curr.line;
    int f = 1;

    while (*c) {
        if (*c == '\t') {
            if (f ==  1) curr.flag  = (uint16_t)atoi(c + 1);
            if (f ==  2) curr.rname = c + 1;
            if (f ==  3) curr.pos   = (uint32_t)atoi(c + 1);
            if (f ==  4) curr.mapq  = (uint8_t)atoi(c + 1);
            if (f ==  5) curr.cigar = c + 1;
            if (f ==  6) curr.rnext = c + 1;
            if (f ==  7) curr.pnext = (uint32_t)atoi(c + 1);
            if (f ==  8) curr.tlen  = (int64_t)atoi(c + 1);
            if (f ==  9) curr.seq   = c + 1;
            if (f == 10) curr.qual  = c + 1;
            if (f == 11) curr.opt   = c + 1;
            f++;
            *c = '\0';
            if (f == 12) break;
        }
        c++;
    }

    if (f == 11) curr.opt = c;
}

bool SAMParser::hasNext(void)
{
    std::streampos current = ifs.tellg();
    ifs.seekg(0, ifs.end);
    std::streampos end = ifs.tellg();
    bool hasNext = current == end ? false : true;
    ifs.seekg(current, ifs.beg); //restore stream position
    return hasNext;
}

void SAMParser::next(void)
{
    if (ifs.getline(curr.line, sizeof(curr.line))) {
        parse();
    } else {
        throwErrorException("Tried to read record behind EOF");
    }
}


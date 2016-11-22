/** @file SAMFile.cc
 *  @brief This file contains the implementation of the SAMFile class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "SAMFile.h"
#include "Common/constants.h"
#include "Common/Exceptions.h"
#include <stdio.h>

SAMFile::SAMFile(const std::string &path, const char *mode, const size_t &blockSize)
    : File(path, mode)
    , blockSize(blockSize)
    , line(NULL)
    , lineSize(sizeof(char) * (1 * MB))
{
    // For fgets and (f)seek to work together properly, a SAMFile instance
    // has to be opened in binary mode
    if (strcmp(mode, "rb") != 0) {
        throwErrorException("rb mode required");
    }

    // 1 million chars should be enough
    line = (char *)malloc(lineSize);

    bool foundHeader = false;

    // Read SAM header
    size_t alignmentSectionBegin = tell();
    for (;;) {
        alignmentSectionBegin = tell();
        if (fgets(line, lineSize, fp) != NULL) {
            if (line[0] == '@') {
                header += line;
            } else {
                
                break;
            }
        } else {
            throwErrorException("Could not read SAM header");
        }
    }
    seek(alignmentSectionBegin);
    if (header.empty() == true) {
        throwErrorException("SAM header is missing");
    }
    LOG("SAM header: %s", header.c_str());

    readBlock();
}

SAMFile::~SAMFile(void)
{
    free(line);
}

void SAMFile::readBlock(void) {
    for (size_t i = 0; i < blockSize; i++) {
        if (fgets(line, lineSize, fp) != NULL) {
            std::string lineString = line;
            LOG("line %d: %s", i, lineString.c_str());
        } else {
            LOG("Truncated block");
            break;
        }
    }
}


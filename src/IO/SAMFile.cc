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
#include "IO/SAMRecord.h"
#include <string.h>

static void parseLine(char *fields[SAMRecord::NUM_FIELDS], char *line)
{
    char *c = line;
    char *pc = c;
    int f = 0;

    while (*c) {
        if (*c == '\t') {
            *c = '\0';
            if (f ==  0) fields[f] = pc;
            if (f ==  1) fields[f] = pc;
            if (f ==  2) fields[f] = pc;
            if (f ==  3) fields[f] = pc;
            if (f ==  4) fields[f] = pc;
            if (f ==  5) fields[f] = pc;
            if (f ==  6) fields[f] = pc;
            if (f ==  7) fields[f] = pc;
            if (f ==  8) fields[f] = pc;
            if (f ==  9) fields[f] = pc;
            if (f == 10) fields[f] = pc;
            if (f == 11) fields[f] = pc;
            f++;
            if (f == 12) { break; }
            pc = c + 1;
        }
        c++;
    }

    if (f == 11) { fields[f] = pc; }
}

SAMFile::SAMFile(const std::string &path,
                 const SAMFile::Mode &mode)
    : File(path, mode)
    , block()
    , header("")
    , line(NULL)
    , lineSize(sizeof(char) * (1*MB))
{
    // Check arguments
    if (path.empty() == true) {
        throwErrorException("path is empty");
    }
    if (mode != SAMFile::MODE_READ) {
        throwErrorException("Currently only MODE_READ supported");
    }

    // 1 million chars should be enough
    line = (char *)malloc(lineSize);

    // Read SAM header
    size_t fpos = tell();
    for (;;) {
        fpos = tell();
        if (fgets(line, lineSize, fp) != NULL) {
            // Trim line
            size_t l = strlen(line) - 1;
            while (l && (line[l] == '\r' || line[l] == '\n')) { line[l--] = '\0'; }

            if (line[0] == '@') {
                header += line;
                header += "\n";
            } else {
                break;
            }
        } else {
            throwErrorException("Could not read SAM header");
        }
    }
    seek(fpos); // rewind to the begin of the alignment section
    if (header.empty() == true) {
        LOG("No SAM header found");
    }
}

SAMFile::~SAMFile(void)
{
    free(line);
}

size_t SAMFile::readBlock(const size_t &blockSize)
{
    if (blockSize < 1) {
        throwErrorException("blockSize must be greater than zero");
    }

    block.clear();

    size_t lineCnt = 0;
    for (size_t i = 0; i < blockSize; i++) {
        if (fgets(line, lineSize, fp) != NULL) {
            // Trim line
            size_t l = strlen(line) - 1;
            while (l && (line[l] == '\r' || line[l] == '\n')) { line[l--] = '\0'; }

            // Parse line, construct SAMRecord and push it to the current block
            char *fields[SAMRecord::NUM_FIELDS];
            parseLine(fields, line);
            SAMRecord samRecord(fields);
            block.push_back(samRecord);

            lineCnt++;
        } else {
            LOG("Truncated block - only read %d lines (%d were requested)", lineCnt, blockSize);
            break;
        }
    }

    return lineCnt;
}


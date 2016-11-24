/** @file SAMFile.cc
 *  @brief This file contains the implementation of the SAMFile class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "IO/SAMFile.h"
#include "Common/constants.h"
#include "Common/Exceptions.h"
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
    , currentBlock()
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

void SAMFile::readBlock(const size_t &blockSize)
{
    if (blockSize < 1) {
        throwErrorException("blockSize must be greater than zero");
    }

    currentBlock.reset();

    std::string rnamePrev("");
    uint32_t posPrev = 0;

    for (size_t i = 0; i < blockSize; i++) {
        size_t fpos = tell();
        if (fgets(line, lineSize, fp) != NULL) {
            // Trim line
            size_t l = strlen(line) - 1;
            while (l && (line[l] == '\r' || line[l] == '\n')) { line[l--] = '\0'; }

            // Parse line and construct samRecord
            char *fields[SAMRecord::NUM_FIELDS];
            parseLine(fields, line);
            SAMRecord samRecord(fields);

            if (samRecord.isMapped() == true) {
                if (rnamePrev.empty() == true) {
                    // This is the first mapped record in this block; just store
                    // its RNAME and POS and add it to the current block
                    rnamePrev = samRecord.rname;
                    posPrev = samRecord.pos;
                    currentBlock.records.push_back(samRecord);
                    currentBlock.numMappedRecords++;
                } else {
                    // We already have a mapped record in this block
                    if (rnamePrev == samRecord.rname) {
                        // RNAME didn't change, check POS
                        if (samRecord.pos >= posPrev) {
                            // Everything fits, just update posPrev and push
                            // the samRecord to the current block
                            posPrev = samRecord.pos;
                            currentBlock.records.push_back(samRecord);
                            currentBlock.numMappedRecords++;
                        } else {
                            throwErrorException("SAM file is not sorted");
                        }
                    } else {
                        // RNAME changed, seek back and break
                        seek(fpos);
                        LOG("RNAME changed - read only %zu lines (%zu were requested)", currentBlock.numRecords, blockSize);
                        break;
                    }
                }
            } else {
                currentBlock.records.push_back(samRecord);
                currentBlock.numUnmappedRecords++;
            }
        } else {
            LOG("Truncated block - read only %zu lines (%zu were requested)", currentBlock.numRecords, blockSize);
            break;
        }
        currentBlock.numRecords++;
    }
}


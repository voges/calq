/** @file SAMFile.cc
 *  @brief This file contains the implementation of the SAMFile class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "IO/SAM/SAMFile.h"
#include "Common/Exceptions.h"

static void parseLine(char *fields[cq::SAMRecord::NUM_FIELDS], char *line)
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

cq::SAMFile::SAMFile(const std::string &path, const Mode &mode)
    : File(path, mode)
    , currentBlock()
    , header("")
    , m_line(NULL)
    , m_numBlocksRead(0)
    , m_numMappedRecordsRead(0)
    , m_numUnmappedRecordsRead(0)
{
    // Check arguments
    if (path.empty() == true) {
        throwErrorException("path is empty");
    }
    if (mode != MODE_READ) {
        throwErrorException("Currently only MODE_READ supported");
    }

    // 1 million chars should be enough
    m_line = (char *)malloc(LINE_SIZE);

    // Read SAM header
    size_t fpos = tell();
    for (;;) {
        fpos = tell();
        if (fgets(m_line, LINE_SIZE, m_fp) != NULL) {
            // Trim line
            size_t l = strlen(m_line) - 1;
            while (l && (m_line[l] == '\r' || m_line[l] == '\n')) {
                m_line[l--] = '\0';
            }

            if (m_line[0] == '@') {
                header += m_line;
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
        CQ_LOG("No SAM header found");
    }
}

cq::SAMFile::~SAMFile(void)
{
    free(m_line);
}

size_t cq::SAMFile::numBlocksRead(void) const
{
    return m_numBlocksRead;
}

size_t cq::SAMFile::numMappedRecordsRead() const
{
    return m_numMappedRecordsRead;
}

size_t cq::SAMFile::numUnmappedRecordsRead() const
{
    return m_numUnmappedRecordsRead;
}

size_t cq::SAMFile::numRecordsRead() const
{
    return (m_numMappedRecordsRead + m_numUnmappedRecordsRead);
}

size_t cq::SAMFile::readBlock(const size_t &blockSize)
{
    if (blockSize < 1) {
        throwErrorException("blockSize must be greater than zero");
    }

    currentBlock.reset();

    std::string rnamePrev("");
    uint32_t posPrev = 0;

    for (size_t i = 0; i < blockSize; i++) {
        size_t fpos = tell();
        if (fgets(m_line, LINE_SIZE, m_fp) != NULL) {
            // Trim line
            size_t l = strlen(m_line) - 1;
            while (l && (m_line[l] == '\r' || m_line[l] == '\n')) {
                m_line[l--] = '\0';
            }

            // Parse line and construct samRecord
            char *fields[SAMRecord::NUM_FIELDS];
            parseLine(fields, m_line);
            SAMRecord samRecord(fields);

            if (samRecord.isMapped() == true) {
                if (rnamePrev.empty() == true) {
                    // This is the first mapped record in this block; just store
                    // its RNAME and POS and add it to the current block
                    rnamePrev = samRecord.rname;
                    posPrev = samRecord.pos;
                    currentBlock.records.push_back(samRecord);
                    currentBlock.m_numMappedRecords++;
                } else {
                    // We already have a mapped record in this block
                    if (rnamePrev == samRecord.rname) {
                        // RNAME didn't change, check POS
                        if (samRecord.pos >= posPrev) {
                            // Everything fits, just update posPrev and push
                            // the samRecord to the current block
                            posPrev = samRecord.pos;
                            currentBlock.records.push_back(samRecord);
                            currentBlock.m_numMappedRecords++;
                        } else {
                            throwErrorException("SAM file is not sorted");
                        }
                    } else {
                        // RNAME changed, seek back and break
                        seek(fpos);
                        CQ_LOG("RNAME changed - read only %zu record(s) (%zu requested)", currentBlock.numRecords(), blockSize);
                        break;
                    }
                }
            } else {
                currentBlock.records.push_back(samRecord);
                currentBlock.m_numUnmappedRecords++;
            }
        } else {
            CQ_LOG("Truncated block - read only %zu record(s) (%zu requested) - reached EOF", currentBlock.numRecords(), blockSize);
            break;
        }
    }

    if (currentBlock.numRecords() > 0) {
        m_numBlocksRead++;
        m_numMappedRecordsRead += currentBlock.numMappedRecords();
        m_numUnmappedRecordsRead += currentBlock.numUnmappedRecords();
    }
    
    return currentBlock.numRecords();
}


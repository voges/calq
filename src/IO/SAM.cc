/** @file SAM.cc
 *  @brief This file contains the implementation of the classes defined in
 *         SAM.h.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "IO/SAM.h"
#include "Common/constants.h"
#include "Common/Exceptions.h"
#include <string.h>

cq::SAMRecord::SAMRecord(char *fields[NUM_FIELDS])
    : qname(fields[0])
    , flag((uint16_t)atoi(fields[1]))
    , rname(fields[2])
    , pos((uint32_t)atoi(fields[3]))
    , mapq((uint8_t)atoi(fields[4]))
    , cigar(fields[5])
    , rnext(fields[6])
    , pnext((uint32_t)atoi(fields[7]))
    , tlen((int64_t)atoi(fields[8]))
    , seq(fields[9])
    , qual(fields[10])
    , opt(fields[11])
    , posMin(0)
    , posMax(0)
    , m_mapped(false)
{
    check();

    if (m_mapped == true) {
        // Compute 0-based first position and 0-based last position this record
        // is mapped to on the reference used for alignment
        posMin = pos - 1;
        posMax = pos - 1;

        size_t cigarIdx = 0;
        size_t cigarLen = cigar.length();
        uint32_t opLen = 0; // length of current CIGAR operation

        for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++) {
            if (isdigit(cigar[cigarIdx])) {
                opLen = opLen*10 + (uint32_t)cigar[cigarIdx] - (uint32_t)'0';
                continue;
            }
            switch (cigar[cigarIdx]) {
            case 'M':
            case '=':
            case 'X':
                posMax += opLen;
                break;
            case 'I':
            case 'S':
                break;
            case 'D':
            case 'N':
                posMax += opLen;
                break;
            case 'H':
            case 'P':
                break; // these have been clipped
            default:
                throwErrorException("Bad CIGAR string");
            }
            opLen = 0;
        }

        posMax -= 1;
    }
}

cq::SAMRecord::~SAMRecord(void)
{
    // empty
}

bool cq::SAMRecord::isMapped(void) const
{
    return m_mapped;
}

void cq::SAMRecord::print(void) const
{
    printf("%s\t", qname.c_str());
    printf("%d\t", flag);
    printf("%s\t", rname.c_str());
    printf("%d\t", pos);
    printf("%d\t", mapq);
    printf("%s\t", cigar.c_str());
    printf("%s\t", rnext.c_str());
    printf("%d\t", pnext);
    printf("%ld\t", tlen);
    printf("%s\t", seq.c_str());
    printf("%s\t", qual.c_str());
    printf("%s\t", opt.c_str());
    printf("\n");
    printf("isMapped: %d, ", m_mapped);
    printf("posMin: %d, ", posMin);
    printf("posMax: %d\n", posMax);
}

void cq::SAMRecord::check(void)
{
    // Check all fields
    if (qname.empty() == true) { throwErrorException("qname is empty"); }
    // flag
    if (rname.empty() == true) { throwErrorException("rname is empty"); }
    // pos
    // mapq
    if (cigar.empty() == true) { throwErrorException("cigar is empty"); }
    if (rnext.empty() == true) { throwErrorException("rnext is empty"); }
    // pnext
    // tlen
    if (seq.empty() == true) { throwErrorException("seq is empty"); }
    if (qual.empty() == true) { throwErrorException("qual is empty"); }
    if (opt.empty() == true) { LOG("opt is empty"); }

    // Check if this record is mapped
    if ((flag & 0x4) != 0) {
        m_mapped = false;
    } else {
        m_mapped = true;
        if (rname == "*" || pos == 0 || cigar == "*" || seq == "*" || qual == "*") {
            throwErrorException("Corrupted record");
        }
    }
}

cq::SAMBlock::SAMBlock(void)
    : records()
    , m_numMappedRecords(0)
    , m_numUnmappedRecords(0)
    , m_numRecords(0)
{
    // empty
}

cq::SAMBlock::~SAMBlock(void)
{
    // empty
}

size_t cq::SAMBlock::numMappedRecords(void) const
{
    return m_numMappedRecords;
}

size_t cq::SAMBlock::numUnmappedRecords(void) const
{
    return m_numUnmappedRecords;
}

size_t cq::SAMBlock::numRecords(void) const
{
    return m_numRecords;
}

void cq::SAMBlock::reset(void)
{
    records.clear();
    m_numUnmappedRecords = 0;
    m_numMappedRecords = 0;
    m_numRecords = 0;
}

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
        LOG("No SAM header found");
    }
}

cq::SAMFile::~SAMFile(void)
{
    free(m_line);
}

void cq::SAMFile::readBlock(const size_t &blockSize)
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
                        LOG("RNAME changed - read only %zu lines (%zu were requested)", currentBlock.numRecords(), blockSize);
                        break;
                    }
                }
            } else {
                currentBlock.records.push_back(samRecord);
                currentBlock.m_numUnmappedRecords++;
            }
        } else {
            LOG("Truncated block - read only %zu lines (%zu were requested)", currentBlock.numRecords(), blockSize);
            break;
        }
        currentBlock.m_numRecords++;
    }
}


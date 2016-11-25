/** @file SAMRecord.cc
 *  @brief This file contains the implementation of the SAMRecord class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "IO/SAM/SAMRecord.h"
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
                opLen = opLen * 10 + (uint32_t)cigar[cigarIdx] - (uint32_t)'0';
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

void cq::SAMRecord::extractPileup(const uint32_t &pileupPosMin,
                                  const uint32_t &pileupPosMax,
                                  std::deque<std::string> &seqPileup,
                                  std::deque<std::string> &qualPileup) const
{
    if ((seqPileup.size() == 0) || (qualPileup.size() == 0)) {
        throwErrorException("Pileups are empty");
    }
    if (seqPileup.size() != qualPileup.size()) {
        throwErrorException("Pileup sizes differ");
    }
    if ((pileupPosMin > posMin) || (pileupPosMax < posMax)) {
        throwErrorException("Pileup does not overlap record");
    }

    size_t cigarIdx = 0;
    size_t cigarLen = cigar.length();
    size_t opLen = 0; // length of current CIGAR operation
    size_t idx = 0;
    size_t pileupIdx = posMin - pileupPosMin;

    for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++) {
        if (isdigit(cigar[cigarIdx])) {
            opLen = opLen*10 + (size_t)cigar[cigarIdx] - (size_t)'0';
            continue;
        }

        switch (cigar[cigarIdx]) {
        case 'M':
        case '=':
        case 'X':
            for (size_t i = 0; i < opLen; i++) {
                seqPileup[pileupIdx] += seq[idx];
                qualPileup[pileupIdx] += qual[idx];
                idx++; pileupIdx++;
            }
            break;
        case 'I':
        case 'S':
            idx += opLen;
            break;
        case 'D':
        case 'N':
            pileupIdx += opLen;
            break;
        case 'H':
        case 'P':
            break; // these have been clipped
        default:
            throwErrorException("Bad CIGAR string");
        }

        opLen = 0;
    }
}

bool cq::SAMRecord::isMapped(void) const
{
    return m_mapped;
}

void cq::SAMRecord::printLong(void) const
{
    printShort();
    printf("isMapped: %d, ", m_mapped);
    printf("posMin: %d, ", posMin);
    printf("posMax: %d\n", posMax);
}

void cq::SAMRecord::printShort(void) const
{
    printf("%s\t", qname.c_str());
    printf("%d\t", flag);
    printf("%s\t", rname.c_str());
    printf("%d\t", pos);
    printf("%d\t", mapq);
    printf("%s\t", cigar.c_str());
    printf("%s\t", rnext.c_str());
    printf("%d\t", pnext);
    printf("%" PRId64 "\n", tlen);
    printf("%s\t", seq.c_str());
    printf("%s\t", qual.c_str());
    printf("%s\t", opt.c_str());
    printf("\n");
}

void cq::SAMRecord::printSeqWithPositionOffset(void) const
{
    printf("%6d-%6d|", posMin, posMax);
    for (unsigned int i = 0; i < posMin; i++) { printf(" "); }
    printf("%s\n", seq.c_str());
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
    if (opt.empty() == true) { CQ_LOG("opt is empty"); }

    // Check if this record is mapped
    if ((flag & 0x4) != 0) {
        m_mapped = false;
    }
    else {
        m_mapped = true;
        if (rname == "*" || pos == 0 || cigar == "*" || seq == "*" || qual == "*") {
            throwErrorException("Corrupted record");
        }
    }
}


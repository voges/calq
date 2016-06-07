/** @file MappedRecord.cc
 *  @brief This file contains the implementation of the MappedRecord class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: what (who)
 */

#include "MappedRecord.h"
#include "Exceptions.h"

MappedRecord::MappedRecord(const uint32_t &pos,
                           const std::string &cigar,
                           const std::string &seq,
                           const std::string &qual)
    : firstPos(pos)
    , lastPos(pos)
    , alignedNucleotides("")
    , insertedNucleotides("")
    , alignedQualityValues("")
    , insertedQualityValues("")
    , cigar(cigar)
{
    size_t cigarIdx = 0;
    size_t cigarLen = cigar.length();
    size_t opLen = 0; // length of current CIGAR operation
    size_t idx = 0;

    for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++) {
        if (isdigit(cigar[cigarIdx])) {
            opLen = opLen * 10 + (size_t)cigar[cigarIdx] - (size_t)'0';
            continue;
        }

        size_t i = 0;
        switch (cigar[cigarIdx]) {
        case 'M':
        case '=':
        case 'X':
            // add matching parts
            alignedNucleotides.append(seq, idx, opLen);
            alignedQualityValues.append(qual, idx, opLen);
            idx += opLen;
            lastPos += opLen;
            break;
        case 'I':
        case 'S':
            // add inserted bases and quality values
            insertedNucleotides.append(seq, idx, opLen);
            insertedQualityValues.append(qual, idx, opLen);
            idx += opLen; // skip inserted part
            break;
        case 'D':
        case 'N': {
            // inflate sequence and quality scores
            for (i = 0; i < opLen; i++) { 
                alignedNucleotides += "_";
                alignedQualityValues += "_";
            }
            lastPos += opLen;
            break;
        }
        case 'H':
        case 'P':
            break; // these have been clipped
        default: 
            throwErrorException("Bad CIGAR string");
        }

        opLen = 0;
    }
}

MappedRecord::~MappedRecord(void)
{
    // empty
}


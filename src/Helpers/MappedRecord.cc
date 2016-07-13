/** @file MappedRecord.cc
 *  @brief This file contains the implementation of the MappedRecord class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: what (who)
 */

#include "Helpers/MappedRecord.h"
#include "Exceptions.h"

MappedRecord::MappedRecord(const SAMRecord &samRecord)
    : posMin(samRecord.pos - 1) // SAM format counts from 1
    , posMax(samRecord.pos - 1)  // SAM format counts from 1
    , nucleotides(samRecord.seq)
    , qualityValues(samRecord.qual)
    , cigar(samRecord.cigar)
{
    // compute last position on the reference this record is mapping to
    size_t cigarIdx = 0;
    size_t cigarLen = cigar.length();
    size_t opLen = 0; // length of current CIGAR operation

    for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++) {
        if (isdigit(cigar[cigarIdx])) {
            opLen = opLen * 10 + (size_t)cigar[cigarIdx] - (size_t)'0';
            continue;
        }
        switch (cigar[cigarIdx]) {
        case 'M':
        case '=':
        case 'X':
            posMax += (uint32_t)opLen;
            break;
        case 'I':
        case 'S':
            break;
        case 'D':
        case 'N':
            posMax += (uint32_t)opLen;
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

MappedRecord::~MappedRecord(void)
{
    // empty
}

void MappedRecord::extractObservations(const uint32_t &observedPosMin,
                                       std::vector<std::string> &observedNucleotides,
                                       std::vector<std::string> &observedQualityValues)
{
    size_t cigarIdx = 0;
    size_t cigarLen = cigar.length();
    size_t opLen = 0; // length of current CIGAR operation
    size_t idx = 0;
    size_t observedIdx = posMin - observedPosMin;

    for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++) {
        if (isdigit(cigar[cigarIdx])) {
            opLen = opLen * 10 + (size_t)cigar[cigarIdx] - (size_t)'0';
            continue;
        }

        switch (cigar[cigarIdx]) {
        case 'M':
        case '=':
        case 'X':
            // add matching parts
            for (size_t i = 0; i < opLen; i++) {
                observedNucleotides[observedIdx] += nucleotides[idx];
                observedQualityValues[observedIdx] += qualityValues[idx];
                idx++;
                observedIdx++;
            }
            break;
        case 'I':
        case 'S':
            idx += opLen;
            break;
        case 'D':
        case 'N':
            observedIdx += opLen;
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


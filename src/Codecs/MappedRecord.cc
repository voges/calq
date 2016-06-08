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

MappedRecord::MappedRecord(const SAMRecord &samRecord, const uint32_t &positionOffset)
    : firstPosition(samRecord.pos - 1) // SAM format counts from 1
    , lastPosition(samRecord.pos - 1)  // SAM format counts from 1
    , positionOffset(positionOffset)
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
            lastPosition += (uint32_t)opLen;
            break;
        case 'I':
        case 'S':
            break;
        case 'D':
        case 'N':
            lastPosition += (uint32_t)opLen;
            break;
        case 'H':
        case 'P':
            break; // these have been clipped
        default: 
            throwErrorException("Bad CIGAR string");
        }
        opLen = 0;
    }

    lastPosition -= 1;
    //std::cout << "Last mapping position: " << lastPosition << std::endl;
}

MappedRecord::~MappedRecord(void)
{
    // empty
}

void MappedRecord::extractObservations(std::vector<std::string> &observedNucleotides,
                                       std::vector<std::string> &observedQualityValues)
{
    size_t cigarIdx = 0;
    size_t cigarLen = cigar.length();
    size_t opLen = 0; // length of current CIGAR operation

    size_t idx = 0;
    size_t observedIdx = firstPosition - positionOffset;

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


/** @file SAMRecord.cc
 *  @brief This file contains the implementation of the SAMRecord class.
 */

// Copyright 2015-2017 Leibniz Universitaet Hannover

#include "IO/SAM/SAMRecord.h"
#include "IO/FASTA/FASTAFile.h"

#include <string.h>
#include <queue>
#include <limits>

#include "Common/Exceptions.h"
#include "Common/log.h"

namespace calq {

SAMRecord::SAMRecord(char *fields[NUM_FIELDS])
    : qname(fields[0]),
      flag((uint16_t)atoi(fields[1])),
      rname(fields[2]),
      pos((uint32_t)atoi(fields[3])),
      mapq((uint8_t)atoi(fields[4])),
      cigar(fields[5]),
      rnext(fields[6]),
      pnext((uint32_t)atoi(fields[7])),
      tlen((int64_t)atoi(fields[8])),
      seq(fields[9]),
      qual(fields[10]),
      opt(fields[11]),
      posMin(0),
      posMax(0),
      mapped_(false) {
    check();

    if (mapped_ == true) {
        // Compute 0-based first position and 0-based last position this record
        // is mapped to on the reference used for alignment
        posMin = pos - 1;
        posMax = pos - 1;

        size_t cigarIdx = 0;
        size_t cigarLen = cigar.length();
        uint32_t opLen = 0;  // length of current CIGAR operation

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
                break;  // these have been clipped
            default:
                throwErrorException("Bad CIGAR string");
            }
            opLen = 0;
        }

        posMax -= 1;
    }
}

SAMRecord::~SAMRecord(void) {}

size_t SAMRecord::calcIndelScore(const FASTAFile& f, size_t offsetRef, size_t offsetRead, size_t abortScore) const{
    size_t score = 0; // Current mismatch score

    size_t cigarIdx = 0; //Current CIGAR position
    size_t cigarLen = cigar.length();
    size_t opLen = 0;  // length of current CIGAR operation
    size_t idx = 0; // Current base position relative to read start
    size_t pileupIdx = posMin; // current base position relatie to reference start

    //Buffers to deal with delay in case of different position offsets
    std::queue<char> seqBuffer;
    std::queue<char> qualBuffer;
    std::queue<char> refBuffer;

    //Loop inspired by addToPileupQueue
    for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++) {

        //Process numbers
        if (isdigit(cigar[cigarIdx])) {
            opLen = opLen*10 + (size_t)cigar[cigarIdx] - (size_t)'0';
            continue;
        }

        switch (cigar[cigarIdx]) {
        case 'M':
        case '=':
        case 'X':
            for (size_t i = 0; i < opLen; i++) {
                //Store in buffer, if offset passed
                if((pileupIdx-posMin) >= offsetRef){
                    refBuffer.push(f.references.at(rname).at(pileupIdx));
                }
                if((pileupIdx-posMin) >= offsetRead){
                    seqBuffer.push(seq[idx]);
                    qualBuffer.push(qual[idx]);
                }
                idx++; pileupIdx++;
            }
            break;
        case 'S':
        case 'I':
            idx += opLen;
            break;
        case 'D':
        case 'N':
            pileupIdx += opLen;
            break;
        case 'H':
        case 'P':
            break;  // these have been clipped
        default:
            throwErrorException("Bad CIGAR string");
        }
        opLen = 0;

    }

    //Process buffers
    while(refBuffer.size() && seqBuffer.size()){
        //Get front of queues
        char seq = seqBuffer.front();
        char qual = qualBuffer.front();
        char ref = refBuffer.front();
        seqBuffer.pop();
        qualBuffer.pop();
        refBuffer.pop();

        //Add qualities of mismatching bases, but ignore undefined reference bases
        if(seq != ref && ref != 'N'){
            score += qual;
            if(score > abortScore)
                return score;
        }

        //Add additional reference bases if available and read bases left
        if(!refBuffer.size() && seqBuffer.size() && f.references.at(rname).size() > pileupIdx){
            refBuffer.push(f.references.at(rname).at(pileupIdx));
            idx++; pileupIdx++;
        }
    }

    return score;
}

bool SAMRecord::isIndelEvidence(size_t maxIndelSize, size_t readOffset, const FASTAFile& f) const{

    if(tlen < readOffset + maxIndelSize || f.references.at(rname).length() < readOffset + maxIndelSize)
        return false;
    //Calc mismatching quality sum for original alignment
    size_t baseline = calcIndelScore(f, readOffset, readOffset, std::numeric_limits<std::size_t>::max());

    //Try insertion and deletion of all sizes smaller than maxIndelSize
    for(size_t i = 1; i <= maxIndelSize; ++i) {

        //Check if deletion aligns better than original
        if(calcIndelScore(f, readOffset, readOffset+i, baseline) <= baseline) {
            return true;
        }

        //Check if insertion aligns better than original
        if(calcIndelScore(f, readOffset+i, readOffset, baseline) <= baseline) {
            return true;
        }

    }

    //Original is best options
    return false;

}

void SAMRecord::addToPileupQueue(SAMPileupDeque *samPileupDeque_, const FASTAFile& f) const {
    if (samPileupDeque_->empty() == true) {
        throwErrorException("samPileupQueue is empty");
    }
    if ((samPileupDeque_->posMin() > posMin) || (samPileupDeque_->posMax() < posMax)) {
        throwErrorException("samPileupQueue does not overlap record");
    }

    size_t cigarIdx = 0;
    size_t cigarLen = cigar.length();
    size_t opLen = 0;  // length of current CIGAR operation
    size_t idx = 0;
    size_t pileupIdx = posMin - samPileupDeque_->posMin();

    const size_t maxIndelSize = 10; //TODO: Command line arguments

    bool justfoundEvidence=false; //Memorize if base before is evidence of indel. Used to ignore evidence if there is already an indel after.
    //TODO: just check in the cigar string and save some time by skippingisIndelEvidence

    for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++) {
        if (isdigit(cigar[cigarIdx])) {
            opLen = opLen*10 + (size_t)cigar[cigarIdx] - (size_t)'0';
            continue;
        }

        size_t front_hq_ctr=0; //Score before clip
        size_t back_hq_ctr=0;  //Score after clip

        //TODO: use offset parameter and add new cmd parameter for phred value
        const char HQ_SOFTCLIP_THRESHOLD = 29 + 64;


        switch (cigar[cigarIdx]) {
        case 'M':
        case '=':
        case 'X':
            for (size_t i = 0; i < opLen; i++) {
                samPileupDeque_->pileups_[pileupIdx].pos = samPileupDeque_->posMin() + pileupIdx;
                samPileupDeque_->pileups_[pileupIdx].seq += seq[idx];
                samPileupDeque_->pileups_[pileupIdx].qual += qual[idx];
                samPileupDeque_->pileups_[pileupIdx].ref = f.references.at(rname)[pileupIdx];

                //Check indel and update evidence counter
                if(isIndelEvidence(maxIndelSize, pileupIdx, f)){
                    samPileupDeque_->pileups_[pileupIdx].indelEvidence +=1;
                    justfoundEvidence = true;
                } else {
                    justfoundEvidence = false;
                }

                idx++; pileupIdx++;
            }
            break;
        case 'S':

            //Process softclips

            //std::cerr << "Softclips" << " detected!" << std::endl;

            //Clips to the right
            for(int l=0;l<opLen;++l){
                if(this->qual[idx+l] >= HQ_SOFTCLIP_THRESHOLD){
                    ++front_hq_ctr;
                } else {
                    break;
                }

            }

            //Clips to the left
            for(int l=opLen-1;l>=0;--l){
                if(this->qual[idx+l] >= HQ_SOFTCLIP_THRESHOLD){
                    ++back_hq_ctr;
                } else {
                    break;
                }

            }

            //Decide which to use
            if(idx-1 > 0)
                samPileupDeque_->pileups_[pileupIdx].hq_softcounter += front_hq_ctr;
            if(idx+opLen < this->qual.size())
                samPileupDeque_->pileups_[pileupIdx+1].hq_softcounter += back_hq_ctr;

        case 'I':

            //Take back evidence increment if insertion follows directly
            if(justfoundEvidence){
                samPileupDeque_->pileups_[pileupIdx-1].indelEvidence -= 1;
                justfoundEvidence = false;
            }
            idx += opLen;
            break;
        case 'D':
        case 'N':
            //Take back evidence increment if deletion follows directly
            if(justfoundEvidence){
                samPileupDeque_->pileups_[pileupIdx-1].indelEvidence -= 1;
                justfoundEvidence = false;
            }
            idx += opLen;
            pileupIdx += opLen;
            break;
        case 'H':
        case 'P':
            break;  // these have been clipped
        default:
            throwErrorException("Bad CIGAR string");
        }

        opLen = 0;
    }
}

bool SAMRecord::isMapped(void) const {
    return mapped_;
}

void SAMRecord::printLong(void) const {
    printShort();
    printf("isMapped: %d, ", mapped_);
    printf("posMin: %d, ", posMin);
    printf("posMax: %d\n", posMax);
}

void SAMRecord::printShort(void) const {
    printf("%s\t", qname.c_str());
    printf("%d\t", flag);
    printf("%s\t", rname.c_str());
    printf("%d\t", pos);
    printf("%d\t", mapq);
    printf("%s\t", cigar.c_str());
    printf("%s\t", rnext.c_str());
    printf("%d\t", pnext);
    printf("%" PRId64 "\t", tlen);
    printf("%s\t", seq.c_str());
    printf("%s\t", qual.c_str());
    printf("%s\t", opt.c_str());
    printf("\n");
}

void SAMRecord::printSeqWithPositionOffset(void) const {
    printf("%s %6d-%6d|", rname.c_str(), posMin, posMax);
    for (unsigned int i = 0; i < posMin; i++) { printf(" "); }
    printf("%s\n", seq.c_str());
}

void SAMRecord::check(void) {
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
    if (opt.empty() == true) { CALQ_LOG("opt is empty"); }

    // Check if this record is mapped
    if ((flag & 0x4) != 0) {
        mapped_ = false;
    } else {
        mapped_ = true;
        if (rname == "*" || pos == 0 || cigar == "*" || seq == "*" || qual == "*") {
            throwErrorException("Corrupted record");
        }
    }
}

}  // namespace calq


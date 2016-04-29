/** @file QualCodec.cc
 *  @brief This file contains the implementations of the QualEncoder and
 *         QualDecoder classes, respectively.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#include "QualCodec.h"
#include "Exceptions.h"

QualEncoder::QualEncoder(ofbitstream &ofbs)
    : ofbs(ofbs)
    , recordCnt(0)
{

}

QualEncoder::~QualEncoder(void)
{
    // Empty
}

static void extract(const std::string &seq,
                    const std::string &qual,
                    const std::string &cigar,
                    std::string &matchedSeq,
                    std::string &matchedQual,
                    std::string &insertedSeq,
                    std::string &insertedQual
                   )
{
    matchedSeq.clear();
    matchedQual.clear();
    insertedSeq.clear();
    insertedQual.clear();

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
            // add matching part to matchedSeq and matchedQual
            matchedSeq.append(seq, idx, opLen);
            matchedQual.append(qual, idx, opLen);
            idx += opLen;
            break;
        case 'I':
        case 'S':
            // add inserted bases and quality scores to insertedSeq and
            // insertedQual, respectively
            insertedSeq.append(seq, idx, opLen);
            insertedQual.append(qual, idx, opLen);
            idx += opLen; // skip inserted part
            break;
        case 'D':
        case 'N': {
            // inflate sequence and quality scores
            for (i = 0; i < opLen; i++) { 
                matchedSeq += "d";
                matchedQual += "d";
            }
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

void QualEncoder::encodeRecord(const SAMRecord &samRecord)
{
    recordCnt++;

    const uint32_t pos = samRecord.pos;
    const std::string cigar(samRecord.cigar);
    const std::string seq(samRecord.seq);
    const std::string qual(samRecord.qual);

    // check if this alignment is complete
    if (   (pos == 0)
        || (cigar.length() == 0 || cigar.compare("*") == 0)
        || (seq.length() == 0 || seq.compare("*") == 0)
        || (qual.length() == 0 || qual.compare("*") == 0)) {
        throwErrorException("Incomplete alignment");
    }

    // expand current sequence
    std::string matchedSeq("");
    std::string matchedQual("");
    std::string insertedSeq("");
    std::string insertedQual("");
    extract(seq, qual, cigar, matchedSeq, matchedQual, insertedSeq, insertedQual);
    //std::cout << "matchedSeq: " << matchedSeq << std::endl;
    //std::cout << "matchedQual: " << matchedQual << std::endl;
    //std::cout << "insertedSeq: " << insertedSeq << std::endl;
    //std::cout << "insertedQual: " << insertedQual << std::endl;

    // TODO: if this is the first record in a new block, simply allocate depths
    // vector; otherwise reallocate depths vector to cover the new region

    // TODO: accumulate depths

    // TODO: design/update nested quantizers for genomic columns

    // TODO: quantize the added quality score vector

    // TODO: apply and update quality score mask for current vector

    // TODO: perform Markov-chain prediction on current vector
    /*
    size_t alphabetSize = 40;
    size_t memorySize = 1;
    Predictor predictor(alphabetSize, memorySize);
    foreach (char q in qual) {
        ret = predictor.predict(std::string memory, std::string &predictedValue);
        if (ret == -1) std::cout << "cannot predict beginning of QS vector" << std::endl;
        predictor.update(std::string memory, std::string q);
        e = q - predictedValue;
    }
    */

    // TODO: adaptive clustering of QS vectors

    // TODO: pass QS vector to arithmetic or Huffman coder

    // ALTERNATIVE: try to predict quality score at a certain position and
    // then quantize the prediction error
}

size_t QualEncoder::finishBlock(void)
{
    size_t ret = 0;
    return ret;
}

QualDecoder::QualDecoder(ifbitstream &ifbs): ifbs(ifbs)
{

}

QualDecoder::~QualDecoder(void)
{

}

void QualDecoder::decodeBlock(std::vector<std::string> &qual)
{
    qual.push_back("*");
}


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
    , posMin(1)
    , posMax(1)
    , depths()
{

}

QualEncoder::~QualEncoder(void)
{
    // Empty
}

static std::string expand(const std::string &seq, const std::string &cigar)
{
    size_t cigarIdx = 0;
    size_t cigarLen = cigar.length();
    size_t opLen = 0; // length of current CIGAR operation
    size_t seqIdx = 0;
    std::string expandedSeq("");

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
            // add matching part to expanded sequence
            expandedSeq.append(seq, seqIdx, opLen);
            seqIdx += opLen;
            break;
        case 'I':
        case 'S':
            seqIdx += opLen; // skip inserted part
            break;
        case 'D':
        case 'N':
            // inflate expanded sequence
            for (i = 0; i < opLen; i++) { expandedSeq += "D"; }
            break;
        case 'H':
        case 'P':
            break; // these have been clipped
        default: 
            throwErrorException("Bad CIGAR string");
        }

        opLen = 0;
    }

    return expandedSeq;
}

void QualEncoder::encodeRecord(const SAMRecord &samRecord)
{
    const uint32_t pos = samRecord.pos;
    const std::string cigar(samRecord.cigar);
    const std::string seq(samRecord.seq);
    const std::string qual(samRecord.qual);

    recordCnt++;

    // check if this alignment is complete
    if (   (pos == 0)
        || (cigar.length() == 0 || cigar.compare("*") == 0)
        || (seq.length() == 0 || seq.compare("*") == 0)
        || (qual.length() == 0 || qual.compare("*") == 0)) {
        throwErrorException("Incomplete alignment");
    }

    // expand current sequence
    std::string exp = expand(seq, cigar);

    // TODO: if this is the first record in a new block, simply allocate depths
    // vector; otherwise reallocate depths vector to cover the new region

    // TODO: accumulate depths

    // TODO: design/update quantizers for genomic columns

    // TODO: quantize the added quality score vector

    // BEGIN meholli
    // TODO: perform Markov-chain prediction on current vector
    /*
    size_t alphabetSize = 40;
    size_t memorySize = 1;
    cq::Predictor predictor = new Predictor(alphabetSize, memorySize);
    foreach (char q in qual) {
        ret = predictor.predict(std::string memory, std::string &predictedValue);
        if (ret == -1) std::cout << "cannot predict beginning of QS vector" << std::endl;
        predictor.update(std::string memory, std::string q);
        e = q - predictedValue;
    }
    delete predictor;
    */
    // END meholli

    // TODO: pass QS vector to arithmetic or Huffman coder
}

size_t QualEncoder::finishBlock(void)
{
    size_t ret = 0;
    return ret;
}

QualDecoder::QualDecoder(ifbitstream &ifbs)
    : ifbs(ifbs)
{

}

QualDecoder::~QualDecoder(void)
{

}

void QualDecoder::decodeBlock(std::vector<std::string> &qual)
{
    qual.push_back("qual record 0");
}


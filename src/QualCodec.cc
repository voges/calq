/** @file QualCodec.cc
 *  @brief This file contains the implementations of the QualEncoder and
 *         QualDecoder classes.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#include "QualCodec.h"
#include "Exceptions.h"
#include "Predictor.h"

static const int ALPHABET_SIZE = 50;
static const int MEMORY_SIZE = 2;
static const int OFFSET = 33;
static const int QMAX = 83;
static const int QMIN = OFFSET;

QualEncoder::QualEncoder(ofbitstream &ofbs)
    : numEncodedRecords(0)
    , numInputRecords(0)
    , ofbs(ofbs)
    , predictor(ALPHABET_SIZE, MEMORY_SIZE, OFFSET)
{
    // empty
}

QualEncoder::~QualEncoder(void)
{
    // empty
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
    numInputRecords++;

    const uint32_t pos = samRecord.pos;
    const std::string cigar(samRecord.cigar);
    const std::string seq(samRecord.seq);
    const std::string qual(samRecord.qual);
    std::cout << "record " << numInputRecords << ": " << pos << " " << cigar << " " << seq << " " << qual << std::endl;

    // check if this alignment is complete
    if (   (pos == 0)
        || (cigar.length() == 0 || cigar.compare("*") == 0)
        || (seq.length() == 0 || seq.compare("*") == 0)
        || (qual.length() == 0 || qual.compare("*") == 0)) {
        std::cout << "Warning: Incomplete alignment; skipping record " << numInputRecords << std::endl;
        return;
    }

    // expand current sequence
    std::string matchedSeq("");
    std::string matchedQual("");
    std::string insertedSeq("");
    std::string insertedQual("");
    extract(seq, qual, cigar, matchedSeq, matchedQual, insertedSeq, insertedQual);

    // TODO: accumulate depths
    //quantizer1.updateDepths(pos, cigar);

    // TODO: quantize the added quality score vector
    // TODO: design/update nested quantizers for genomic columns
    /*for(std::string::size_type i = 0; i < qual.size(); ++i) {
        int q = (int)qual[i];
        int qTilde1 = quantizer1.quantize(q, pos+i);
        qualTilde1.push_back(qTilde1);
    }*/

    // TODO: apply and update quality score mask for current vector
    /*for(size_t i = 0; i < qualTilde1.size(); ++i) {
        int qTilde1 = (int)qual[i];
        int qTilde1Bar = mask.apply(qTilde1);
        qualTilde1Bar.push_back(qTilde1Bar);
    }*/

    // Markov-chain prediction for current current qual vector
    std::vector<int> err;
    std::vector<int> memory(MEMORY_SIZE, -1);

    for(std::string::size_type i = 0; i < qual.size(); ++i) {
        int q = (int)qual[i];

        if (i < MEMORY_SIZE) {
            std::fill(memory.begin(), memory.end(), -1); // this is not necessary
            err.push_back(q);
        } else {
            // fill memory
            for (size_t m = 0; m < MEMORY_SIZE; m++) {
                memory[m] = ((int)qual[i-1-m]);
            }

            // predict current q value, update predictor, and compute prediction
            // error
            int qHat = predictor.predict(memory);
            predictor.update(memory, q);
            int e = q - qHat;
            err.push_back(e);
        }
    }

    // TODO: adaptive clustering of QS vectors

    // TODO: run-length encoding?

    // TODO: pass QS vector to arithmetic or Huffman coder
    /*for(size_t i = 0; i < tmp.size(); ++i) {
        int q = (int)tmp[i];
        compress();
    }*/

    // ALTERNATIVE: try to predict quality score at a certain position and
    // then quantize the prediction error

    numEncodedRecords++;
}

size_t QualEncoder::finishBlock(void)
{
    size_t ret = 0;
    predictor.createCSVFile();
    predictor.reset();
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


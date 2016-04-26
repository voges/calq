#include "QualCodec.h"

QualEncoder::QualEncoder(ofbitstream &ofbs)
    : recordCnt(0)
    , posMin(1)
    , posMax(1)
    , ofbs(ofbs)
{

}

QualEncoder::~QualEncoder(void)
{

}

void QualEncoder::encodeRecord(const SAMRecord &samRecord)
{
    const uint32_t pos = samRecord.pos;
    std::string cigar(samRecord.cigar);
    std::string seq(samRecord.seq);
    std::string qual(samRecord.qual);

    std::cout << "line " << recordCnt << ": " << qual << std::endl;

    recordCnt++;

    // TODO: check if this alignment is complete
    if (   (pos == 0)
        || (cigar.length() == 0 || (cigar[0] == '*' && cigar[1] == '\0'))
        || (seq.length() == 0 || (seq[0] == '*' && seq[1] == '\0'))
        || (qual.length() == 0 || (qual[0] == '*' && qual[1] == '\0'))) {
        // alignment is incomplete
    }

    // TODO: expand current sequence
    //str_t *exp = str_new();
    //expand(exp, cigar, seq)) {

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
{

}

QualDecoder::~QualDecoder(void)
{

}

void QualDecoder::decodeBlock(std::vector<std::string> &qual)
{

}


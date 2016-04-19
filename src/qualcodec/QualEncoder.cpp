#include "qualcodec/QualEncoder.hpp"
#include <iostream>

cq::QualEncoder::QualEncoder(void)
{
    recordCnt = 0;
    posMin = 1;
    posMax = 1;
}

cq::QualEncoder::~QualEncoder(void)
{

}

void cq::QualEncoder::encodeRecord(const samrec_t * const samrec)
{
    const uint32_t pos = samrec->pos;
    std::string cigar(samrec->cigar);
    std::string seq(samrec->seq);
    std::string qual(samrec->qual);

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

    // TODO: perform Markov-chain prediction on current vector

    // TODO: pass QS vector to arithmetic or Huffman coder
}

size_t cq::QualEncoder::finishBlock(FILE *fp)
{
    size_t ret = 0;
    return ret;
}


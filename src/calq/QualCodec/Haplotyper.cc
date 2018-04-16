#include "Haplotyper.h"


//Returns score, takes seq and qual pileup and position
Haplotyper::Haplotyper(size_t sigma, size_t ploidy, size_t qualOffset, size_t nrQuantizers, size_t maxHQSoftclip_propagation, size_t minHQSoftclip_streak, size_t gaussRadius)
    : SIGMA(sigma), kernel(sigma), spreader(maxHQSoftclip_propagation, minHQSoftclip_streak), genotyper(ploidy, qualOffset, nrQuantizers){
    double THRESHOLD = 0.0000001;
    size_t size = this->kernel.calcMinSize(THRESHOLD, gaussRadius*2+1);
 
    buffer = FilterBuffer([this](size_t pos, size_t size) -> double{return this->kernel.calcValue(pos,size);}, size);
}

size_t Haplotyper::getOffset() const {
    return buffer.getOffset() + spreader.getOffset();
}

double Haplotyper::push(const std::string& seqPile, const std::string& qualPile, size_t hq_softclips, char reference){
    std::string refstr;
    refstr += reference;
    refstr += reference;
    double refProb = 1 - genotyper.getGenotypelikelihoods(seqPile, qualPile)[refstr];
    buffer.push(spreader.push(refProb, hq_softclips/qualPile.size()));
    return buffer.filter();
}

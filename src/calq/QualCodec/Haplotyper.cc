#include "Haplotyper.h"

#include <iostream>
#include <iomanip>
#include <sstream>

//Returns score, takes seq and qual pileup and position
Haplotyper::Haplotyper(size_t sigma, size_t ploidy, size_t qualOffset, size_t nrQuantizers, size_t maxHQSoftclip_propagation, size_t minHQSoftclip_streak, size_t gaussRadius)
    : SIGMA(sigma), kernel(sigma), spreader(maxHQSoftclip_propagation, minHQSoftclip_streak), genotyper(ploidy, qualOffset, nrQuantizers), nr_quantizers(nrQuantizers){
    double THRESHOLD = 0.0000001;
    size_t size = this->kernel.calcMinSize(THRESHOLD, gaussRadius*2+1);
 
    buffer = FilterBuffer([this](size_t pos, size_t size) -> double{return this->kernel.calcValue(pos,size);}, size);
}

size_t Haplotyper::getOffset() const {
    return buffer.getOffset() + spreader.getOffset();
}

size_t Haplotyper::push(const std::string& seqPile, const std::string& qualPile, size_t hq_softclips, char reference){
    static CircularBuffer<std::string> debug(this->getOffset(), "\n");

    if(seqPile.empty()) {
        buffer.push(spreader.push(0.0,0));
        return 0;
    }
    std::string refstr;
    refstr += reference;
    refstr += reference;
    double refProb;
    if(reference == 'N')
        refProb = 1.0;
    else
        refProb= 1 - genotyper.getGenotypelikelihoods(seqPile, qualPile).at(refstr);
    buffer.push(spreader.push(refProb, hq_softclips/qualPile.size()));
    double activity = buffer.filter();

    std::stringstream s;

    s << refstr << " " << seqPile << " " << qualPile << " " << hq_softclips << " ";

    s << std::fixed << std::setw( 6 ) << std::setprecision( 4 )
              << std::setfill( '0' ) << refProb;

    std::cerr << debug.push(s.str()) << " " << std::fixed << std::setw( 6 ) << std::setprecision( 4 )
              << std::setfill( '0' ) << activity << std::endl;


    return activity * nr_quantizers;
}


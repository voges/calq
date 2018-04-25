#include "Haplotyper.h"

#include <iostream>
#include <iomanip>
#include <sstream>

#include <cmath>

//Returns score, takes seq and qual pileup and position
Haplotyper::Haplotyper(size_t sigma, size_t ploidy, size_t qualOffset, size_t nrQuantizers, size_t maxHQSoftclip_propagation, size_t minHQSoftclip_streak, size_t gaussRadius)
    : SIGMA(sigma), kernel(sigma), spreader(maxHQSoftclip_propagation, minHQSoftclip_streak), genotyper(ploidy, qualOffset, nrQuantizers), nr_quantizers(nrQuantizers), polyploidy(ploidy){
    double THRESHOLD = 0.0000001;
    size_t size = this->kernel.calcMinSize(THRESHOLD, gaussRadius*2+1);
 
    buffer = FilterBuffer([this](size_t pos, size_t size) -> double{return this->kernel.calcValue(pos,size);}, size);

    NO_INDEL_LIKELIHOOD = log10(1.0-pow(10,-4.5));
    INDEL_LIKELIHOOD = pow(10,-4.5);
}

size_t Haplotyper::getOffset() const {
    return buffer.getOffset() + spreader.getOffset();
}

std::map<size_t, double> Haplotyper::calcIndelLikelihoods(size_t numberOfEvidence){
    double denominator = -log10(polyploidy);
    std::map<size_t, double> likelihoods;


    if(numberOfEvidence == 0) {
        for(int variantCount = 0; variantCount <= polyploidy; ++variantCount) {
            likelihoods[variantCount] = 1.0;
        }
        return likelihoods;
    }

    likelihoods[0] = numberOfEvidence * NO_INDEL_LIKELIHOOD;
    for(int variantCount = 1; variantCount <= polyploidy; ++variantCount) {
        double refLikelihood = NO_INDEL_LIKELIHOOD + log10(polyploidy - numberOfEvidence);
        double altLikelihood = INDEL_LIKELIHOOD + log10(numberOfEvidence);
        likelihoods[variantCount] = numberOfEvidence * (log10(exp10(refLikelihood)+exp10(altLikelihood)) + denominator);
        likelihoods[variantCount] = pow(10,(likelihoods[variantCount]/-10.0));
    }
    return likelihoods;
}

size_t Haplotyper::push(const std::string& seqPile, const std::string& qualPile, size_t hq_softclips, size_t indels, char reference){
    static CircularBuffer<std::string> debug(this->getOffset(), "\n");

    if(seqPile.empty()) {
        buffer.push(spreader.push(0.0,0));
        return 0;
    }
    std::string refstr;
    refstr += reference;
    refstr += reference;
    double refProb;
    if(reference == 'N'){
        auto likelihoods = genotyper.getGenotypelikelihoods(seqPile, qualPile);
        refProb = 1-(likelihoods.at("AA") + likelihoods.at("CC") + likelihoods.at("GG") + likelihoods.at("TT"))/4.0; //Probability of hom
    }
    else {
        double snpProb = 1 - genotyper.getGenotypelikelihoods(seqPile, qualPile).at(refstr);
        double indelProb = 1 - calcIndelLikelihoods(indels).at(0);

        refProb = std::max(snpProb, indelProb);
    }
    buffer.push(spreader.push(refProb, hq_softclips/qualPile.size()));
    double activity = buffer.filter();

  /*  std::stringstream s;

    s << refstr << " " << seqPile << " " << qualPile << " " << hq_softclips << " ";

    s << std::fixed << std::setw( 6 ) << std::setprecision( 4 )
              << std::setfill( '0' ) << refProb;

 //   std::cerr << debug.push(s.str()) << " " << std::fixed << std::setw( 6 ) << std::setprecision( 4 )
  //            << std::setfill( '0' ) << activity << std::endl;*/

    if(hq_softclips > 0)
        std::cerr << hq_softclips << " detected!" << std::endl;


    return activity * nr_quantizers;
}


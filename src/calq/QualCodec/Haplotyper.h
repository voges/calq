#include <functional>
#include <cmath>

#include "FilterBuffer.h"
#include "Genotyper.h"
#include "BaseSpreader.h"

class Haplotyper {
private:
    //Saving kernel and old raw quality scores
    FilterBuffer buffer;

    GaussKernel kernel;

    BaseSpreader spreader;

    calq::Genotyper genotyper;

    const size_t SIGMA;

    size_t nr_quantizers;

public:

    //Init
    Haplotyper(size_t sigma, size_t ploidy, size_t qualOffset, size_t nrQuantizers, size_t maxHQSoftclip_propagation, size_t minHQSoftclip_streak, size_t gaussRadius);

    //Returns offset between activity scores' position and front
    size_t getOffset() const;

    //Pushes new activity score calculated using parameters and returns filtered acticityscore for (pos-offset)
    size_t push(const std::string& seqPile, const std::string& qualPile, size_t hq_softclips, char reference);

};

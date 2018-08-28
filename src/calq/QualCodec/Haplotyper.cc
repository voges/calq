/** @file Haplotyper.cc
 *  @brief This file contains the implementation of the Haplotyper class.
 */

// Copyright 2015-2018 Leibniz Universitaet Hannover

#include <iostream>
#include <iomanip>
#include <sstream>
#include <limits>
#include <algorithm>

#include "QualCodec/Haplotyper.h"

#include "Common/ErrorExceptionReporter.h"

// ----------------------------------------------------------------------------------------------------------------------

namespace calq {

// Returns score, takes seq and qual pileup and position
Haplotyper::Haplotyper(size_t sigma, size_t ploidy, size_t qualOffset, size_t nrQuantizers, size_t maxHQSoftclip_propagation,
                       size_t minHQSoftclip_streak, size_t gaussRadius,
                       bool debug, bool squashed, Options::FilterType filterType) : SIGMA(sigma),
                                                                                    spreader(maxHQSoftclip_propagation, minHQSoftclip_streak, squashed),
                                                                                    genotyper(static_cast<const int &>(ploidy),
                                                                                              static_cast<const int &>(qualOffset),
                                                                                              static_cast<const int &>(nrQuantizers)),
                                                                                    nr_quantizers(nrQuantizers), polyploidy(ploidy), DEBUG(debug),
                                                                                    squashedActivity(squashed) {
    if (filterType == Options::FilterType::GAUSS) {
        GaussKernel kernel(sigma);
        double THRESHOLD = 0.0000001;
        size_t size = kernel.calcMinSize(THRESHOLD, gaussRadius * 2 + 1);

        buffer = FilterBuffer([kernel](size_t pos, size_t size) -> double { return kernel.calcValue(pos, size); }, size);
        localDistortion = kernel.calcValue((size - 1) / 2, size);
    } else if (filterType == Options::FilterType::RECTANGLE) {
        RectangleKernel kernel(sigma);
        size_t size = kernel.calcMinSize(gaussRadius * 2 + 1);

        buffer = FilterBuffer([kernel](size_t pos, size_t size) -> double { return kernel.calcValue(pos, size); }, size);
        localDistortion = kernel.calcValue((size - 1) / 2, size);
    } else {
        throwErrorException("FilterType not supported by haplotyper");
    }
}

// ----------------------------------------------------------------------------------------------------------------------

size_t Haplotyper::getOffset() const {
    return buffer.getOffset() + spreader.getOffset() - 1;
}

// ----------------------------------------------------------------------------------------------------------------------


double log10sum(double a, double b) {
    if (a > b) {
        log10sum(b, a);
    } else if (a == -std::numeric_limits<double>::infinity()) {
        return b;
    }
    return b + log10(1 + pow(10.0, -(b - a)));
}

// ----------------------------------------------------------------------------------------------------------------------

// Calc priors like in GATK
std::vector<double> Haplotyper::calcPriors(double hetero) {
    hetero = log10(hetero);
    std::vector<double> result(polyploidy + 1, hetero);
    double sum = -std::numeric_limits<double>::infinity();
    for (size_t i = 1; i < polyploidy + 1; ++i) {
        result[i] -= log10(i);
        sum = log10sum(sum, result[i]);
    }
    result[0] = log10(1.0 - exp10(sum));

    return result;
}

// ----------------------------------------------------------------------------------------------------------------------

std::vector<double> Haplotyper::calcNonRefLikelihoods(char ref, const std::string &seqPile, const std::string &qualPile) {
    std::vector<double> result(polyploidy + 1, 0.0);
    std::map<std::string, double> SNPlikelihoods = genotyper.getGenotypelikelihoods(seqPile, qualPile);

    for (const auto &m : SNPlikelihoods) {
        size_t altCount = polyploidy - std::count(m.first.begin(), m.first.end(), ref);
        result[altCount] += m.second;
    }

    for (double &i : result) {
        i = log10(i);
    }

    return result;
}

// ----------------------------------------------------------------------------------------------------------------------

double Haplotyper::calcActivityScore(char ref, const std::string &seqPile, const std::string &qualPile, double heterozygosity) {
    if (ref == 'N') {
        return 1.0;
    }
    std::vector<double> likelihoods = calcNonRefLikelihoods(ref, seqPile, qualPile);
    static std::vector<double> priors;
    if (priors.empty()) {
        priors = calcPriors(heterozygosity);
    }

    // --------------Calc Posteriors like in GATK ------------------
    double posteriori0 = likelihoods[0] + priors[0];
    bool map0 = true;
    for (size_t i = 1; i < polyploidy + 1; ++i) {
        if (likelihoods[i] + priors[i] > posteriori0) {
            map0 = false;
            break;
        }
    }

    if (map0) {
        return 0.0;
    }

    double altLikelihoodSum = -std::numeric_limits<double>::infinity();
    double altPriorSum = -std::numeric_limits<double>::infinity();

    for (size_t i = 1; i < polyploidy + 1; ++i) {
        altLikelihoodSum = log10sum(altLikelihoodSum, likelihoods[i]);
        altPriorSum = log10sum(altPriorSum, priors[i]);
    }

    double altPosteriorSum = altLikelihoodSum + altPriorSum;
    // Normalize
    posteriori0 = posteriori0 - log10sum(altPosteriorSum, posteriori0);

    // --------------------------------------------------------

    /*  const double CUTOFF = -1.0;

      if (posteriori0 > CUTOFF)
          return 0.0;*/

    return 1.0 - exp10(posteriori0);
}

// ----------------------------------------------------------------------------------------------------------------------

size_t Haplotyper::push(const std::string &seqPile, const std::string &qualPile, size_t hq_softclips, char reference) {
    // Empty input
    if (seqPile.empty()) {
        buffer.push(spreader.push(0.0, 0));
        return 0;
    }

    // Build reference string

    const double HETEROZYGOSITY = 1.0 / 1000.0;

    double altProb = calcActivityScore(reference, seqPile, qualPile, HETEROZYGOSITY);


    // Filter activity score
    buffer.push(spreader.push(altProb, hq_softclips / qualPile.size()));
    double activity = buffer.filter();
    if (squashedActivity) {
        activity = std::min(activity, 1.0);
    }

    size_t quant = getQuantizerIndex(activity);

    if (DEBUG) {
        static CircularBuffer<std::string> debug(this->getOffset(), "\n");
        std::stringstream s;

        s << reference << " " << seqPile << " ";

        s << std::fixed << std::setw(6) << std::setprecision(4)
          << std::setfill('0') << altProb;

        std::string out = debug.push(s.str());
        if (out != "\n") {
            std::cerr << out << " " << std::fixed << std::setw(6) << std::setprecision(4)
                      << std::setfill('0') << activity << " " << quant << std::endl;
        }

        if (hq_softclips > 0)
            std::cerr << hq_softclips << " detected!" << std::endl;
    }

    return quant;
}

// ----------------------------------------------------------------------------------------------------------------------

size_t Haplotyper::getQuantizerIndex(double activity) {
    return (size_t) std::min(std::floor((activity / localDistortion) * nr_quantizers), static_cast<double>(nr_quantizers - 1));
}
}  // namespace calq

// ----------------------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------------------


/**
 * @file haplotyper.cc
 */

#include "haplotyper.h"
#include <algorithm>
#include <iostream>
#include <limits>
#include "data-structures.h"
#include "exceptions.h"

namespace calq {

Haplotyper::Haplotyper(const size_t sigma, const size_t ploidy, const size_t qualOffset, const size_t numQuantizers,
                       const size_t maxHqSoftclipPropagation, const size_t hqSoftclipStreak, size_t const filterCutOff,
                       const bool squashed, const FilterType filterType)
    : spreader_(maxHqSoftclipPropagation, hqSoftclipStreak, squashed),
      genotyper_(static_cast<const int&>(ploidy), static_cast<const int&>(qualOffset),
                 static_cast<const int&>(numQuantizers)),
      numQuantizers_(numQuantizers),
      polyploidy_(ploidy),
      squashedActivity_(squashed) {
    if (filterType == FilterType::GAUSS) {
        GaussKernel kernel(sigma);
        double THRESHOLD = 0.0000001;
        size_t minSize = kernel.calcMinSize(THRESHOLD, filterCutOff * 2 + 1);

        buffer_ =
            FilterBuffer([kernel](size_t pos, size_t size) -> double { return kernel.calcValue(pos, size); }, minSize);
        localDistortion_ = kernel.calcValue((minSize - 1) / 2, minSize);
    } else if (filterType == FilterType::RECT) {
        RectangleKernel kernel(sigma);
        size_t minSize = kernel.calcMinSize(filterCutOff * 2 + 1);

        buffer_ =
            FilterBuffer([kernel](size_t pos, size_t size) -> double { return kernel.calcValue(pos, size); }, minSize);
        localDistortion_ = kernel.calcValue((minSize - 1) / 2, minSize);
    } else {
        throwErrorException("FilterType not supported by haplotyper");
    }
}

size_t Haplotyper::getOffset() const { return buffer_.getOffset() + spreader_.getOffset() - 1; }

double log10sum(const double a, const double b) {
    if (a > b) {
        return log10sum(b, a);
    } else if (a == -std::numeric_limits<double>::infinity()) {
        return b;
    }
    return b + log10(1 + pow(10.0, -(b - a)));
}

std::vector<double> Haplotyper::calcPriors(double hetero) {
    // Calculate priors like in GATK
    hetero = log10(hetero);
    std::vector<double> result(polyploidy_ + 1, hetero);
    double sum = -std::numeric_limits<double>::infinity();
    for (size_t i = 1; i < polyploidy_ + 1; ++i) {
        result[i] -= log10(i);
        sum = log10sum(sum, result[i]);
    }
    result[0] = log10(1.0 - pow(10, sum));

    return result;
}

std::vector<double> Haplotyper::calcNonRefLikelihoods(const char ref, const std::string& seqPile,
                                                      const std::string& qualPile) {
    std::vector<double> result(polyploidy_ + 1, 0.0);
    std::map<std::string, double> snpLikelihoods = genotyper_.getGenotypeLikelihoods(seqPile, qualPile);

    for (const auto& m : snpLikelihoods) {
        size_t altCount = polyploidy_ - std::count(m.first.begin(), m.first.end(), ref);
        result[altCount] += m.second;
    }

    for (double& i : result) {
        i = log10(i);
    }

    return result;
}

double Haplotyper::calcActivityScore(const char ref, const std::string& seqPile, const std::string& qualPile,
                                     const double heterozygosity) {
    if (ref == 'N') {
        return 1.0;
    }

    std::vector<double> likelihoods = calcNonRefLikelihoods(ref, seqPile, qualPile);
    static std::vector<double> priors;
    if (priors.empty()) {
        priors = calcPriors(heterozygosity);
    }

    // Calculate posteriors like in GATK
    double posterior = likelihoods[0] + priors[0];
    bool map0 = true;
    for (size_t i = 1; i < polyploidy_ + 1; ++i) {
        if (likelihoods[i] + priors[i] > posterior) {
            map0 = false;
            break;
        }
    }

    if (map0) {
        return 0.0;
    }

    double altLikelihoodSum = -std::numeric_limits<double>::infinity();
    double altPriorSum = -std::numeric_limits<double>::infinity();

    for (size_t i = 1; i < polyploidy_ + 1; ++i) {
        altLikelihoodSum = log10sum(altLikelihoodSum, likelihoods[i]);
        altPriorSum = log10sum(altPriorSum, priors[i]);
    }

    double altPosteriorSum = altLikelihoodSum + altPriorSum;

    // Normalize
    posterior = posterior - log10sum(altPosteriorSum, posterior);

    return 1.0 - pow(10, posterior);
}

size_t Haplotyper::push(const std::string& seqPile, const std::string& qualPile, const size_t hqSoftclips,
                        const char reference) {
    // Empty input
    if (seqPile.empty()) {
        buffer_.push(spreader_.push(0.0, 0));
        return 0;
    }

    // Build reference string

    const double HETEROZYGOSITY = 1.0 / 1000.0;

    double altProb = calcActivityScore(reference, seqPile, qualPile, HETEROZYGOSITY);

    // Filter activity score
    buffer_.push(spreader_.push(altProb, hqSoftclips / qualPile.size()));
    double activity = buffer_.filter();
    if (squashedActivity_) {
        activity = std::min(activity, 1.0);
    }

    size_t quant = getQuantizerIndex(activity);

    return quant;
}

size_t Haplotyper::getQuantizerIndex(const double activity) {
    return (size_t)std::min(std::floor((activity / localDistortion_) * numQuantizers_),
                            static_cast<double>(numQuantizers_ - 1));
}

}  // namespace calq

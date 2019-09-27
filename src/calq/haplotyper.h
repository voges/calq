/**
 * @file haplotyper.h
 */

#ifndef CALQ_HAPLOTYPER_H_
#define CALQ_HAPLOTYPER_H_

#include <cmath>
#include <functional>
#include <map>
#include <string>
#include <vector>
#include "filter-buffer.h"
#include "genotyper.h"
#include "softclip-spreader.h"

namespace calq {

enum struct FilterType;

class Haplotyper {
   public:
    Haplotyper() = delete;
    Haplotyper(size_t sigma, size_t ploidy, size_t qualOffset, size_t numQuantizers, size_t maxHqSoftclipPropagation,
               size_t hqSoftclipStreak, size_t filterCutOff, bool squashed, FilterType filterType);
    Haplotyper(const Haplotyper&) = delete;
    Haplotyper& operator=(const Haplotyper&) = delete;
    Haplotyper(Haplotyper&&) = delete;
    Haplotyper& operator=(Haplotyper&&) = delete;
    ~Haplotyper() = default;

    /**
     * Returns offset between activity scores position and front
     */
    size_t getOffset() const;

    /**
     * Pushes new activity score calculated using parameters and returns filtered activity score for (pos - offset)
     */
    size_t push(const std::string& seqPile, const std::string& qualPile, size_t hqSoftclips, char reference);

    std::vector<double> calcNonRefLikelihoods(char ref, const std::string& seqPile, const std::string& qualPile);

    double calcActivityScore(char ref, const std::string& seqPile, const std::string& qualPile, double heterozygosity);

    std::vector<double> calcPriors(double hetero);

   private:
    FilterBuffer buffer_;  /// Saving kernel and old raw quality scores
    SoftclipSpreader spreader_;
    Genotyper genotyper_;
    size_t numQuantizers_;
    size_t polyploidy_;
    const bool squashedActivity_;
    double localDistortion_;

    size_t getQuantizerIndex(double activity);
};

}  // namespace calq

#endif  // CALQ_HAPLOTYPER_H_

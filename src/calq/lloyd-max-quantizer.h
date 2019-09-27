/**
 * @file lloyd-max-quantizer.h
 */

#ifndef CALQ_LLOYD_MAX_QUANTIZER_H_
#define CALQ_LLOYD_MAX_QUANTIZER_H_

#include <cmath>
#include <vector>
#include "exceptions.h"
#include "probability-distribution.h"
#include "quantizer.h"

namespace calq {

class LloydMaxQuantizer : public Quantizer {
   public:
    LloydMaxQuantizer() = delete;
    explicit LloydMaxQuantizer(size_t steps);
    LloydMaxQuantizer(const LloydMaxQuantizer &) = delete;
    LloydMaxQuantizer &operator=(const LloydMaxQuantizer &) = delete;
    LloydMaxQuantizer(LloydMaxQuantizer &&) = delete;
    LloydMaxQuantizer &operator=(LloydMaxQuantizer &&) = delete;
    ~LloydMaxQuantizer() override = default;

    /**
     * Create quantization boundaries from a probability distribution and fill LUT
     */
    void build(const ProbabilityDistribution &pdf);

   private:
    std::vector<double> borders_;  /// Decision thresholds
    std::vector<double> values_;   /// Representative values for each interval
    size_t steps_;

    /**
     * Fill LUT of base class
     */
    void fillLUT(const ProbabilityDistribution &pdf);

    /**
     * Calculate the centroid in a region of a PDF
     */
    static double centroid(size_t left, size_t right, const ProbabilityDistribution &pdf);

    /**
     * Calculate quantization borders using PDF
     */
    void calcBorders(const ProbabilityDistribution &pdf);
};

}  // namespace calq

#endif  // CALQ_LLOYD_MAX_QUANTIZER_H_

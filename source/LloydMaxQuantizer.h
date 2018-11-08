/** @file LloydMaxQuantizer.h
 *  @brief This file contains the definition of a Lloyd quantizer
 */

// Copyright 2015-2018 Leibniz Universitaet Hannover

#ifndef CALQ_QUALCODEC_QUANTIZERS_LLOYDMAXQUANTIZER_H_
#define CALQ_QUALCODEC_QUANTIZERS_LLOYDMAXQUANTIZER_H_

#include <cmath>

#include <vector>

#include "Quantizer.h"
#include "ProbabilityDistribution.h"
#include "Exceptions.h"

// ----------------------------------------------------------------------------------------------------------------------

namespace calq {

class LloydMaxQuantizer : public Quantizer {
 private:
    // Decision thresholds
    std::vector<double> borders;

    // Representative values for each interval
    std::vector<double> values;

    size_t steps;

    // Fill LUT of base class
    void fillLUT(const ProbabilityDistribution &pdf);

    // Calculates the centroid in a region of a pdf
    double calcCentroid(size_t left, size_t right, const ProbabilityDistribution &pdf);

    // Calculates quantization borders using pdf
    void calcBorders(const ProbabilityDistribution &pdf);
 public:
    explicit LloydMaxQuantizer(size_t steps);

    // Creates quantization boundaries from a probability distribution and fills LUT
    void build(const ProbabilityDistribution &pdf);
};
}  // namespace calq

// ----------------------------------------------------------------------------------------------------------------------

#endif  // CALQ_QUALCODEC_QUANTIZERS_LLOYDMAXQUANTIZER_H_

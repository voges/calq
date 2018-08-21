/** @file ProbabilityDistribution.h
 *  @brief This file contains the definition of a probability distribution structure
 */

// Copyright 2015-2018 Leibniz Universitaet Hannover

#ifndef CALQ_QUALCODEC_QUANTIZERS_PROBABILITYDISTRIBUTION_H_
#define CALQ_QUALCODEC_QUANTIZERS_PROBABILITYDISTRIBUTION_H_

#include <stddef.h>

#include <vector>


// ----------------------------------------------------------------------------------------------------------------------

namespace calq {

// Samples a probability distribution
class ProbabilityDistribution {
 private:
    std::vector<size_t> pdf;

    // Offset
    size_t rangeMin;

 public:
    ProbabilityDistribution(size_t rangeMin, size_t rangeMax);

    // Increases counter of value
    void addToPdf(size_t qualScore, size_t number = 1);

    // Reset all counters to zero
    void resetPdf();

    // Returns size of interval of valid values
    size_t size() const;

    // Returns current count of value
    size_t getCount(size_t value) const;

    // Returns current count of index (starting at zero)
    size_t operator[](size_t index) const;

    // Lowest possible value
    size_t getRangeMin() const;

    // Greatest possible value
    size_t getRangeMax() const;
};
}  // namespace calq

// ----------------------------------------------------------------------------------------------------------------------

#endif  // CALQ_QUALCODEC_QUANTIZERS_PROBABILITYDISTRIBUTION_H_

/** @file ProbabilityDistribution.h
 *  @brief This file contains the definition of a probability distribution structure
 */

// Copyright 2015-2018 Leibniz Universitaet Hannover

#ifndef CALQ_QUALCODEC_QUANTIZERS_PROBABILITYDISTRIBUTION_H_
#define CALQ_QUALCODEC_QUANTIZERS_PROBABILITYDISTRIBUTION_H_

#include <vector>


// ----------------------------------------------------------------------------------------------------------------------

namespace calq {

class ProbabilityDistribution {
 private:
    std::vector<size_t> pdf;
    size_t rangeMin;

 public:
    ProbabilityDistribution(size_t rangeMin, size_t rangeMax);

    void addToPdf(size_t qualScore, size_t number = 1);

    void resetPdf();

    size_t getCount(size_t value);
};
}  // namespace calq

// ----------------------------------------------------------------------------------------------------------------------

#endif  // CALQ_QUALCODEC_QUANTIZERS_PROBABILITYDISTRIBUTION_H_

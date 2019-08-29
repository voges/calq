/**
 * @file probability-distribution.h
 */

#ifndef CALQ_PROBABILITY_DISTRIBUTION_H_
#define CALQ_PROBABILITY_DISTRIBUTION_H_

#include <cstddef>
#include <vector>

namespace calq {

class ProbabilityDistribution {
   public:
    ProbabilityDistribution(size_t rangeMin, size_t rangeMax);
    void add(size_t value, size_t count = 1);
    void reset();
    size_t size() const;
    size_t count(size_t value) const;
    size_t operator[](size_t index) const;
    size_t rangeMin() const;
    size_t rangeMax() const;

   private:
    std::vector<size_t> pdf_;
    size_t rangeMin_;
};

}  // namespace calq

#endif  // CALQ_PROBABILITY_DISTRIBUTION_H_

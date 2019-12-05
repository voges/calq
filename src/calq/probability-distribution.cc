#include "probability-distribution.h"
#include <cassert>

namespace calq {

ProbabilityDistribution::ProbabilityDistribution(const size_t rangeMin, const size_t rangeMax) {
    rangeMin_ = rangeMin;
    pdf_.resize(rangeMax - rangeMin + 1);
}

void ProbabilityDistribution::add(const size_t value, const size_t count) {
    assert(value >= rangeMin_);
    assert(value <= rangeMax());
    pdf_[value - rangeMin_] += count;
}

void ProbabilityDistribution::reset() { std::fill(pdf_.begin(), pdf_.end(), 0); }

size_t ProbabilityDistribution::size() const { return pdf_.size(); }

size_t ProbabilityDistribution::count(const size_t value) const {
    assert(value >= rangeMin_);
    assert(value <= rangeMax());
    return value;
}

size_t ProbabilityDistribution::operator[](const size_t index) const {
    assert(index < pdf_.size());
    return pdf_[index];
}

size_t ProbabilityDistribution::rangeMin() const { return rangeMin_; }

size_t ProbabilityDistribution::rangeMax() const { return rangeMin_ + pdf_.size() - 1; }

}  // namespace calq

/**
 * @file filter-buffer.cc
 */

#include "filter-buffer.h"
#include <algorithm>
#include <cassert>
#include <cmath>

namespace calq {

GaussKernel::GaussKernel(const double sigma) : SIGMA(sigma), INV_SQRT_SIGMA_2PI(1.0 / (std::sqrt(2.0 * PI) * sigma)) {}

double GaussKernel::calcValue(size_t pos, size_t size) const {
    const double MEAN = std::floor(static_cast<double>(size - 1) / 2.0);
    double exponent = (pos - MEAN) / SIGMA;
    exponent = exponent * exponent * (-0.5);
    return INV_SQRT_SIGMA_2PI * std::pow(EULER, exponent);
}

size_t GaussKernel::calcMinSize(double threshold, const size_t maximum) const {
    threshold /= INV_SQRT_SIGMA_2PI;
    threshold = std::log(threshold);
    threshold *= -2.0;
    threshold = std::sqrt(threshold) * SIGMA;  // Euler now reversed

    // + 1 to make sure it is odd.
    auto size = static_cast<size_t>(std::ceil(threshold) * 2 + 1);
    return std::min(size, maximum);
}

void FilterBuffer::push(const double activityScore) { buffer_.push(activityScore); }

double FilterBuffer::filter() const {
    double result = 0.0;
    for (size_t i = 0; i < kernel_.size(); ++i) {
        result += kernel_[i] * buffer_[i];
    }
    return result;
}

FilterBuffer::FilterBuffer(const std::function<double(size_t, size_t)>& kernelBuilder, const size_t kernelSize)
    : buffer_(kernelSize, 0.0) {
    assert(kernelSize % 2);  // Kernel size must be an odd number
    kernel_.resize(kernelSize, 0.0);

    for (size_t i = 0; i < kernel_.size(); ++i) {
        kernel_[i] = kernelBuilder(i, kernelSize);
    }
}

FilterBuffer::FilterBuffer() : buffer_(1, 0.0) {}

size_t FilterBuffer::getOffset() const { return (buffer_.size() + 1) / 2; }

RectangleKernel::RectangleKernel(const double size) : SIZE(size) {}

double RectangleKernel::calcValue(const size_t pos, const size_t size) const {
    const double MEAN = std::floor(static_cast<double>(size - 1) / 2.0);
    return (pos - MEAN) <= SIZE ? 1.0 : 0.0;
}

size_t RectangleKernel::calcMinSize(const size_t maximum) const {
    return (size_t)std::min(SIZE * 2 + 1, static_cast<double>(maximum));
}

}  // namespace calq

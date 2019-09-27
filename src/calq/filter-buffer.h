/**
 * @file filter-buffer.h
 */

#ifndef CALQ_FILTER_BUFFER_H_
#define CALQ_FILTER_BUFFER_H_

#include <cstddef>
#include <functional>
#include <vector>
#include "circular-buffer.h"

namespace calq {

/**
 * Samples a normal distribution
 */
class GaussKernel {
   public:
    GaussKernel() = delete;
    explicit GaussKernel(double sigma = 1.0);
    GaussKernel(const GaussKernel &) = default;
    GaussKernel &operator=(const GaussKernel &) = delete;
    GaussKernel(GaussKernel &&) = delete;
    GaussKernel &operator=(GaussKernel &&) = delete;
    ~GaussKernel() = default;

    /**
     * Get Gauss value at position 'pos' and buffer size 'size' with mean = size / 2
     */
    double calcValue(size_t pos, size_t size) const;

    /**
     * Calculates how big a buffer must be to contain all values above threshold. No size greater than maximum is
     * returned.
     */
    size_t calcMinSize(double threshold, size_t maximum = 101) const;

   private:
    const double SIGMA;
    constexpr static double PI = 3.14159265359;
    constexpr static double EULER = 2.71828182846;
    const double INV_SQRT_SIGMA_2PI;
};

/**
 * Samples a uniform distribution
 */
class RectangleKernel {
   public:
    RectangleKernel() = delete;
    explicit RectangleKernel(double size = 1.0);
    RectangleKernel(const RectangleKernel &) = default;
    RectangleKernel &operator=(const RectangleKernel &) = delete;
    RectangleKernel(RectangleKernel &&) = delete;
    RectangleKernel &operator=(RectangleKernel &&) = delete;
    ~RectangleKernel() = default;

    /**
     * Get rect value at position 'pos' and buffer size 'size' with mean = size / 2
     */
    double calcValue(size_t pos, size_t size) const;

    /**
     * Calculates how big a buffer must be to contain all values above threshold. No size greater than maximum is
     * returned.
     */
    size_t calcMinSize(size_t maximum = 101) const;

   private:
    const double SIZE;
};

/**
 * Filters an input signal using a filter kernel
 */
class FilterBuffer {
   public:
    FilterBuffer();
    FilterBuffer(const std::function<double(size_t, size_t)> &kernelBuilder, size_t kernelSize);
    FilterBuffer(const FilterBuffer &) = delete;
    FilterBuffer &operator=(const FilterBuffer &) = default;
    FilterBuffer(FilterBuffer &&) = delete;
    // TODO(Fabian): The move assignment operator should be marked as deleted. However, it is still used in the code.
    FilterBuffer &operator=(FilterBuffer &&) = default;
    ~FilterBuffer() = default;

    /**
     * New activity score in pipeline
     */
    void push(double activityScore);

    /**
     * Calculate filter score at offset position
     */
    double filter() const;

    /**
     * Distance between buffer center and borders
     */
    size_t getOffset() const;

   private:
    std::vector<double> kernel_;
    CircularBuffer<double> buffer_;
};

}  // namespace calq

#endif  // CALQ_FILTER_BUFFER_H_

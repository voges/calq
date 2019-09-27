/**
 * @file uniform-min-max-quantizer.cc
 */

#include "uniform-min-max-quantizer.h"

namespace calq {

UniformMinMaxQuantizer::UniformMinMaxQuantizer(const int minValue, const int maxValue, const int numSteps)
    : UniformQuantizer(minValue, maxValue, numSteps) {
    // Change the smallest and largest reconstruction values

    int smallestIndex = 0;
    int largestIndex = numSteps - 1;

    for (auto& lutElem : lut_) {
        int currentIndex = lutElem.second.first;
        if (currentIndex == smallestIndex) {
            lutElem.second.second = minValue;
        }
        if (currentIndex == largestIndex) {
            lutElem.second.second = maxValue;
        }
    }

    for (auto& inverseLutElem : inverseLut_) {
        int currentIndex = inverseLutElem.first;
        if (currentIndex == smallestIndex) {
            inverseLutElem.second = minValue;
        }
        if (currentIndex == largestIndex) {
            inverseLutElem.second = maxValue;
        }
    }
}

}  // namespace calq

/**
 * @file uniform-quantizer.cc
 */

#include "uniform-quantizer.h"
#include <cassert>
#include <cmath>
#include <queue>
#include <utility>

namespace calq {

UniformQuantizer::UniformQuantizer(const int minValue, const int maxValue, const int numSteps) : Quantizer() {
    assert(minValue <= maxValue);
    assert(numSteps > 1);

    // Compute the step size
    auto stepSize =
        static_cast<int>(floor((static_cast<double>(maxValue - minValue)) / (static_cast<double>(numSteps))));

    // Compute the borders and the representative values
    std::queue<int> borders;
    std::queue<int> reconstructionValues;
    int newBorder = minValue;
    borders.push(minValue);
    reconstructionValues.push(minValue + static_cast<int>(round(static_cast<double>(stepSize) / 2.0)));
    for (int i = 0; i < (numSteps - 1); i++) {
        newBorder += stepSize;
        borders.push(newBorder);
        reconstructionValues.push(newBorder + static_cast<int>(round((static_cast<double>(stepSize) / 2.0))));
    }
    borders.push(maxValue);

    // Fill the quantization table
    borders.pop();
    int currentIndex = 0;
    int currentReconstructionValue = reconstructionValues.front();
    int currentBorder = borders.front();
    for (int value = minValue; value <= maxValue; ++value) {
        if (value > currentBorder) {
            currentIndex++;
            reconstructionValues.pop();
            borders.pop();
            currentReconstructionValue = reconstructionValues.front();
            currentBorder = borders.front();
        }
        std::pair<int, int> curr(currentIndex, currentReconstructionValue);
        lut_.insert(std::pair<int, std::pair<int, int>>(value, curr));
        inverseLut_.insert(curr);
    }
}

}  // namespace calq

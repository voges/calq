#include "calq/uniform_quantizer.h"

#include <math.h>

#include <queue>
#include <utility>

#include "calq/exceptions.h"

namespace calq {

UniformQuantizer::UniformQuantizer(const int &valueMin, const int &valueMax, const int &nrSteps)
    : Quantizer() {
    if ((valueMin > valueMax) || (nrSteps <= 1)) {
        throwErrorException("Error in quantizer initialization");
    }

    // Compute the step size
    int stepSize = static_cast<int>(floor(static_cast<double>(valueMax - valueMin) / static_cast<double>(nrSteps)));

    // Compute the borders and the representative values
    std::queue<int> borders;
    std::queue<int> reconstructionValues;
    int newBorder = valueMin;
    borders.push(valueMin);
    reconstructionValues.push(
        valueMin + static_cast<int>(round(static_cast<double>(stepSize) / static_cast<double>(2)))
    );
    for (int i = 0; i < (nrSteps-1); i++) {
        newBorder += stepSize;
        borders.push(newBorder);
        reconstructionValues.push(
            newBorder + static_cast<int>(round(static_cast<double>(stepSize) / static_cast<double>(2)))
        );
    }
    borders.push(valueMax);

    // Fill the quantization table
    borders.pop();
    int currentIndex = 0;
    int currentReconstructionValue = reconstructionValues.front();
    int currentBorder = borders.front();
    for (int value = valueMin; value <= valueMax; ++value) {
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

UniformQuantizer::~UniformQuantizer(void) {}

}  // namespace calq

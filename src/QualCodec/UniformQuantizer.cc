/** @file UniformQuantizer.cc
 *  @brief This file contains the implementation of the UniformQuantizer class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#include "QualCodec/UniformQuantizer.h"

#include <math.h>

#include <queue>

#include "Common/Exceptions.h"

namespace calq {

UniformQuantizer::UniformQuantizer(const int &valueMin, const int &valueMax, const unsigned int &nrSteps)
    : Quantizer()
{
    if ((valueMin > valueMax) || (nrSteps <= 1)) {
        throwErrorException("Error in quantizer initialization");
    }

    // Compute the step size
    double stepSize = (valueMax - valueMin) / nrSteps;

    // Compute the borders and the representative values
    std::queue<double> borders;
    std::queue<int> reconstructionValues;
    double newBorder = valueMin;
    borders.push(valueMin);
    reconstructionValues.push(valueMin + (int)round(stepSize/2));
    for (unsigned i = 0; i < (nrSteps-1); i++) {
        newBorder += stepSize;
        borders.push(newBorder);
        reconstructionValues.push((int)newBorder + (int)round(stepSize/2));
    }
    borders.push(valueMax);

    // Fill the quantization table
    borders.pop();
    int currentIndex = 0;
    int currentReconstructionValue = reconstructionValues.front();
    double currentBorder = borders.front();
    for (int value = valueMin; value <= valueMax; value++) {
        if (value > currentBorder) {
            currentIndex++;
            reconstructionValues.pop();
            borders.pop();
            currentReconstructionValue = reconstructionValues.front();
            currentBorder = borders.front();
        }
        std::pair<int,int> curr(currentIndex, currentReconstructionValue);
        lut_.insert(std::pair<int,std::pair<int,int>>(value, curr));
        inverseLut_.insert(curr);
    }
}

UniformQuantizer::~UniformQuantizer(void) {}

} // namespace calq


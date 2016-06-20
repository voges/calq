/** @file UniformQuantizer.cc
 *  @brief This file contains the implementation of the UniformQuantizer class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: what (who)
 */

#include "Codecs/UniformQuantizer.h"
#include "Exceptions.h"
#include <cmath>
#include <queue>
#include <vector>

UniformQuantizer::UniformQuantizer(const int &minimumValue, 
                                   const int &maximumValue, 
                                   const unsigned int &numberOfSteps)
{
    // sanity checks
    if ((minimumValue >= maximumValue) || (numberOfSteps <= 1)) {
        throwErrorException("Error in quantizer initialization");
    }

    // compute the step size
    double stepSize = (maximumValue - minimumValue) / numberOfSteps;

    // compute the borders and the representative values
    std::queue<double> borders;
    std::queue<int> representativeValues;
    double newBorder = minimumValue;

    borders.push(minimumValue);
    representativeValues.push(minimumValue + round(stepSize/2));
    for (unsigned i = 0; i < numberOfSteps-1; i++) {
        newBorder += stepSize;
        borders.push(newBorder);
        representativeValues.push(newBorder + round(stepSize/2));
    }
    borders.push(maximumValue);

    // fill the quantization table
    borders.pop();
    int currentIndex = 0;
    int currentRepresentativeValue = representativeValues.front();
    double currentBorder = borders.front();
    for (int value = minimumValue; value <= maximumValue; value++) {
        if (value > currentBorder) {
            currentIndex++;
            representativeValues.pop();
            borders.pop();
            currentRepresentativeValue = representativeValues.front();
            currentBorder = borders.front();
        }
        std::pair<int,int> curr(currentIndex, currentRepresentativeValue);
        lut.insert(std::pair<int,std::pair<int,int>>(value, curr));
    }
}

UniformQuantizer::~UniformQuantizer(void)
{
    // empty
}

int UniformQuantizer::getIndex(const int &value)
{
    if (lut.find(value) == lut.end()) {
        throwErrorException("Value out of range for quantizer");
    }

    return lut[value].first;
}

int UniformQuantizer::getRepresentativeValue(const int &value)
{
    if (lut.find(value) == lut.end()) {
        throwErrorException("Value out of range for quantizer");
    }

    return lut[value].second;
}

void UniformQuantizer::print(void)
{
    for (auto const &lutEntry : lut) {
        std::cout << lutEntry.first << ": ";
        std::cout << lutEntry.second.first << ",";
        std::cout << lutEntry.second.second << std::endl;
    }
}


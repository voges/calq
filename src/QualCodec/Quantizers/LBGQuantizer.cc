/** @file LBGQuantizer.cc
 *  @brief This file contains the implementation of the LBGQuantizer class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#include "QualCodec/Quantizers/LBGQuantizer.h"

#include <algorithm>
#include <cmath>
#include <queue>
#include <sstream>
#include <vector>

#include "Common/Exceptions.h"
#include "Common/log.h"

namespace calq {

LBGQuantizer::LBGQuantizer(const int &valueMax,
                           const int &valueMin,
                           const int &nrSteps,
                           const std::string &sampleValues)
    : Quantizer(),
      centers_(),
      valueMax_(valueMax),
      valueMin_(valueMin)
{
    if ((valueMin_ >= valueMax_) || (nrSteps <= 1)) {
        throwErrorException("Error in quantizer initialization");
    }

    // Compute the step size
    double stepSize = (valueMax_ - valueMin_) / nrSteps;

    // Compute the borders and the representative values
    std::queue<double> borders;
    std::queue<int> reconstructionValues;
    double newBorder = valueMin_;
    borders.push(valueMin_);
    reconstructionValues.push(valueMin_ + round(stepSize/2));
    centers_.push_back(newBorder + round(stepSize/2));
    for (int i = 0; i < (nrSteps-1); i++) {
        newBorder += stepSize;
        borders.push(newBorder);
        reconstructionValues.push(newBorder + round(stepSize/2));
        centers_.push_back(newBorder + round(stepSize/2));
    }
    borders.push(valueMax_);

    // Fill the quantization table
    borders.pop();
    int currentIndex = 0;
    int currentReconstructionValue = reconstructionValues.front();
    double currentBorder = borders.front();
    for (int value = valueMin_; value <= valueMax_; value++) {
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

    // Train the quantizer
    train(sampleValues);
}

LBGQuantizer::~LBGQuantizer(void) {}

static double dist(const double &a, const double &b)
{
    if (a > b) {
        return a-b;
    } else {
        return b-a;
    }
}

static double nearest(const std::tuple<int, size_t, double> &elem, const std::vector<double> &centers)
{
    double distance = std::numeric_limits<double>::max();
    double newCenter = std::get<2>(elem);

    for (double const &center : centers) {
        double distanceToCenter = dist(center, (double)(std::get<0>(elem)));
        if (distanceToCenter < distance) {
            distance = distanceToCenter;
            newCenter = center;
        }
    }

    return newCenter;
}

void LBGQuantizer::train(const std::string &sampleValues)
{
    std::vector<std::tuple<int, size_t, double>> sampleValueDistribution;

    // Initialize QV distribution
    CALQ_LOG("Initializing QV distribution");
    for(int value = valueMin_; value <= valueMax_; ++value) {
        sampleValueDistribution.emplace_back(std::make_tuple(value, 0, (double)valueToReconstructionValue(value)));
    }

    // Generate cumulated absolute frequency distribution
    CALQ_LOG("Generating cumulated absolute frequency distribution");
    for (const char &sampleValue : sampleValues) {
        for (auto &elem : sampleValueDistribution){
            if(sampleValue == std::get<0>(elem)) {
                std::get<1>(elem)++;
                break;
            }
        }
    }

    // Do the clustering
    CALQ_LOG("Clustering");
    bool changed = false;
    size_t nrIterations = 0;

    do {
        changed = false;

        // Compute new cluster centers
        for (double &center : centers_) {
            double newCenter = 0;
            size_t freq = 0;

            for(auto const &elem : sampleValueDistribution) {
                if (std::get<2>(elem) == center) {
                    newCenter += (double)(std::get<0>(elem)) * (double)(std::get<1>(elem));
                    freq += std::get<1>(elem);
                }
            }

            if (freq != 0) {
                center = (double)newCenter / (double)freq;
            }
        }

        // Assign each QV to new cluster center
        for (auto &elem : sampleValueDistribution){
            double center = nearest(elem, centers_);
            if (std::get<2>(elem) != center) {
                std::get<2>(elem) = center;
                changed = true;
            }
        }

        nrIterations++;
    } while (changed == true);

    CALQ_LOG("Finished clustering after %zu iterations", nrIterations);

    // Fill LUT and inverse LUT
    CALQ_LOG("Filling LUT and inverse LUT");
    int index = 0;
    int reconstructionValuePrev = 0;
    bool first = true;

    for (auto const &elem : sampleValueDistribution) {
        int value = std::get<0>(elem);

        auto lutIt = lut_.find(value);
        if (lutIt == lut_.end()) {
            throwErrorException("Did not find value");
        }

        int reconstructionValue = round(std::get<2>(elem));
        if (first == true) {
            reconstructionValuePrev = reconstructionValue;
            first = false;
        }

        if (reconstructionValue != reconstructionValuePrev) {
            index++;
        }

        auto inverseLutIt = inverseLut_.find(index);
        if (inverseLutIt == inverseLut_.end()) {
            throwErrorException("Did not find value");
        }

        // Fill LUT
        lutIt->second.first = index;
        lutIt->second.second = reconstructionValue;

        // Fill inverse LUT
        inverseLutIt->second = reconstructionValue;

        reconstructionValuePrev = reconstructionValue;
    }
}

} // namespace calq


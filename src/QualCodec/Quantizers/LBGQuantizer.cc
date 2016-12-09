/** @file LBGQuantizer.cc
 *  @brief This file contains the implementation of the LBGQuantizer class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#include "QualCodec/Quantizers/LBGQuantizer.h"

#include <random>

#include "Common/Exceptions.h"
#include "Common/log.h"

namespace calq {

LBGQuantizer::LBGQuantizer(const int &k,
                           const std::map<int, size_t> &sampleValueDistribution)
    : Quantizer()
{
    if (k <= 1) {
        throwErrorException("More than 1 cluster centers required");
    }

    if (sampleValueDistribution.empty() == true) {
        throwErrorException("sampleValueDistribution is empty");
    }

    int sampleValueMax = sampleValueDistribution.end()->first;
    int sampleValueMin = sampleValueDistribution.begin()->first;

    // Initialize k random cluster centers
    //   i -> (centerValue,centerValueFrequency,centerSampleValues)
    std::map<int, std::tuple<int, size_t, std::vector<int>>> centers; 
    std::random_device seeder;
    std::mt19937 rng(seeder());
    std::uniform_int_distribution<int> gen(sampleValueMin, sampleValueMax);
    for (int i = 0; i < k; ++i) {
        for (;;) {
            int centerValue = gen(rng);
            std::cout << "rng: " << centerValue << std::endl;
            if (centers.find(centerValue) == centers.end()) {
                size_t centerValueFrequency =  sampleValueDistribution.at(centerValue);
                std::vector<int> centerSampleValues;
                std::pair<int, std::tuple<int, size_t, std::vector<int>>> center;
                center.first = i;
                std::get<0>(center.second) = centerValue;
                std::get<1>(center.second) = centerValueFrequency;
                std::get<2>(center.second) = centerSampleValues;
                centers.insert(center);
                break;
            }
        }
    }

    // Do the clustering
    CALQ_LOG("Clustering");
//     bool changed = false;
//     size_t nrIterations = 0;
// 
//     do {
//         changed = false;
// 
//         // Assign sample values to cluster centers
//         for(int sampleValue = sampleValueMin_; sampleValue <= sampleValueMax_; ++sampleValue) {
//             double smallestDistance = std::numeric_limits<double>::max();
//             int nearestCenterValue = -1;
//             size_t sampleValueFrequency = sampleValueDistribution.at(sampleValue);
//             for (auto const &center : centers) {
//                 int centerValue = center.first;
//                 size_t centerValueFrequency = center.second.first;
//                 double distance = (sampleValue-centerValue)*(sampleValue-centerValue)
//                                  +(sampleValueFrequency-centerValueFrequency)*(sampleValueFrequency-centerValueFrequency);
//                 if (distance < smallestDistance) {
//                     smallestDistance = distance;
//                     nearestCenterValue = centerValue;
//                 }
//             }
//             auto nearestCenter = centers.at(nearestCenterValue);
//             nearestCenter.second.push_back(sampleValue);
//         }
// 
//         // Compute new cluster centers
//         for (auto &center : centers) {
//             auto sampleValues = center.second.second;
//             double newCenterValue = 0.0;
// 
//             for (auto const &sampleValue : sampleValues) {
//                 newCenterValue += (double)sampleValue;
//             }
//             newCenterValue /= sampleValues.size();
// 
//             size_t newCenterValueFrequency = sampleValueDistribution.at(newCenterValue);
// 
//             center.first = newCenterValue;
//             center.second.first = newCenterValueFrequency;
//             center.second.second.clear();
//         }
// 
//         nrIterations++;
//     } while (changed == true);
// 
//     CALQ_LOG("Finished clustering after %zu iterations", nrIterations);

    // Fill LUT and inverse LUT
//     CALQ_LOG("Filling LUT and inverse LUT");
//     int index = 0;
//     int reconstructionValuePrev = 0;
//     bool first = true;
// 
//     for (auto const &elem : sampleValueDistribution) {
//         int value = std::get<0>(elem);
// 
//         auto lutIt = lut_.find(value);
//         if (lutIt == lut_.end()) {
//             throwErrorException("Did not find value");
//         }
// 
//         int reconstructionValue = round(std::get<2>(elem));
//         if (first == true) {
//             reconstructionValuePrev = reconstructionValue;
//             first = false;
//         }
// 
//         if (reconstructionValue != reconstructionValuePrev) {
//             index++;
//         }
// 
//         auto inverseLutIt = inverseLut_.find(index);
//         if (inverseLutIt == inverseLut_.end()) {
//             throwErrorException("Did not find value");
//         }
// 
//         // Fill LUT
//         lutIt->second.first = index;
//         lutIt->second.second = reconstructionValue;
// 
//         // Fill inverse LUT
//         inverseLutIt->second = reconstructionValue;
// 
//         reconstructionValuePrev = reconstructionValue;
//     }
}

LBGQuantizer::~LBGQuantizer(void) {}

} // namespace calq


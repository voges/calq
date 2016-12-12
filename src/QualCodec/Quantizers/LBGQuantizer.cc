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

static void initCenters(const int &sampleValueMin,
                        const int &sampleValueMax,
                        const std::map<int, size_t> &sampleValueDistribution,
                        const int &k,
                        std::map<int, std::tuple<int, size_t, std::vector<int>>> *centers)
{
    CALQ_LOG("Initializing %d random cluster centers", k);

    for (auto &center : *centers) {
        std::get<2>(center.second).clear();
        //center.second.clear();
    }
    centers->clear();

    std::random_device seeder;
    std::mt19937 rng(seeder());
    std::uniform_int_distribution<int> gen(sampleValueMin, sampleValueMax);
    std::map<int, int> centerValues;

    for (int i = 0; i < k; ++i) {
        bool proceed = false;
        do {
            proceed = false;
            int centerValue = gen(rng);
            if (centerValues.find(centerValue) == centerValues.end()) {
                size_t centerValueFrequency =  sampleValueDistribution.at(centerValue);
                std::vector<int> centerSampleValues;
                std::pair<int, std::tuple<int, size_t, std::vector<int>>> center;
                center.first = i;
                std::get<0>(center.second) = centerValue;
                std::get<1>(center.second) = centerValueFrequency;
                std::get<2>(center.second) = centerSampleValues;
                centers->insert(center);
                centerValues.insert(std::pair<int, int>(centerValue, i));
                proceed = true;
            }
            printf("generating...\n");
        } while (proceed == false);
    }

    // Bubble sort
//     bool swapped = false;
//     do {
//         swapped = false;
//         for (int i = 0; i < k-1; ++i) {
//             int first = std::get<0>(centers->at(i));
//             int second = std::get<0>(centers->at(i+1));
//             if (first > second) {
//                 swapped = true;
// 
//                 std::tuple<int, size_t, std::vector<int>> tmp;
// 
//                 std::get<0>(tmp) = std::get<0>(centers->at(i));
//                 std::get<1>(tmp) = std::get<1>(centers->at(i));
//                 std::get<2>(tmp) = std::get<2>(centers->at(i));
// 
//                 std::get<0>(centers->at(i)) = std::get<0>(centers->at(i+1));
//                 std::get<1>(centers->at(i)) = std::get<1>(centers->at(i+1));
//                 std::get<2>(centers->at(i)) = std::get<2>(centers->at(i+1));
// 
//                 std::get<0>(centers->at(i+1)) = std::get<0>(tmp);
//                 std::get<1>(centers->at(i+1)) = std::get<1>(tmp);
//                 std::get<2>(centers->at(i+1)) = std::get<2>(centers->at(i+1));
//             }
//         }
//     } while (swapped == true);
}

static void printCenters(const std::map<int, std::tuple<int, size_t, std::vector<int>>> &centers)
{
    for (auto const &center : centers) {
        printf("%d -> %d, %zu, (", center.first, std::get<0>(center.second), std::get<1>(center.second));
        for (auto const &sampleValue : std::get<2>(center.second)) {
            printf(" %d ", sampleValue);
        }
        printf(")\n");
    }
}

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

    int sampleValueMax = sampleValueDistribution.rbegin()->first;
    int sampleValueMin = sampleValueDistribution.begin()->first;

    // Initialize k random cluster centers
    //   i -> (centerValue,centerValueFrequency,centerSampleValues)
    std::map<int, std::tuple<int, size_t, std::vector<int>>> centers; 
    initCenters(sampleValueMin, sampleValueMax, sampleValueDistribution, k, &centers);

    // Do the clustering
    int changed = 0;
    size_t nrIterations = 0;
    size_t nrIterationsMax = 100;
    size_t nrReinits = 0;
    CALQ_LOG("Clustering (max #iterations: %zu)", nrIterationsMax);

    do {
        changed = 0;

        // Assign sample values to cluster centers
        for(int sampleValue = sampleValueMin; sampleValue <= sampleValueMax; ++sampleValue) {

            // Search nearest center for sampleValue
            double smallestDistance = std::numeric_limits<double>::max();
            int nearestCenter = -1;
            size_t sampleValueFrequency = sampleValueDistribution.at(sampleValue);
            for (auto const &center : centers) {
                int centerValue = std::get<0>(center.second);
                size_t centerValueFrequency = std::get<1>(center.second);
                double distance = (sampleValue-centerValue)*(sampleValue-centerValue)
                                 +(sampleValueFrequency-centerValueFrequency)*(sampleValueFrequency-centerValueFrequency);
                if (distance < smallestDistance) {
                    smallestDistance = distance;
                    nearestCenter = center.first;
                }
            }

            // Add sampleValue to nearest center
            std::get<2>(centers.at(nearestCenter)).push_back(sampleValue);
        }

        // Compute new cluster centers
        for (auto &center : centers) {
            auto centerSampleValues = std::get<2>(center.second);
            double newCenterValue = 0.0;

            for (auto const &centerSampleValue : centerSampleValues) {
                newCenterValue += (double)centerSampleValue;
            }
            newCenterValue /= centerSampleValues.size();
            newCenterValue = round(newCenterValue);

            size_t newCenterValueFrequency = sampleValueDistribution.at(newCenterValue);

            int oldCenterValue = std::get<0>(center.second);
            if (oldCenterValue != newCenterValue) {
                changed++;
            }

            std::get<0>(center.second) = newCenterValue;
            std::get<1>(center.second) = newCenterValueFrequency;
            std::get<2>(center.second).clear();
        }

        // Check if two cluster centers merged; if so, re-initialize randomly
        bool reinit = false;
        do {
            reinit = false;
            for (auto const &center : centers) {
                int centerValue = std::get<0>(center.second);
                int identicalCenterValues = 0;
                for (auto const &centerToCompare : centers) {
                    int centerValueToCompare = std::get<0>(centerToCompare.second);
                    if (centerValue == centerValueToCompare) {
                        identicalCenterValues++;
                    }
                    // If we find the same cluster center twice, re-initialize!
                    if (identicalCenterValues > 1) {
                        reinit = true;
                    }
                }
            }

            if (reinit == true) {
                initCenters(sampleValueMin, sampleValueMax, sampleValueDistribution, k, &centers);
                nrReinits++;
            }
        } while (reinit == true);

        CALQ_LOG("Iteration %zu: %d cluster centers changed", nrIterations, changed);
        nrIterations++;
    } while (changed > 1 && nrIterations < nrIterationsMax);

    CALQ_LOG("Finished clustering after %zu iteration(s) (did %zu re-inits)", nrIterations, nrReinits);

    // Fill LUT and inverse LUT
    for(int sampleValue = sampleValueMin; sampleValue <= sampleValueMax; ++sampleValue) {

        // Find nearest center for this sampleValue
        double smallestDistance = std::numeric_limits<double>::max();
        int nearestCenter = -1;
        size_t sampleValueFrequency = sampleValueDistribution.at(sampleValue);
        for (auto const &center : centers) {
            int centerValue = std::get<0>(center.second);
            size_t centerValueFrequency = std::get<1>(center.second);
            double distance = (sampleValue-centerValue)*(sampleValue-centerValue)
                              +(sampleValueFrequency-centerValueFrequency)*(sampleValueFrequency-centerValueFrequency);
            if (distance < smallestDistance) {
                smallestDistance = distance;
                nearestCenter = center.first;
            }
        }

        int value = sampleValue;
        int index = nearestCenter;
        int reconstructionValue = std::get<0>(centers.at(nearestCenter));

        // Insert in LUT and inverse LUT
        std::pair<int, int> elem(index, reconstructionValue);
        lut_.insert(std::pair<int, std::pair<int, int>>(value, elem));
        inverseLut_.insert(elem);
    }

    // Regularize LUT
//     auto previousLutEntry = lut_.at(0);
// 
//     for (auto &lutEntry : lut_) {
//         int currentIndex = lutEntry.second.first;
//         int previousIndex = previousLutEntry.first;
// 
//         if (currentIndex < previousIndex) {
//             lutEntry.second.first = previousLutEntry.first;
//             lutEntry.second.second = previousLutEntry.second;
//         }
//         previousLutEntry = lutEntry.second;
//     }
}

LBGQuantizer::~LBGQuantizer(void) {}

} // namespace calq


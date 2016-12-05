/** @file KMEANQuantizers.h
 *  @brief This file contains the definitions of the KMEANQuantizer class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef KMEANQUANTIZER_H
#define KMEANQUANZIZER_H

#include <map>
#include <string>
#include <vector>
#include <fstream>

/** @brief Class: KMEANQuantizer
 *
 *  The KMEANQuantizer class quantizes integer values.
 */
class KMEANQuantizer {
public:
    KMEANQuantizer(const int &minimumValue,
                     const int &maximumValue,
                     const unsigned int &numberOfSteps);
    ~KMEANQuantizer(void);

    int valueToIndex(const int &value);
    double indexToReconstructionValue(const int &index);
    double valueToReconstructionValue(const int &value);

    void print(void) const;
    void train(const std::string &samplePoints,const int &blockId, const int &quantId);

private:
    std::map<int,std::pair<int,double>> lut; // value->(index,reconstructionValue)
    std::map<int,double> inverseLut; // index->reconstructionValue
    std::vector<int> samplePoints;
    std::vector<double> centers;
    const size_t minValue;
    const size_t maxValue;
};

#endif // KMEANQUANTIZER_H


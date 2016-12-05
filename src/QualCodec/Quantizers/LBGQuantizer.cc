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

KMEANQuantizer::KMEANQuantizer(const int &minimumValue,
                                const int &maximumValue,
                                const unsigned int &numberOfSteps)
    : lut()
    , inverseLut()
    , samplePoints()
    , centers()
    , minValue(minimumValue)
    , maxValue(maximumValue)
{
    // Sanity checks
    if ((minValue >= maxValue) || (numberOfSteps <= 1)) {
        throwErrorException("Error in quantizer initialization");
    }

    // Compute the step size
    double stepSize = (maxValue - minValue) / numberOfSteps;

    // Compute the borders and the representative values
    std::queue<double> borders;
    std::queue<int> reconstructionValues;
    double newBorder = minValue;

    borders.push(minValue);
    reconstructionValues.push(minValue + round(stepSize/2));
    centers.push_back(newBorder + round(stepSize/2));
    for (size_t i = 0; i < numberOfSteps-1; i++) {
        newBorder += stepSize;
        borders.push(newBorder);
        reconstructionValues.push(newBorder + round(stepSize/2));
        centers.push_back(newBorder + round(stepSize/2));
    }
    borders.push(maxValue);

    // Fill the quantization table
    borders.pop();
    int currentIndex = 0;
    int currentReconstructionValue = reconstructionValues.front();
    double currentBorder = borders.front();
    for (size_t value = minValue; value <= maxValue; value++) {
        if (value > currentBorder) {
            currentIndex++;
            reconstructionValues.pop();
            borders.pop();
            currentReconstructionValue = reconstructionValues.front();
            currentBorder = borders.front();
        }
        std::pair<int,int> curr(currentIndex, currentReconstructionValue);
        lut.insert(std::pair<int,std::pair<int,int>>(value, curr));
        inverseLut.insert(curr);
    }
}

KMEANQuantizer::~KMEANQuantizer(void)
{
    // empty
}

int KMEANQuantizer::valueToIndex(const int &value)
{
    if (lut.find(value) == lut.end()) {
        throwErrorException("Value out of range for quantizer");
    }

    return lut[value].first;
}

double KMEANQuantizer::indexToReconstructionValue(const int &index)
{
    if (inverseLut.find(index) == inverseLut.end()) {
        throwErrorException("Quantization index not found");
    }

    return inverseLut[index];
}

double KMEANQuantizer::valueToReconstructionValue(const int &value)
{
    if (lut.find(value) == lut.end()) {
        throwErrorException("Value out of range for quantizer");
    }

    return lut[value].second;
}

void KMEANQuantizer::print(void) const
{
    std::cout << "LUT:" << std::endl;
    for (auto const &lutEntry : lut) {
        std::cout << "  " << lutEntry.first << ": ";
        std::cout << lutEntry.second.first << ",";
        std::cout << lutEntry.second.second << std::endl;
    }

    std::cout << "Inverse LUT:" << std::endl;
    for (auto const &inverseLutEntry : inverseLut) {
        std::cout << "  " << inverseLutEntry.first << ": ";
        std::cout << inverseLutEntry.second << std::endl;
    }
}


struct compclass {
  bool operator() (int i,int j) { return (i<j);}
} compobj;


double dist(double a, double b){
    if (a>b){
        return a-b;
    } else{
        return b-a;
    }
}

double nearest(std::tuple<size_t,size_t,double> elem, const std::vector<double> &centers){
    double distance = std::numeric_limits<double>::max();
    double newCenter = std::get<2>(elem);
    for (auto &center : centers){
        if(dist(center,std::get<0>(elem)) < distance){
            distance = dist(center,std::get<0>(elem));
            newCenter = center;
        }
    }
    return newCenter;
}


void KMEANQuantizer::train(const std::string &samplePoints,const int &blockId, const int &quantId)
{
    //std::map<int,std::vector<int>> ClusterIDSamplePoints;
    //std::vector<std::pair<int,double>> samplePointsVec;
    std::vector<std::tuple<size_t,size_t,double>> sPtupleVec; // dist
    
    std::ofstream ofs;
    std::ostringstream os;
    os << "Block_" << blockId << "_ID_" << quantId << ".csv";
    std::string fileName = os.str();
    ofs.open(fileName,std::ofstream::out);
    ofs << "Iteration,QV,Anzahl,Center" << std::endl;
    
    // Initialize QV distribution
    for(size_t i = minValue; i<=maxValue; ++i){
        sPtupleVec.emplace_back(std::make_tuple(i,0,valueToReconstructionValue(i)));
    }

    // Iterate through all QVs and generate cumulated absolute frequency distribution
    for(const char& c: samplePoints){
        for(auto &elem: sPtupleVec){
            if(c == std::get<0>(elem)){
                std::get<1>(elem)++;
                break;
            }
        }
    }

    // ^^^^ std::tuple<int, int, double>
    //      std::tuple<qv, cumFreq, center>

    size_t changed = 0;
    size_t iteration = 0;
    do{ 
        changed = 0;       

        // Compute new cluster centers
        for(auto &center: centers){
            double newCenter = 0;
            size_t ctr = 0;

            for(auto const &pt: sPtupleVec){
                if(std::get<2>(pt) == center){
                    newCenter += std::get<0>(pt) * std::get<1>(pt);
                    ctr += std::get<1>(pt);
                }
            }
            if (ctr != 0){
                center = (double)newCenter/(double)ctr;
            }
        }

        // Assign each QV to new cluster center
        for(auto &pt: sPtupleVec){
            ofs << iteration << "," << std::get<0>(pt) << "," << std::get<1>(pt) << "," << std::get<2>(pt) << std::endl;
            double center = nearest(pt, centers);
            if(std::get<2>(pt) != center){
                std::get<2>(pt) = center;
                changed++;
            }
        }
        iteration++;
    }while(changed!=0);

    std::map<double, int> reverseInverseLut;
    for(auto &elem: inverseLut){
        elem.second = centers[elem.first];
        reverseInverseLut.emplace(elem.second,elem.first);
    }

    for (auto &elem: lut){
        for(auto &pt: sPtupleVec){
            if(elem.first == std::get<0>(pt)){
                elem.second.first = reverseInverseLut[std::get<2>(pt)];
                elem.second.second = std::get<2>(pt);
                break;
            }
        }
    }
    ofs.close();



    CALQ_LOG(" Quantizer %d for Block %d trained.",quantId, blockId);
}

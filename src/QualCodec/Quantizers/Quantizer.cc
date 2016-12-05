/** @file Quantizer.cc
 *  @brief This file contains the implementation of the Quantizer base class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#include "QualCodec/Quantizers/Quantizer.h"

#include "Common/Exceptions.h"

namespace calq {

Quantizer::Quantizer(void) : lut_() , inverseLut_() {}

Quantizer::Quantizer(const std::map<int, int> &inverseLut)
    : lut_(),
      inverseLut_(inverseLut) {}

Quantizer::~Quantizer(void) {}

int Quantizer::valueToIndex(const int &value) const
{
    if (lut_.find(value) == lut_.end()) {
        throwErrorException("Value out of range");
    }
    return lut_.at(value).first;
}

int Quantizer::indexToReconstructionValue(const int &index) const
{
    if (inverseLut_.find(index) == inverseLut_.end()) {
        throwErrorException("Index not found");
    }
    return inverseLut_.at(index);
}

int Quantizer::valueToReconstructionValue(const int &value) const
{
    if (lut_.find(value) == lut_.end()) {
        throwErrorException("Value out of range");
    }

    return lut_.at(value).second;
}

const std::map<int, int> & Quantizer::inverseLut(void) const
{
    return inverseLut_;
}

void Quantizer::print(void) const
{
    std::cout << "LUT:" << std::endl;
    for (auto const &lutEntry : lut_) {
        std::cout << "  " << lutEntry.first << ": ";
        std::cout << lutEntry.second.first << ",";
        std::cout << lutEntry.second.second << std::endl;
    }

    std::cout << "Inverse LUT:" << std::endl;
    for (auto const &inverseLutEntry : inverseLut_) {
        std::cout << "  " << inverseLutEntry.first << ": ";
        std::cout << inverseLutEntry.second << std::endl;
    }
}

} // namespace calq


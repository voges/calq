/** @file Quantizer.cc
 *  @brief This file contains the implementation of the Quantizer base class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "QualCodec/Quantizer.h"
#include "Common/Exceptions.h"

cq::Quantizer::Quantizer(void)
    : lut_()
    , inverseLut_()
{
    // empty
}

cq::Quantizer::~Quantizer(void)
{
    // empty
}

int cq::Quantizer::valueToIndex(const int &value) const
{
    if (lut_.find(value) == lut_.end()) {
        throwErrorException("Value out of range");
    }
    return lut_.at(value).first;
}

int cq::Quantizer::indexToReconstructionValue(const int &index) const
{
    if (inverseLut_.find(index) == inverseLut_.end()) {
        throwErrorException("Index not found");
    }
    return inverseLut_.at(index);
}

int cq::Quantizer::valueToReconstructionValue(const int &value) const
{
    if (lut_.find(value) == lut_.end()) {
        throwErrorException("Value out of range");
    }

    return lut_.at(value).second;
}

const std::map<int, int> & cq::Quantizer::inverseLut(void) const
{
    return inverseLut_;
}

void cq::Quantizer::print(void) const
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


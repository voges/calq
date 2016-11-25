/** @file UniformQuantizer.h
 *  @brief This file contains the definitions of the UniformQuantizer class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef CQ_UNIFORMQUANTIZER_H
#define CQ_UNIFORMQUANTIZER_H

#include <map>
#include <string>
#include <vector>

namespace cq {

class UniformQuantizer {
public:
    UniformQuantizer(const int &minimumValue,
                     const int &maximumValue,
                     const unsigned int &numberOfSteps);
    ~UniformQuantizer(void);

    int valueToIndex(const int &value);
    int indexToReconstructionValue(const int &index);
    int valueToReconstructionValue(const int &value);

    void print(void) const;

private:
    std::map<int,std::pair<int,int>> lut; // value->(index,reconstructionValue)
    std::map<int,int> inverseLut; // index->reconstructionValue
};

}

#endif // CQ_UNIFORMQUANTIZER_H


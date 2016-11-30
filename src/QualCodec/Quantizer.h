/** @file Quantizer.h
 *  @brief This file contains the definition of the Quantizer base class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#ifndef CALQ_QUALCODEC_QUANTIZER_H_
#define CALQ_QUALCODEC_QUANTIZER_H_

#include <map>

namespace calq {

class Quantizer {
public:
    Quantizer(void);
    virtual ~Quantizer(void);

    int valueToIndex(const int &value) const;
    int indexToReconstructionValue(const int &index) const;
    int valueToReconstructionValue(const int &value) const;

    const std::map<int,int> & inverseLut(void) const;

    void print(void) const;

protected:
    std::map<int,std::pair<int,int>> lut_; // value->(index,reconstructionValue)
    std::map<int,int> inverseLut_; // index->reconstructionValue
};

} // namespace calq

#endif // CALQ_QUALCODEC_QUANTIZER_H_


/** @file UniformQuantizer.h
 *  @brief This file contains the definitions of the UniformQuantizer class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#ifndef CALQ_QUALCODEC_UNIFORMQUANTIZER_H_
#define CALQ_QUALCODEC_UNIFORMQUANTIZER_H_

#include "QualCodec/Quantizer.h"

namespace calq {

class UniformQuantizer : public Quantizer {
public:
    UniformQuantizer(const int &valueMax, const int &valueMin, const int &nrSteps);
    ~UniformQuantizer(void);
};

} // namespace calq

#endif // CALQ_QUALCODEC_UNIFORMQUANTIZER_H_


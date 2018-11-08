/** @file OneToOneMappingQuantizer.h
 *  @brief This file contains the definitions of the OneToOneMappingQuantizer class.
 */

// Copyright 2015-2017 Leibniz Universitaet Hannover

#ifndef CALQ_QUALCODEC_QUANTIZERS_ONETOONEMAPPINGQUANTIZER_H_
#define CALQ_QUALCODEC_QUANTIZERS_ONETOONEMAPPINGQUANTIZER_H_

#include "Quantizer.h"

namespace calq {

class OneToOneMappingQuantizer : public Quantizer {
 public:
    OneToOneMappingQuantizer(int valueMin, int valueMax);
    ~OneToOneMappingQuantizer(void);
};

}  // namespace calq

#endif  // CALQ_QUALCODEC_QUANTIZERS_ONETOONEMAPPINGQUANTIZER_H_

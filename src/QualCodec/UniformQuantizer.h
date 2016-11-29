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

#include "QualCodec/Quantizer.h"

namespace cq {

class UniformQuantizer : public Quantizer {
public:
    UniformQuantizer(const int &valueMin, const int &valueMax, const unsigned int &nrSteps);
    ~UniformQuantizer(void);
};

}

#endif // CQ_UNIFORMQUANTIZER_H


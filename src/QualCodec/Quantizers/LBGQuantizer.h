/** @file LBGQuantizer.h
 *  @brief This file contains the definitions of the LBGQuantizer class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#ifndef CALQ_QUALCODEC_QUANTIZERS_LBGQUANTIZER_H_
#define CALQ_QUALCODEC_QUANTIZERS_LBGQUANTIZER_H_

#include <map>
#include <string>
#include <vector>

#include "QualCodec/Quantizers/Quantizer.h"

namespace calq {

class LBGQuantizer : public Quantizer {
public:
    LBGQuantizer(const int &k,
                 const std::map<int, size_t> &sampleValueDistribution);
    ~LBGQuantizer(void);
};

} // namespace calq

#endif // CALQ_QUALCODEC_QUANTIZERS_LBGQUANTIZER_H_


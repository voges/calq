/** @file CalqEncoder.h
 *  @brief This file contains the definition of the CalqEncoder class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#ifndef CALQ_CALQENCODER_H_
#define CALQ_CALQENCODER_H_

#include <string>
#include <vector>

#include "Common/Options.h"
#include "IO/CQ/CQFile.h"
#include "IO/SAM/SAMFile.h"

namespace calq {

class CalqEncoder {
public:
    explicit CalqEncoder(const Options &options);
    ~CalqEncoder(void);

    void encode(void);

private:
    size_t blockSize_;
    CQFile cqFile_;
    int polyploidy_;
    int qualityValueMax_;
    int qualityValueMin_;
    int qualityValueOffset_;
    std::vector<std::string> referenceFileNames_;
    SAMFile samFile_;
};

} // namespace calq

#endif // CALQ_CALQENCODER_H_


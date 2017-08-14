/** @file CalqEncoder.h
 *  @brief This file contains the definition of the CalqEncoder class.
 */

// Copyright 2015-2017 Leibniz Universitaet Hannover

#ifndef CALQ_CALQENCODER_H_
#define CALQ_CALQENCODER_H_

#include <string>
#include <vector>

#include "Common/Options.h"
#include "config.h"
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
    std::string inputFileName_;
    int polyploidy_;
    int qualityValueMin_;
    int qualityValueMax_;
    int qualityValueOffset_;
    std::vector<std::string> referenceFileNames_;
    SAMFile samFile_;
};

}  // namespace calq

#endif  // CALQ_CALQENCODER_H_


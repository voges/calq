/** @file CalqEncoder.h
 *  @brief This file contains the definition of the CalqEncoder class.
 */

// Copyright 2015-2017 Leibniz Universitaet Hannover

#ifndef CALQ_CALQENCODER_H_
#define CALQ_CALQENCODER_H_

#include <string>
#include <vector>

#include "Common/Options.h"

namespace calq {

class CQFile;
class FASTAFile;
class SAMFile;

class CalqEncoder {
 public:
    explicit CalqEncoder(const Options &options);
    ~CalqEncoder(void);

    void encode(void);

 private:
    CQFile* cqFile_;
    SAMFile* samFile_;
    FASTAFile* fastaFile_;  // Reference

    const Options options;
};

}  // namespace calq

#endif  // CALQ_CALQENCODER_H_


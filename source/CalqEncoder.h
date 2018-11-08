/** @file CalqEncoder.h
 *  @brief This file contains the definition of the CalqEncoder class.
 */

// Copyright 2015-2017 Leibniz Universitaet Hannover

#ifndef CALQ_CALQENCODER_H_
#define CALQ_CALQENCODER_H_

#include <memory>

#include "Options.h"

namespace calq {

class CQFile;

class FASTAFile;

class SAMFile;

class CalqEncoder {
 public:
    explicit CalqEncoder(const Options &options);
    ~CalqEncoder();
    void encode();

 private:
    std::unique_ptr<CQFile> cqFile_;
    std::unique_ptr<SAMFile> samFile_;
    std::unique_ptr<FASTAFile> fastaFile_;  // Reference

    const Options options;
};

}  // namespace calq

#endif  // CALQ_CALQENCODER_H_

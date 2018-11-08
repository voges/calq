/** @file CalqDecoder.h
 *  @brief This file contains the definition of the CalqDecoder class.
 */

// Copyright 2015-2017 Leibniz Universitaet Hannover

#ifndef CALQ_CALQDECODER_H_
#define CALQ_CALQDECODER_H_

#include <memory>

#include "Options.h"

namespace calq {

class CQFile;

class File;

class SAMFile;

class CalqDecoder {
 public:
    explicit CalqDecoder(const Options &options);
    ~CalqDecoder();

    void decode();

 private:
    std::unique_ptr<CQFile> cqFile_;
    std::unique_ptr<File> qualFile_;
    std::unique_ptr<SAMFile> sideInformationFile_;
};

}  // namespace calq

#endif  // CALQ_CALQDECODER_H_

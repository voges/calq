/** @file CalqDecoder.h
 *  @brief This file contains the definition of the CalqDecoder class.
 */

// Copyright 2015-2017 Leibniz Universitaet Hannover

#ifndef CALQ_CALQDECODER_H_
#define CALQ_CALQDECODER_H_

#include "Common/Options.h"

namespace calq {

class CQFile;
class File;
class SAMFile;

class CalqDecoder {
 public:
    explicit CalqDecoder(const Options &pptions);
    ~CalqDecoder(void);

    void decode(void);

 private:
    CQFile* cqFile_;
    File* qualFile_;
    SAMFile* sideInformationFile_;
};

}  // namespace calq

#endif  // CALQ_CALQDECODER_H_


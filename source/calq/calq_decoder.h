#ifndef CALQ_CALQ_DECODER_H_
#define CALQ_CALQ_DECODER_H_

#include <memory>

#include "calq/options.h"

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

#endif  // CALQ_CALQ_DECODER_H_

#ifndef CALQ_CALQ_ENCODER_H_
#define CALQ_CALQ_ENCODER_H_

#include <memory>

#include "calq/options.h"

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

#endif  // CALQ_CALQ_ENCODER_H_

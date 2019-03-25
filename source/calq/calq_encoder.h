#ifndef CALQ_CALQ_ENCODER_H_
#define CALQ_CALQ_ENCODER_H_

#include <string>
#include <vector>

#include "calq/options.h"
#include "calq/cq_file.h"
#include "calq/sam_file.h"

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

#endif  // CALQ_CALQ_ENCODER_H_

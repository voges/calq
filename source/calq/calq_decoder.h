#ifndef CALQ_CALQ_DECODER_H_
#define CALQ_CALQ_DECODER_H_

#include "calq/options.h"
#include "calq/cq_file.h"
#include "calq/file.h"
#include "calq/sam_file.h"

namespace calq {

class CalqDecoder {
 public:
    explicit CalqDecoder(const Options &pptions);
    ~CalqDecoder(void);

    void decode(void);

 private:
    CQFile cqFile_;
    File qualFile_;
    SAMFile sideInformationFile_;
};

}  // namespace calq

#endif  // CALQ_CALQ_DECODER_H_

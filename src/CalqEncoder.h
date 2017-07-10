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
#if MPEG_CE5_DESCRIPTOR_STREAMS_OUTPUT
    File mpegParametersSetFile_;
    File mpegQVCI0File_;
    File mpegQVI0File_;
#endif
#if MPEG_CE5_DESCRIPTOR_STREAMS_COMPRESSION_EXTENSION
    CQFile qvciFile_;
    CQFile qviFile_;
#endif
    std::string inputFileName_;
    int polyploidy_;
    int qualityValueMin_;
    int qualityValueMax_;
    int qualityValueOffset_;
    std::vector<std::string> referenceFileNames_;
    SAMFile samFile_;
};

} // namespace calq

#endif // CALQ_CALQENCODER_H_


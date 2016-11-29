/** @file CalqEncoder.h
 *  @brief This file contains the definition of the CalqEncoder class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef CQ_CALQENCODER_H
#define CQ_CALQENCODER_H

#include "Common/CLIOptions.h"
#include "IO/CQ/CQFile.h"
#include "IO/SAM/SAMFile.h"

namespace cq {

class CalqEncoder {
public:
    explicit CalqEncoder(const CLIOptions &cliOptions);
    ~CalqEncoder(void);

    void encode(void);

private:
    size_t blockSize_;
    CQFile cqFile_;
    unsigned int polyploidy_;
    std::vector<std::string> referenceFileNames_;
    SAMFile samFile_;
};

}

#endif // CQ_CALQENCODER_H


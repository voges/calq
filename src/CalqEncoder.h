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
    size_t m_blockSize;
    CQFile m_cqFile;
    unsigned int m_polyploidy;
    std::vector<std::string> m_referenceFileNames;
    SAMFile m_samFile;
};

}

#endif // CQ_CALQENCODER_H


/** @file CalqCodec.h
 *  @brief This file contains the definitions of the CalqEncoder and
 *         CalqDecoder classes.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef CQ_CALQCODEC_H
#define CQ_CALQCODEC_H

#include "Common/CLIOptions.h"
#include "IO/CQ/CQFile.h"
#include "IO/File.h"
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

class CalqDecoder {
public:
    explicit CalqDecoder(const CLIOptions &cliOptions);
    ~CalqDecoder(void);

    void decode(void);

private:
    CQFile m_cqFile;
    File m_qualFile;
    SAMFile m_sideInformationFile;
};

}

#endif // CQ_CALQCODEC_H


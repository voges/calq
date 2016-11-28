/** @file CalqDecoder.h
 *  @brief This file contains the definition of the CalqDecoder class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef CQ_CALQDECODER_H
#define CQ_CALQDECODER_H

#include "Common/CLIOptions.h"
#include "IO/CQ/CQFile.h"
#include "IO/File.h"
#include "IO/SAM/SAMFile.h"

namespace cq {

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

#endif // CQ_CALQDECODER_H


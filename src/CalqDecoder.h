/** @file CalqDecoder.h
 *  @brief This file contains the definition of the CalqDecoder class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#ifndef CALQ_CALQDECODER_H_
#define CALQ_CALQDECODER_H_

#include "Common/Options.h"
#include "IO/CQ/CQFile.h"
#include "IO/File.h"
#include "IO/SAM/SAMFile.h"

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

} // namespace calq

#endif // CALQ_CALQDECODER_H_


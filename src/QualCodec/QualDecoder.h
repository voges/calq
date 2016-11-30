/** @file QualDecoder.h
 *  @brief This file contains the definition of the QualDecoder class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#ifndef CALQ_QUALCODEC_QUALDECODER_H_
#define CALQ_QUALCODEC_QUALDECODER_H_

#include <string>

#include "IO/CQ/CQFile.h"

namespace calq {

class QualDecoder {
public:
    explicit QualDecoder(void);
    ~QualDecoder(void);

    size_t readBlock(CQFile &cqFile);

private:
    void reset(void);

    std::string unmappedQualityValues_;
    std::string mappedQuantizerIndices_;
    std::string mappedQualityValueIndices_;
};

}

#endif // CALQ_QUALCODEC_QUALDECODER_H_


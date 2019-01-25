#ifndef CALQAPP_PROGRAM_OPTIONS_H_
#define CALQAPP_PROGRAM_OPTIONS_H_

// -----------------------------------------------------------------------------

#include <string>
#include <vector>

// -----------------------------------------------------------------------------

#include "calq/calq_coder.h"

// -----------------------------------------------------------------------------

namespace calqapp {

// -----------------------------------------------------------------------------

class ProgramOptions
{
 public:
    ProgramOptions(
            int argc,
            char *argv[]
    );

    ~ProgramOptions();

    void validate();

    calq::EncodingOptions options;

    bool force;
    bool debugStreams;
    bool test;
    std::string inputFilePath;
    std::string outputFilePath;
    size_t blockSize;
    std::string qualityValueType;
    std::string referenceFilePath;
    std::string filterTypeStr;
    std::string quantizerTypeStr;
    std::string versionStr;
    bool decompress;
    std::string sideInformationFilePath;


    size_t quantizationMin;
    size_t quantizationMax;
    size_t polyploidy;
    size_t hqSoftClipThreshold;

 private:
    void processCommandLine(
            int argc,
            char *argv[]
    );
};

// -----------------------------------------------------------------------------

}  // namespace calqapp

// -----------------------------------------------------------------------------

#endif  // CALQAPP_PROGRAM_OPTIONS_H_

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
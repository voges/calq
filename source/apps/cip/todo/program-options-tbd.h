#ifndef CIP_PROGRAM_OPTIONS_H_
#define CIP_PROGRAM_OPTIONS_H_

#include <string>

namespace cip {

class ProgramOptions {
   public:
    ProgramOptions(int argc, char *argv[]);

    void validate();

    bool force;
    bool debugStreams;
    bool test;
    bool help;
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

    size_t minNumQuantSteps;
    size_t manNumQuantSteps;
    size_t polyploidy;
    size_t hqSoftClipThreshold;

   private:
    void processCommandLine(int argc, char *argv[]);

    void validateV1();
    void validateV2();
    static void validateCompress();
    void validateDecompress();
    void validateCommon();
};

}  // namespace cip

#endif  // CIP_PROGRAM_OPTIONS_H_

#ifndef CALQAPP_PROGRAM_OPTIONS_H_
#define CALQAPP_PROGRAM_OPTIONS_H_


#include <string>
#include <vector>


namespace calqapp {


class ProgramOptions
{
 public:
    ProgramOptions(
            int argc,
            char *argv[]
    );

    ~ProgramOptions();

    void validate();
    
    enum struct QuantizerType {
        NONE,
        UNIFORM,
        LLOYD_MAX
    };

    enum struct FilterType {
        NONE,
        GAUSS,
        RECTANGLE
    };

    enum struct Version {
        NONE,
        V1,
        V2
    };

    bool force;
    bool verbose;
    bool test;
    bool squash;
    std::string inputFilePath;
    std::string outputFilePath;
    size_t blockSize;
    size_t filterSize;
    size_t quantizationMin;
    size_t quantizationMax;
    size_t polyploidy;
    size_t qualityValueMax;
    size_t qualityValueMin;
    size_t qualityValueOffset;
    std::string qualityValueType;
    std::string referenceFilePath;
    FilterType filterType;
    std::string filterTypeStr;
    QuantizerType quantizerType;
    std::string quantizerTypeStr;
    Version version;
    std::string versionStr;
    bool decompress;
    std::string sideInformationFilePath;

 private:
    void processCommandLine(
            int argc,
            char *argv[]
    );
};


}  // namespace calqapp


#endif  // CALQAPP_PROGRAM_OPTIONS_H_

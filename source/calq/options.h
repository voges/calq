#ifndef CALQ_OPTIONS_H_
#define CALQ_OPTIONS_H_

#include <string>
#include <vector>

namespace calq {

struct Options {
    Options();
    ~Options();

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

    // Options for both compression and decompression
    bool force;
    bool debug;
    bool test;
    bool squash;
    std::string inputFileName;
    std::string outputFileName;
    // Options for only compression
    size_t blockSize;
    size_t filterSize;
    size_t quantizationMin;
    size_t quantizationMax;
    size_t polyploidy;
    size_t qualityValueMax;
    size_t qualityValueMin;
    size_t qualityValueOffset;
    std::string qualityValueType;
    FilterType filterType;
    std::string filterTypeStr;
    QuantizerType quantizerType;
    std::string quantizerTypeStr;
    std::string referenceFileNames;
    // Options for only decompression
    bool decompress;
    std::string sideInformationFileName;
    Version version;
    std::string versionStr;
};

}  // namespace calq

#endif  // CALQ_OPTIONS_H_

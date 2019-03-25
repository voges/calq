#ifndef CALQ_OPTIONS_H_
#define CALQ_OPTIONS_H_

#include <string>
#include <vector>

namespace calq {

struct Options {
    Options(void);
    ~Options(void);

    void validate(void);

    // Options for both compression and decompression
    bool force;
    std::string inputFileName;
    std::string outputFileName;
    // Options for only compression
    int blockSize;
    int polyploidy;
    int qualityValueMax;
    int qualityValueMin;
    int qualityValueOffset;
    std::string qualityValueType;
    std::vector<std::string> referenceFileNames;
    // Options for only decompression
    bool decompress;
    std::string sideInformationFileName;
};

}  // namespace calq

#endif  // CALQ_OPTIONS_H_

/** @file Options.h
 *  @brief This file contains the definition of the Options struct.
 */

// Copyright 2015-2017 Leibniz Universitaet Hannover

#ifndef CALQ_COMMON_OPTIONS_H_
#define CALQ_COMMON_OPTIONS_H_

#include <string>
#include <vector>

namespace calq {

struct Options {
    Options(void);
    ~Options(void);

    void validate(void);

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

    // Options for both compression and decompression
    bool force;
    bool debug;
    bool test;
    bool squash;
    std::string inputFileName;
    std::string outputFileName;
    // Options for only compression
    int blockSize;
    int filterSize;
    int quantizationMin;
    int quantizationMax;
    int polyploidy;
    int qualityValueMax;
    int qualityValueMin;
    int qualityValueOffset;
    std::string qualityValueType;
    FilterType filterType;
    std::string filterTypeStr;
    QuantizerType quantizerType;
    std::string quantizerTypeStr;
    std::string referenceFileNames;
    // Options for only decompression
    bool decompress;
    std::string sideInformationFileName;
};

}  // namespace calq

#endif  // CALQ_COMMON_OPTIONS_H_


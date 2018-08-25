/** @file Options.cc
 *  @brief This file contains the implementation of the Options struct.
 */

// Copyright 2015-2017 Leibniz Universitaet Hannover

#include "Common/Options.h"
#include "Common/log.h"
#include "Common/ErrorExceptionReporter.h"

namespace calq {

Options::Options()
// Options for both compression and decompression
        : force(false),
          debug(false),
          test(false),
          inputFileName(""),
          outputFileName(""),
        // Options for only compression
          blockSize(0),
          filterSize(0),
          quantizationMin(0),
          quantizationMax(0),
          polyploidy(0),
          qualityValueMax(0),
          qualityValueMin(0),
          qualityValueOffset(0),
          qualityValueType(""),
          filterType(FilterType::NONE),
          quantizerType(QuantizerType::NONE),
          referenceFileNames(),
        // Options for only decompression
          decompress(false),
          sideInformationFileName(""),
          version(Version::NONE),
          versionStr("") {
}

Options::~Options() = default;

void Options::validate() {
    if (versionStr == "v1") {
        version = Version::V1;
        CALQ_LOG("Using CALQ version 1");
    } else if (versionStr == "v2") {
        version = Version::V2;
        CALQ_LOG("Using CALQ version 2");
    } else {
        throwErrorException("Unknown CALQ version");
    }

    // force
    if (force) {
        CALQ_LOG("Force switch set - overwriting output file(s)");
    }

    if (debug) {
        CALQ_LOG("Debug switch set - verbose output");
    }

    if (test) {
        CALQ_LOG("Test switch set - running test cases instead of compression");
    }

    if (squash) {
        CALQ_LOG("Acitivity scores are squashed between 0.0 and 1.0");
    } else {
        CALQ_LOG("Acitivity scores are !NOT! squashed between 0.0 and 1.0");
    }

    // inputFileName
    CALQ_LOG("Input file name: %s", inputFileName.c_str());
    if (inputFileName.empty()) {
        throwErrorException("No input file name provided");
    }
    if (!decompress) {
        if (fileNameExtension(inputFileName) != std::string("sam")) {
            throwErrorException("Input file name extension must be 'sam'");
        }
    } else {
        if (fileNameExtension(inputFileName) != std::string("cq")) {
            CALQ_LOG("Warning: Input file name extension is not 'cq'");
//             throwErrorException("Input file name extension must be 'cq'");
        }
    }
    if (!fileExists(inputFileName)) {
        throwErrorException("Cannot access input file");
    }

    // outputFileName
    if (!decompress) {
        if (outputFileName.empty()) {
            CALQ_LOG("No output file name provided - constructing output file name from input file name");
            outputFileName += inputFileName + ".cq";
        }
    } else {
        if (outputFileName.empty()) {
            CALQ_LOG("No output file name provided - constructing output file name from input file name");
            outputFileName += inputFileName + ".qual";
        }
    }
    CALQ_LOG("Output file name: %s", outputFileName.c_str());
    if (fileExists(outputFileName)) {
        if (!force) {
            throwErrorException("Not overwriting output file (use option 'f' to force overwriting)");
        }
    }

    // blockSize
    if (!decompress) {
        CALQ_LOG("Block size: %d", static_cast<int>(blockSize));
        if (blockSize < 1) {
            throwErrorException("Block size must be greater than 0");
        }
    }

    // Haplotyper filter size
    if (!decompress) {
        CALQ_LOG("Filter size: %d", static_cast<int>(filterSize));
        if (filterSize < 1) {
            throwErrorException("Filter size must be greater than 0");
        }
    }

    // Quantization
    if (!decompress) {
        CALQ_LOG("Quantization min steps: %d", static_cast<int>(quantizationMin));
        if (quantizationMin < 2) {
            throwErrorException("Quantization must be greater than 1");
        }
    }

    // Quantization
    if (!decompress) {
        CALQ_LOG("Quantization max steps: %d", static_cast<int>(quantizationMax));
        if (quantizationMax < 2 || quantizationMax < quantizationMin) {
            throwErrorException("Quantization must be greater than 1 and quantizationMin");
        }
    }

    // polyploidy
    if (!decompress) {
        CALQ_LOG("Polyploidy: %d", static_cast<int>(polyploidy));
        if (polyploidy < 1) {
            throwErrorException("Polyploidy must be greater than 0");
        }
    }

    // qualityValueType
    if (!decompress) {
        CALQ_LOG("Filter type: %s", filterTypeStr.c_str());
        if (filterTypeStr == "Gauss") {
            filterType = FilterType::GAUSS;
        } else if (filterTypeStr == "Rectangle") {
            filterType = FilterType::RECTANGLE;
        } else {
            throwErrorException("Filter type not supported");
        }
    }

    // qualityValueType
    if (!decompress) {
        CALQ_LOG("Quantizer type: %s", quantizerTypeStr.c_str());
        if (quantizerTypeStr == "Uniform") {
            quantizerType = QuantizerType::UNIFORM;
        } else if (quantizerTypeStr == "Lloyd") {
            quantizerType = QuantizerType::LLOYD_MAX;
        } else {
            throwErrorException("Quantizer type not supported");
        }
    }

    // qualityValueType
    if (!decompress) {
        CALQ_LOG("Quality value type: %s", qualityValueType.c_str());
        if (qualityValueType == "Sanger") {
            // Sanger: Phred+33 [0,40]
            qualityValueOffset = 33;
            qualityValueMin = 0;
            qualityValueMax = 40;
        } else if (qualityValueType == "Illumina-1.3+") {
            // Illumina 1.3+: Phred+64 [0,40]
            qualityValueOffset = 64;
            qualityValueMin = 0;
            qualityValueMax = 40;
        } else if (qualityValueType == "Illumina-1.5+") {
            // Illumina 1.5+: Phred+64 [0,40] with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator ('B')
            CALQ_LOG("Warning: Read Segment Quality Control Indicator will not be treated specifically by CALQ");
            qualityValueOffset = 64;
            qualityValueMin = 0;
            qualityValueMax = 40;
        } else if (qualityValueType == "Illumina-1.8+") {
            // Illumina 1.8+ Phred+33 [0,41]
            qualityValueOffset = 33;
            qualityValueMin = 0;
            qualityValueMax = 41;
        } else if (qualityValueType == "Max33") {
            // Max33 Phred+33 [0,93]
            qualityValueOffset = 33;
            qualityValueMin = 0;
            qualityValueMax = 93;
        } else if (qualityValueType == "Max64") {
            // Max64 Phred+64 [0,62]
            qualityValueOffset = 64;
            qualityValueMin = 0;
            qualityValueMax = 62;
        } else {
            throwErrorException("Quality value type not supported");
        }
        CALQ_LOG("Quality value offset: %d", static_cast<int>(qualityValueOffset));
        CALQ_LOG("Quality value range: [%d,%d]", static_cast<int>(qualityValueMin), static_cast<int>(qualityValueMax));

        // referenceFiles
        if (!decompress) {
            if (referenceFileNames.empty()) {
                CALQ_LOG("Operating without reference file");
            } else {
                CALQ_LOG("Operating with reference file:");
                CALQ_LOG("  %s", referenceFileNames.c_str());
                if (referenceFileNames.empty()) {
                    throwErrorException("Reference file name not provided");
                }
                if (!fileExists(referenceFileNames)) {
                    throwErrorException("Cannot access reference file");
                }
                if (fileNameExtension(referenceFileNames) != std::string("fa")
                    && fileNameExtension(referenceFileNames) != std::string("fasta")) {
                    throwErrorException("Reference file name extension must be 'fa' or 'fasta'");
                }
            }
        }
    }

    // decompress
    if (!decompress) {
        CALQ_LOG("Compressing");
    } else {
        CALQ_LOG("Decompressing");
    }

    // sideInformationFileName
    if (decompress) {
        CALQ_LOG("Side information file name: %s", sideInformationFileName.c_str());
        if (sideInformationFileName.empty()) {
            throwErrorException("No side information file name provided");
        }
        if (fileNameExtension(sideInformationFileName) != std::string("sam")) {
            throwErrorException("Side information file name extension must be 'sam'");
        }
        if (!fileExists(sideInformationFileName)) {
            throwErrorException("Cannot access side information file");
        }
    }
}

}  // namespace calq


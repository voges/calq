/** @file CLIOptions.h
 *  @brief This file contains the CLIOptions class which is used to store
 *         command line interface options provided by the user or computed
 *         from user input.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef CLIOPTIONS_H
#define CLIOPTIONS_H

#include <string>
#include <vector>

class CLIOptions {
public:
    CLIOptions(void)
        // Options for both compression and decompression
        : force(false)
        , inputFile("")
        , outputFile("")
        // Options for only compression
        , blockSize(0)
        , polyploidy(0)
        , qualityValueMax(0)
        , qualityValueMin(0)
        , qualityValueType("")
        , referenceFiles()
        // Options for only decompression
        , decompress(false)
        , sideInformationFile("")

    {}
    ~CLIOptions(void) {}

public:
    // Options for both compression and decompression
    bool force;
    std::string inputFile;
    std::string outputFile;
    // Options for only compression
    int blockSize;
    int polyploidy;
    int qualityValueMax;
    int qualityValueMin;
    std::string qualityValueType;
    std::vector<std::string> referenceFiles;
    // Options for only decompression
    bool decompress;
    std::string sideInformationFile;
};

#endif // CLIOPTIONS_H


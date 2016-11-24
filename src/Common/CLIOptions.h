/** @file CLIOptions.h
 *  @brief This file contains the definition of the CLIOptions class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef CQ_CLIOPTIONS_H
#define CQ_CLIOPTIONS_H

#include <string>
#include <vector>

namespace cq {

class CLIOptions {
public:
    CLIOptions(void);
    ~CLIOptions(void);

public:
    // Options for both compression and decompression
    bool force;
    std::string inputFileName;
    std::string outputFileName;
    // Options for only compression
    int blockSize;
    int polyploidy;
    std::string qualityValueType;
    std::vector<std::string> referenceFileNames;
    // Options for only decompression
    bool decompress;
    std::string sideInformationFileName;
};

}

#endif // CQ_CLIOPTIONS_H


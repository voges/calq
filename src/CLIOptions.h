/** @file CLIOptions.h
 *  @brief This file contains the CLIOptions class. Usually, an instance of
 *         class is globally available to all linked modules.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#ifndef CLIOPTIONS_H
#define CLIOPTIONS_H

#include <string>
#include <vector>

class CLIOptions {
public:
    CLIOptions(void)
        : decompress(false)
        , force(false)
        , inFileName("")
        , outFileName("")
        , referenceFileNames()
    {}
    ~CLIOptions(void) {}

public:
    bool decompress;
    bool force;
    std::string inFileName;
    std::string outFileName;
    std::vector<std::string> referenceFileNames;
};

extern CLIOptions cliOptions;

#endif // CLIOPTIONS_H


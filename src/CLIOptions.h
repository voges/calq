/** @file CLIOptions.h
 *  @brief This file contains the CLIOptions class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#ifndef CLIOPTIONS_H
#define CLIOPTIONS_H

#include <string>
#include <vector>

/** @brief Class: CLIOptions
 *
 *  The CLIOptions class is a container for command line options provided by
 *  the user.
 */
class CLIOptions {
public:
    CLIOptions(void)
        : decompress(false)
        , force(false)
        , inFileName("")
        , outFileName("")
        , refFileNames()
    {}
    ~CLIOptions(void) {}

public:
    bool decompress;
    bool force;
    std::string inFileName;
    std::string outFileName;
    std::vector<std::string> refFileNames;
};

extern CLIOptions cliOptions; // this instance is globally accessible

#endif // CLIOPTIONS_H


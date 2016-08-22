/** @file CLIOptions.h
 *  @brief This file contains the CLIOptions class.
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
        : blockSize(0)
        , decompress(false)
        , force(false)
        , inFileName("")
        , outFileName("")
        , polyploidy(0)
        , quantizedPrintout(false)
        , refFileNames()
        , samFileName("")
        , type("")
        , verbose(false)
    {}
    ~CLIOptions(void) {}

public:
    int blockSize;
    bool decompress;
    bool force;
    std::string inFileName;
    std::string outFileName;
    int polyploidy;
    bool quantizedPrintout;
    std::vector<std::string> refFileNames;
    std::string samFileName;
    std::string type;
    bool verbose;
};

extern CLIOptions cliOptions; // this instance is globally accessible

#endif // CLIOPTIONS_H


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
        : alignmentsFileName("")
        , blockSize(0)
        , decompress(false)
        , force(false)
        , inFileName("")
        , outFileName("")
        , polyploidy(0)
        , quantizedPrintout(false)
        , qvMin(0)
        , qvMax(0)
        , refFileNames()
        , stats(false)
        , type("")
        , verbose(false)
    {}
    ~CLIOptions(void) {}

public:
    std::string alignmentsFileName;
    int blockSize;
    bool decompress;
    bool force;
    std::string inFileName;
    std::string outFileName;
    int polyploidy;
    bool quantizedPrintout;
    int qvMin;
    int qvMax;
    std::vector<std::string> refFileNames;
    bool stats;
    std::string type;
    bool verbose;
};

//extern CLIOptions cliOptions; // this instance is globally accessible

#endif // CLIOPTIONS_H


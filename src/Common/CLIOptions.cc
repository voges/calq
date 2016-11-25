/** @file CLIOptions.cc
 *  @brief This file contains the implementation of the CLIOptions class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "Common/CLIOptions.h"

cq::CLIOptions::CLIOptions(void)
    // Options for both compression and decompression
    : force(false)
    , inputFileName("")
    , outputFileName("")
    // Options for only compression
    , blockSize(0)
    , polyploidy(0)
    , referenceFileNames()
    // Options for only decompression
    , decompress(false)
    , sideInformationFileName("")
{
    // empty
}

cq::CLIOptions::~CLIOptions(void)
{
    // empty
}


/** @file calq.cc
 *  @brief This file contains the main function of the calq compression tool.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  2016-11-21: Restructured and simplified code (voges)
 *  2016-05-30: Added reference file(s) CLI argument (voges)
 *  2016-05-29: Added TCLAP for command line argument parsing; this btw
 *              is then compatible with Windows, too (as there is no getopt.h
 *              on Windows) (voges)
 */

#include "CalqCodec.h"
#include "cmake_config.h"
#include "Common/CLIOptions.h"
#include "Common/Exceptions.h"
#include "Common/helpers.h"
#include "Common/os_config.h"
#include "tclap/CmdLine.h"

#ifdef CQ_OS_WINDOWS
    #define TCLAP_NAMESTARTSTRING "~~"
    #define TCLAP_FLAGSTARTSTRING "/"
#else
    #define TCLAP_NAMESTARTSTRING "--"
    #define TCLAP_FLAGSTARTSTRING "-"
#endif

static void printVersionAndCopyright(void)
{
    printf("-----------------------------------------------\n");
    printf("Program: CALQ\n");
    printf("Version: %s\n", CQ_VERSION);
    printf("Build time: %s\n", CQ_TIMESTAMP_UTC);
    printf("Git commit hash: %s\n", CQ_GIT_COMMIT_HASH_SHORT);
    printf("-----------------------------------------------\n");
    printf("Copyright (c) 2015-%d\n", CQ_BUILD_YEAR);
    printf("Leibniz Universitaet Hannover\n");
    printf("Institut fuer Informationsverarbeitung (TNT)\n");
    printf("Contact: Jan Voges <voges@tnt.uni-hannover.de>\n");
    printf("-----------------------------------------------\n");
    printf("Please report bugs to voges@tnt.uni-hannover.de\n");
    printf("-----------------------------------------------\n");
}

static void checkAndProcessCLIOptions(cq::CLIOptions &cliOptions)
{
    // force
    if (cliOptions.force == true) {
        LOG("Force switch set - overwriting output file(s)");
    }

    // inputFile
    LOG("Input file name: %s", cliOptions.inputFileName.c_str());
    if (cliOptions.inputFileName.empty() == true) {
        throwErrorException("No input file name provided");
    }
    if (cliOptions.decompress == false) {
        if (cq::fileNameExtension(cliOptions.inputFileName) != std::string("sam")) {
            throwErrorException("Input file name extension must be 'sam'");
        }
    } else {
        if (cq::fileNameExtension(cliOptions.inputFileName) != std::string("cq")) {
            throwErrorException("Input file name extension must be 'cq'");
        }
    }
    if (cq::fileExists(cliOptions.inputFileName) == false) {
        throwErrorException("Cannot access input file");
    }

    // outputFile
    if (cliOptions.decompress == false) {
        if (cliOptions.outputFileName.empty() == true) {
            LOG("No output file name provided - constructing output file name from input file name");
            cliOptions.outputFileName += cliOptions.inputFileName + ".cq";
            LOG("Output file name: %s", cliOptions.outputFileName.c_str());
        }
    } else {
        if (cliOptions.outputFileName.empty() == true) {
            LOG("No output file name provided - constructing output file name from input file name");
            cliOptions.outputFileName += cliOptions.inputFileName + ".qual";
            LOG("Output file name: %s", cliOptions.outputFileName.c_str());
        }
    }
    if (cq::fileExists(cliOptions.outputFileName) == true) {
        if (cliOptions.force == false) {
            throwErrorException("Not overwriting output file (use option 'f' to force overwriting)");
        }
    }

    // blockSize
    if (cliOptions.decompress == false) {
        LOG("Block size: %d", cliOptions.blockSize);
        if (cliOptions.blockSize < 1) {
            throwErrorException("Block size must be greater than 0");
        }
    }

    // polyploidy
    if (cliOptions.decompress == false) {
        LOG("Polyploidy: %d", cliOptions.polyploidy);
        if (cliOptions.polyploidy < 1) {
            throwErrorException("Polyploidy must be greater than 0");
        }
    }

    // referenceFiles
    if (cliOptions.decompress == false) {
        if (cliOptions.referenceFileNames.empty() == true) {
            LOG("Operating without reference file(s)");
        } else {
            LOG("Operating with %zu reference file(s):", cliOptions.referenceFileNames.size());
            for (std::string const &referenceFileName : cliOptions.referenceFileNames) {
                LOG("  %s", referenceFileName.c_str());
                if (referenceFileName.empty() == true) {
                    throwErrorException("Reference file name not proviced");
                }
                if (cq::fileExists(referenceFileName) == false) {
                    throwErrorException("Cannot access reference file");
                }
                if (cq::fileNameExtension(referenceFileName) != std::string("fa")
                    && cq::fileNameExtension(referenceFileName) != std::string("fasta")) {
                    throwErrorException("Reference file name extension must be 'fa' or 'fasta'");
                }
            }
        }
    }

    // decompress
    if (cliOptions.decompress == false) {
        LOG("Compressing");
    } else {
        LOG("Decompressing");
    }

    // sideInformationFile
    if (cliOptions.decompress == true) {
        LOG("Side information file name: %s", cliOptions.sideInformationFileName.c_str());
        if (cliOptions.sideInformationFileName.empty() == true) {
            throwErrorException("No side information file name provided");
        }
        if (cq::fileNameExtension(cliOptions.sideInformationFileName) != std::string("si")) {
            throwErrorException("Side information file name extension must be 'si'");
        }
        if (cq::fileExists(cliOptions.sideInformationFileName) == false) {
            throwErrorException("Cannot access side information file");
        }
    }
}

int main(int argc, char *argv[])
{
    printVersionAndCopyright();

    try {
        cq::CLIOptions cliOptions;

        // TCLAP class
        TCLAP::CmdLine cmd("CALQ", ' ', CQ_VERSION);

        // TCLAP arguments (both compression and decompression)
        TCLAP::SwitchArg forceSwitch("f", "force", "Force overwriting of output files etc.", cmd, false);
        TCLAP::UnlabeledValueArg<std::string> inputFileNameArg("inputFileName", "Input file name", true, "", "string", cmd);
        TCLAP::ValueArg<std::string> outputFileNameArg("o", "outputFileName", "Output file name", false, "", "string", cmd);

        // TCLAP arguments (only compression)
        TCLAP::ValueArg<int> blockSizeArg("b", "blockSize", "Block size (in number of SAM records)", false, 10000, "int", cmd);
        TCLAP::ValueArg<int> polyploidyArg("p", "polyploidy", "Polyploidy", false, 2, "int", cmd);
        TCLAP::MultiArg<std::string> referenceFileNamesArg("r", "referenceFileNames", "Reference file name(s) (FASTA format)", false, "string", cmd);

        // TCLAP arguments (only decompression)
        TCLAP::SwitchArg decompressSwitch("d", "decompress", "Decompress", cmd, false);
        TCLAP::ValueArg<std::string> sideInformationFileNameArg("s", "sideInformationFilName", "Side information file name", false, "", "string", cmd);

        // Let the TCLAP class parse the provided arguments
        cmd.parse(argc, argv);

        // Check for sanity in compression mode
        if (decompressSwitch.isSet() == false) {
            //if (referenceFileNamesArg.isSet() == false) {
            //    throwErrorException("Argument 'r' required in compression mode");
            //}
            if (sideInformationFileNameArg.isSet() == true) {
                throwErrorException("Argument 's' forbidden in compression mode");
            }
        }

        // Check for sanity in decompression mode
        if (decompressSwitch.isSet() == true) {
            if (blockSizeArg.isSet() == true) {
                throwErrorException("Argument 'b' forbidden in decompression mode");
            }
            if (polyploidyArg.isSet() == true) {
                throwErrorException("Argument 'p' forbidden in decompression mode");
            }
            if (referenceFileNamesArg.isSet() == true) {
                throwErrorException("Argument 'r' forbidden in decompression mode");
            }
            if (sideInformationFileNameArg.isSet() == false) {
                throwErrorException("Argument 's' required in decompression mode");
            }
        }

        // Get the value parsed by each arg
        cliOptions.force = forceSwitch.getValue();
        cliOptions.inputFileName = inputFileNameArg.getValue();
        cliOptions.outputFileName = outputFileNameArg.getValue();
        cliOptions.blockSize = blockSizeArg.getValue();
        cliOptions.polyploidy = polyploidyArg.getValue();
        cliOptions.referenceFileNames = referenceFileNamesArg.getValue();
        cliOptions.decompress = decompressSwitch.getValue();
        cliOptions.sideInformationFileName = sideInformationFileNameArg.getValue();

        // Check CLIOptions
        checkAndProcessCLIOptions(cliOptions);

        // Compress or decompress
        if (cliOptions.decompress == false) {
            cq::CalqEncoder calqEncoder(cliOptions);
            calqEncoder.encode();
        } else {
            cq::CalqDecoder calqDecoder(cliOptions);
            calqDecoder.decode();
        }

    } catch (TCLAP::ArgException &tclapException) {
        ERROR("%s (argument: %s)", tclapException.error().c_str(), tclapException.argId().c_str());
        return EXIT_FAILURE;
    } catch (const cq::ErrorException &errorException) {
        ERROR("%s", errorException.what());
        return EXIT_FAILURE;
    } catch (const std::exception &stdException) {
        ERROR("Fatal: %s", stdException.what());
        return EXIT_FAILURE;
    } catch (...) {
        ERROR("Unkown error occured");
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


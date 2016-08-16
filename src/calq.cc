/** @file calq.cc
 *  @brief This file contains the main function of the calq compression tool.
 *
 *  This file contains the main function (and some help functions) of the calq
 *  compression tool.
 *  Encoding: The program takes a SAM file as input and encodes the quality
 *            values. The encoded bitstream is written to a CQ file (*.cq).
 *  Decoding: The program reads in the compressed CQ file and reconstructs
 *            the quality scores. The quality scores are written to a QUAL
 *            file (*.qual). The SAM file is -not- reconstructed.
 *
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  2016-05-30: Added reference file(s) (voges)
 *  2016-05-29: Added TCLAP for command line argument parsing; this btw
 *              is then compatible with Windows, too (as there is no getopt.h
 *              on Windows) (voges)
 */

#include "Common/debug.h"
#include "Common/os_config.h"

#ifdef OS_WINDOWS
    #define TCLAP_NAMESTARTSTRING "~~"
    #define TCLAP_FLAGSTARTSTRING "/"
#else
    #define TCLAP_NAMESTARTSTRING "--"
    #define TCLAP_FLAGSTARTSTRING "-"
#endif

#include "CalqCodec.h"
#include "cmake_config.h"
#include "Common/CLIOptions.h"
#include "Common/Exceptions.h"
#include "Common/fileSystemHelpers.h"
#include "tclap/CmdLine.h"
#include <iostream>

CLIOptions cliOptions;

static void printVersion(void)
{
    std::cout << "Program: calq" << std::endl;
    std::cout << "Version: " << VERSION << std::endl;
    std::cout << "Build time: " << TIMESTAMP_UTC << std::endl;
    //std::cout << "Git branch: " << GIT_BRANCH << std::endl;
    std::cout << "Git commit hash: " << GIT_COMMIT_HASH_SHORT << std::endl;
    std::cout << std::endl;
}

static void printCopyright(void)
{
    std::cout << "Copyright (c) 2015-" << BUILD_YEAR << std::endl;
    std::cout << "Leibniz Universitaet Hannover" << std::endl;
    std::cout << "Institut fuer Informationsverarbeitung (TNT)" << std::endl;
    std::cout << "Contact: Jan Voges <voges@tnt.uni-hannover.de>" << std::endl;
    std::cout << std::endl;
}

int main(int argc, char *argv[])
{
    std::cout << std::setprecision(2) << std::fixed;

    printVersion();
    printCopyright();

    try {
        // TCLAP class
        TCLAP::CmdLine cmd("calq - Lossy compression of next-generation sequencing quality values", ' ', VERSION);

        // TCLAP arguments
        TCLAP::ValueArg<int> blockSizeArg("b", "blockSize", "Block size (in number of SAM records)", false, 10000, "int", cmd);
        TCLAP::SwitchArg decompressSwitch("d", "decompress", "Decompress CQ file", cmd, false);
        TCLAP::SwitchArg forceSwitch("f", "force", "Force overwriting of output files, etc.", cmd, false);
        TCLAP::UnlabeledValueArg<std::string> infileArg("infile", "Input file", true, "", "string", cmd);
        TCLAP::ValueArg<std::string> outfileArg("o", "outfile", "Output file", false, "", "string", cmd);
        TCLAP::ValueArg<int> polyploidyArg("p", "polyploidy", "Polyploidy", false, 2, "int", cmd);
        TCLAP::MultiArg<std::string> referenceArg("r", "reference", "Reference file(s) (FASTA format)", true, "string", cmd);
        TCLAP::ValueArg<std::string> typeArg("t", "type", "Type of quality values (sanger, illumina-1.3+, illumina-1.5+, illumina-1.8+)", false, "illumina-1.8+", "string", cmd);
        TCLAP::SwitchArg verboseSwitch("v", "verbose", "Verbose output", cmd, false);

        // Let the TCLAP class parse the provided arguments
        cmd.parse(argc, argv);

        // Check for sanity
        if (blockSizeArg.isSet() && decompressSwitch.isSet()) {
            throwErrorException("Combining arguments 'b' and 'd' is forbidden");
        }
        if (polyploidyArg.isSet() && decompressSwitch.isSet()) {
            throwErrorException("Combining arguments 'p' and 'd' is forbidden");
        }

        // Get the value parsed by each arg
        cliOptions.blockSize = blockSizeArg.getValue();
        cliOptions.decompress = decompressSwitch.getValue();
        cliOptions.force = forceSwitch.getValue();
        cliOptions.inFileName = infileArg.getValue();
        cliOptions.outFileName = outfileArg.getValue();
        cliOptions.polyploidy = polyploidyArg.getValue();
        cliOptions.refFileNames = referenceArg.getValue();
        cliOptions.type = typeArg.getValue();
        cliOptions.verbose = verboseSwitch.getValue();

        // Check block size
        if (cliOptions.blockSize < 1) {
            throwErrorException("Block size must be greater than 0");
        }
        if (decompressSwitch.isSet() == false) {
            std::cout << ME << "Using block size: " << cliOptions.blockSize << std::endl;
        }

        // Check polyploidy
        if (cliOptions.polyploidy < 1) {
            throwErrorException("Polyploidy must be greater than 0");
        }
        if (cliOptions.polyploidy > 6 && cliOptions.force == false) {
            throwErrorException("Polyploidy very high (use option 'f' to force processing)");
        }
        std::cout << ME << "Using polyploidy: " << cliOptions.polyploidy << std::endl;

        // Check type of quality values
        // Supported quality value ranges:
        //   Sanger         Phred+33   [0,40]
        //   Illumina 1.3+  Phred+64   [0,40]
        //   Illumina 1.5+  Phred+64   [0,40] with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator ('B')
        //   Illumina 1.8+  Phred+33   [0,41]
        int qvOffset = 0;
        int qvMin = 0;
        int qvMax = 0;
        if (cliOptions.type == "sanger") {
            qvOffset = 33;
            qvMin = qvOffset;
            qvMax = qvOffset + 40;
        } else if (cliOptions.type == "illumina-1.3+") {
            qvOffset = 64;
            qvMin = qvOffset;
            qvMax = qvOffset + 40;
        } else if (cliOptions.type == "illumina-1.5+") {
            qvOffset = 64;
            qvMin = qvOffset;
            qvMax = qvOffset + 40;
        } else if (cliOptions.type == "illumina-1.8+") {
            qvOffset = 33;
            qvMin = qvOffset;
            qvMax = qvOffset + 41;
        } else {
            throwErrorException("Quality value type not supported");
        }
        std::cout << ME << "Using quality value type: " << cliOptions.type << " (offset,min,max)=(" << qvOffset << "," << qvMin << "," << qvMax << ")" << std::endl;

        // Check verbosity
        if (cliOptions.verbose == true) {
            std::cout << ME << "Verbose output activated" << std::endl;
        }

        // Check if the input file exists and if it has the correct extension
        std::cout << ME << "Checking input file: " << cliOptions.inFileName << std::endl;
        if (!fileExists(cliOptions.inFileName)) {
            throwErrorException("Cannot access input file");
        }
        if (cliOptions.decompress == false) {
            if (fileNameExtension(cliOptions.inFileName) != std::string("sam")) {
                throwErrorException("Input file fileNameExtension must be 'sam'");
            }
        } else {
            if (fileNameExtension(cliOptions.inFileName) != std::string("cq")) {
                throwErrorException("Input file fileNameExtension must be 'cq'");
            }
        }
        std::cout << ME << "OK" << std::endl;

        // Create correct output file name if it was not provided via the
        // command line options and check if output file is already there and
        // if the user wants to overwrite it in this case
        std::cout << ME << "Checking output file: ";
        if (cliOptions.decompress == false) {
            if (cliOptions.outFileName.empty()) {
                cliOptions.outFileName.append(cliOptions.inFileName);
                cliOptions.outFileName.append(".cq");
            }
        } else {
            if (cliOptions.outFileName.empty()) {
                cliOptions.outFileName.append(cliOptions.inFileName);
                cliOptions.outFileName.append(".qual");
            }
        }
        std::cout << cliOptions.outFileName << std::endl;
        if (fileExists(cliOptions.outFileName) && cliOptions.force == false) {
            throwErrorException("Output file already exists (use option 'f' to force overwriting)");
        }
        std::cout << ME << "OK" << std::endl;

        // Check if the reference file(s) exist and check for FASTA file
        // ending(s)
        for (auto const &refFileName : cliOptions.refFileNames) {
            std::cout << ME << "Checking reference file: " << refFileName << std::endl;
            if (   fileNameExtension(refFileName) != std::string("fa")
                && fileNameExtension(refFileName) != std::string("fasta")) {
                throwErrorException("Reference file fileNameExtension must be 'fa' or 'fasta'");
            }
            if (!fileExists(refFileName)) {
                throwErrorException("Cannot access reference file");
            }
            std::cout << ME << "OK" << std::endl;
        }

        // Compress or decompress
        if (cliOptions.decompress == false) {
            CalqEncoder calqEncoder(cliOptions.inFileName, 
                                    cliOptions.outFileName, 
                                    cliOptions.refFileNames, 
                                    (unsigned int)cliOptions.blockSize, 
                                    (unsigned int)cliOptions.polyploidy, 
                                    qvOffset, 
                                    qvMin, 
                                    qvMax);
            calqEncoder.encode();
        } else {
            CalqDecoder calqDecoder(cliOptions.inFileName, 
                                    cliOptions.outFileName, 
                                    cliOptions.refFileNames);
            calqDecoder.decode();
        }
    std::cout << ME << "Finished" << std::endl;
    } catch (TCLAP::ArgException &tclapException) {
        std::cerr << ME << "Argument error: " << tclapException.error() << " (argument: " << tclapException.argId() << ")" << std::endl;
        return EXIT_FAILURE;
    } catch (const ErrorException &errorException) {
        std::cerr << ME << "Error: " << errorException.what() << std::endl;
        return EXIT_FAILURE;
    } catch (const std::exception &stdException) {
        std::cerr << ME << "Fatal error: " << stdException.what() << std::endl;
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << ME << "Unkown error occured" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


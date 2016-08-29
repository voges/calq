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

        // TCLAP arguments (both compression and decompression)
        TCLAP::SwitchArg forceSwitch("f", "force", "Force overwriting of output files, etc.", cmd, false);
        TCLAP::UnlabeledValueArg<std::string> infileArg("infile", "Input file", true, "", "string", cmd);
        TCLAP::ValueArg<std::string> outfileArg("o", "outfile", "Output file", false, "", "string", cmd);
        TCLAP::SwitchArg verboseSwitch("v", "verbose", "Verbose output", cmd, false);

        // TCLAP arguments (only compression)
        TCLAP::ValueArg<int> blockSizeArg("b", "blockSize", "Block size (in number of SAM records)", false, 10000, "int", cmd);
        TCLAP::ValueArg<int> polyploidyArg("p", "polyploidy", "Polyploidy", false, 2, "int", cmd);
        TCLAP::SwitchArg quantizedPrintoutSwitch("q", "quantizedPrintout", "Print quantized quality values", cmd, false);
        TCLAP::MultiArg<std::string> referenceArg("r", "reference", "Reference file(s) (FASTA format)", false, "string", cmd);
        TCLAP::ValueArg<std::string> typeArg("t", "type", "Type of quality values (sanger, illumina-1.3+, illumina-1.5+, illumina-1.8+)", false, "illumina-1.8+", "string", cmd);

        // TCLAP arguments (only decompression)
        TCLAP::SwitchArg decompressSwitch("d", "decompress", "Decompress CQ file", cmd, false);
        TCLAP::ValueArg<std::string> samfileArg("s", "samfile", "SAM file", false, "", "string", cmd);

        // Let the TCLAP class parse the provided arguments
        cmd.parse(argc, argv);

        // Check for sanity in compression mode
        if (decompressSwitch.isSet() == false) {
            if (referenceArg.isSet() == false) {
                throwErrorException("Argument 'r' required in compression mode");
            }
            if (samfileArg.isSet() == true) {
                throwErrorException("Argument 's' forbidden in compression mode");
            }
        }

        // Check for sanity in decompression mode
        if (decompressSwitch.isSet() == true) {
            if (samfileArg.isSet() == false) {
                throwErrorException("Argument 's' required in compression mode");
            }
            if (blockSizeArg.isSet()) {
                throwErrorException("Argument 'b' forbidden in decompression mode");
            }
            if (polyploidyArg.isSet() == true) {
                throwErrorException("Argument 'p' forbidden in decompression mode");
            }
            if (quantizedPrintoutSwitch.isSet() == true) {
                throwErrorException("Argument 'q' forbidden in decompression mode");
            }
            if (referenceArg.isSet() == true) {
                throwErrorException("Argument 'r' forbidden in decompression mode");
            }
            if (typeArg.isSet() == true) {
                throwErrorException("Argument 't' forbidden in decompression mode");
            }
        }

        // Get the value parsed by each arg
        cliOptions.blockSize = blockSizeArg.getValue();
        cliOptions.decompress = decompressSwitch.getValue();
        cliOptions.force = forceSwitch.getValue();
        cliOptions.inFileName = infileArg.getValue();
        cliOptions.outFileName = outfileArg.getValue();
        cliOptions.polyploidy = polyploidyArg.getValue();
        cliOptions.quantizedPrintout = quantizedPrintoutSwitch.getValue();
        cliOptions.refFileNames = referenceArg.getValue();
        cliOptions.samFileName = samfileArg.getValue();
        cliOptions.type = typeArg.getValue();
        cliOptions.verbose = verboseSwitch.getValue();

        std::cout << ME << "Checking command line arguments" << std::endl;

        // Check command line options for both compression and decompression
        if (cliOptions.force == true) {
            std::cout << ME << "Force switch set - overwriting output file(s)" << std::endl;
        }
        if (fileExists(cliOptions.inFileName) == false) {
            throwErrorException("Cannot access input file");
        }
        if (cliOptions.verbose == true) {
            std::cout << ME << "Verbose output activated" << std::endl;
        }

        // Check command line options for compression mode
        int qvMin = 0;
        int qvMax = 0;
        if (cliOptions.decompress == false) {
            if (fileNameExtension(cliOptions.inFileName) != std::string("sam")) {
                throwErrorException("Input file extension must be 'sam'");
            }
            std::cout << ME << "Input file: " << cliOptions.inFileName << std::endl;

            if (cliOptions.outFileName.empty()) {
                cliOptions.outFileName.append(cliOptions.inFileName);
                cliOptions.outFileName.append(".cq");
            }
            if ((fileExists(cliOptions.outFileName) == true) && (cliOptions.force == false)) {
                throwErrorException("Output file already exists (use option 'f' to force overwriting)");
            }
            std::cout << ME << "Output file: " << cliOptions.outFileName << std::endl;

            if (cliOptions.blockSize < 1) {
                throwErrorException("Block size must be greater than 0");
            }
            std::cout << ME << "Block size: " << cliOptions.blockSize << std::endl;

            if (cliOptions.polyploidy < 1) {
                throwErrorException("Polyploidy must be greater than 0");
            }
            if ((cliOptions.polyploidy > 6) && (cliOptions.force == false)) {
                throwErrorException("Polyploidy very high (use option 'f' to force processing)");
            }
            std::cout << ME << "Polyploidy: " << cliOptions.polyploidy << std::endl;

            if (cliOptions.quantizedPrintout == true) {
                std::cout << ME << "Printing out quantized quality values to stderr" << std::endl;
            }

            for (auto const &refFileName : cliOptions.refFileNames) {
                if (   fileNameExtension(refFileName) != std::string("fa")
                    && fileNameExtension(refFileName) != std::string("fasta")) {
                    throwErrorException("Reference file extension must be 'fa' or 'fasta'");
                }
                if (!fileExists(refFileName)) {
                    throwErrorException("Cannot access reference file");
                }
                std::cout << ME << "Reference file: " << refFileName << std::endl;
            }

            // Supported quality value ranges:
            //   Sanger         Phred+33   [0,40]
            //   Illumina 1.3+  Phred+64   [0,40]
            //   Illumina 1.5+  Phred+64   [0,40] with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator ('B')
            //   Illumina 1.8+  Phred+33   [0,41]
            if (cliOptions.type == "sanger") {
                qvMin = 33;
                qvMax = qvMin + 40;
            } else if (cliOptions.type == "illumina-1.3+") {
                qvMin = 64;
                qvMax = qvMin + 40;
            } else if (cliOptions.type == "illumina-1.5+") {
                qvMin = 64;
                qvMax = qvMin + 40;
            } else if (cliOptions.type == "illumina-1.8+") {
                qvMin = 33;
                qvMax = qvMin + 41;
            } else {
                throwErrorException("Quality value type not supported");
            }
            std::cout << ME << "Quality value type: " << cliOptions.type << " [" << qvMin << "," << qvMax << "]" << std::endl;
        }

        // Check command line options for compression mode
        if (cliOptions.decompress == true) {
            if (fileNameExtension(cliOptions.inFileName) != std::string("cq")) {
                throwErrorException("Input file extension must be 'cq'");
            }
            std::cout << ME << "Input file: " << cliOptions.inFileName << std::endl;

            if (cliOptions.outFileName.empty()) {
                cliOptions.outFileName.append(cliOptions.inFileName);
                cliOptions.outFileName.append(".qual");
            }
            if ((fileExists(cliOptions.outFileName) == true) && (cliOptions.force == false)) {
                throwErrorException("Output file already exists (use option 'f' to force overwriting)");
            }
            std::cout << ME << "Output file: " << cliOptions.outFileName << std::endl;

            // Check if the SAM file exists and if it has the correct extension
            if (fileExists(cliOptions.samFileName) == false) {
                throwErrorException("Cannot access SAM file");
            }
            if (fileNameExtension(cliOptions.samFileName) != std::string("sam")) {
                throwErrorException("SAM file extension must be 'sam'");
            }
            std::cout << ME << "SAM file: " << cliOptions.samFileName << std::endl;
        }

        // Compress or decompress
        if (cliOptions.decompress == false) {
            CalqEncoder calqEncoder(cliOptions.inFileName,
                                    cliOptions.outFileName,
                                    cliOptions.refFileNames,
                                    cliOptions.blockSize,
                                    cliOptions.polyploidy,
                                    cliOptions.quantizedPrintout,
                                    qvMin,
                                    qvMax);
            calqEncoder.encode();
        } else {
            CalqDecoder calqDecoder(cliOptions.inFileName,
                                    cliOptions.outFileName,
                                    cliOptions.samFileName);
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


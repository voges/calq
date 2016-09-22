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

#include "Common/helpers.h"
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
#include "Common/helpers.h"
#include "tclap/CmdLine.h"
#include <cctype>
#include <iostream>

CLIOptions cliOptions;

static void printVersionAndCopyright(void)
{
    LOG("----------------------------------------------");
    LOG("Program: calq");
    LOG("Version: %s", VERSION);
    LOG("Build time: %s", TIMESTAMP_UTC);
    LOG("Git commit hash: %s", GIT_COMMIT_HASH_SHORT);
    LOG("----------------------------------------------");
    LOG("Copyright (c) 2015-%d", BUILD_YEAR);
    LOG("Leibniz Universitaet Hannover");
    LOG("Institut fuer Informationsverarbeitung (TNT)");
    LOG("Contact: Jan Voges <voges@tnt.uni-hannover.de>");
    LOG("----------------------------------------------");
}

int main(int argc, char *argv[])
{
    printVersionAndCopyright();

    try {
        CLIOptions cliOptions;

        // TCLAP class
        TCLAP::CmdLine cmd("calq", ' ', VERSION);

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
        TCLAP::ValueArg<std::string> typeArg("t", "type", "Quality value type (sanger, illumina-1.3+, illumina-1.5+, illumina-1.8+, min:max)", false, "illumina-1.8+", "string", cmd);

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
            if (blockSizeArg.isSet() == true) {
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

        // Check command line options for both compression and decompression
        if (cliOptions.force == true) {
            LOG("Force switch set - overwriting output file(s)");
        }
        if (fileExists(cliOptions.inFileName) == false) {
            LOG("Input file name: %s", cliOptions.inFileName.c_str());
            throwErrorException("Cannot access input file");
        }
        if (cliOptions.verbose == true) {
            LOG("Verbose output activated");
        }

        // Check command line options for compression mode
        if (cliOptions.decompress == false) {
            LOG("Input file: %s", cliOptions.inFileName.c_str());
            if (fileNameExtension(cliOptions.inFileName) != std::string("sam")) {
                throwErrorException("Input file extension must be 'sam'");
            }

            if (cliOptions.outFileName.empty()) {
                cliOptions.outFileName.append(cliOptions.inFileName);
                cliOptions.outFileName.append(".cq");
            }
            LOG("Output file: %s", cliOptions.outFileName.c_str());
            if ((fileExists(cliOptions.outFileName) == true) && (cliOptions.force == false)) {
                throwErrorException("Output file already exists (use option 'f' to force overwriting)");
            }

            LOG("Block size: %d", cliOptions.blockSize);
            if (cliOptions.blockSize < 1) {
                throwErrorException("Block size must be greater than 0");
            }

            LOG("Polyploidy: %d", cliOptions.polyploidy);
            if (cliOptions.polyploidy < 1) {
                throwErrorException("Polyploidy must be greater than 0");
            }
            if ((cliOptions.polyploidy > 6) && (cliOptions.force == false)) {
                throwErrorException("Polyploidy very high (use option 'f' to force processing)");
            }

            if (cliOptions.quantizedPrintout == true) {
                LOG("Printing quantized quality values");
            }

            for (auto const &refFileName : cliOptions.refFileNames) {
                LOG("Reference file: %s", refFileName.c_str());
                if (   fileNameExtension(refFileName) != std::string("fa")
                    && fileNameExtension(refFileName) != std::string("fasta")) {
                    throwErrorException("Reference file extension must be 'fa' or 'fasta'");
                }
                if (!fileExists(refFileName)) {
                    throwErrorException("Cannot access reference file");
                }
            }

            // Supported quality value ranges:
            //   Sanger         Phred+33   [0,40]
            //   Illumina 1.3+  Phred+64   [0,40]
            //   Illumina 1.5+  Phred+64   [0,40] with 0=unused, 1=unused, 2=Read Segment Quality Control Indicator ('B')
            //   Illumina 1.8+  Phred+33   [0,41]
            LOG("Quality value type: %s", cliOptions.type.c_str());
            if (cliOptions.type == "sanger") {
                cliOptions.qvMin = 33;
                cliOptions.qvMax = cliOptions.qvMin + 40;
            } else if (cliOptions.type == "illumina-1.3+") {
                cliOptions.qvMin = 64;
                cliOptions.qvMax = cliOptions.qvMin + 40;
            } else if (cliOptions.type == "illumina-1.5+") {
                cliOptions.qvMin = 64;
                cliOptions.qvMax = cliOptions.qvMin + 40;
            } else if (cliOptions.type == "illumina-1.8+") {
                cliOptions.qvMin = 33;
                cliOptions.qvMax = cliOptions.qvMin + 41;
            } else {
                size_t colon = cliOptions.type.find_first_of(':');
                if (colon != std::string::npos) {
                    std::string min = cliOptions.type.substr(0, colon);
                    std::string max = cliOptions.type.substr(colon+1);

                    if (std::all_of(min.begin(), min.end(), ::isdigit) == false) {
                        throwErrorException("Quality value minimum is not numeric");
                    }
                    if (min.empty() == true) {
                        throwErrorException("Quality value minimum is empty");
                    }
                    if (std::all_of(max.begin(), max.end(), ::isdigit) == false) {
                        throwErrorException("Quality value maximum is not numeric");
                    }
                    if (max.empty() == true) {
                        throwErrorException("Quality value maximum is empty");
                    }

                    cliOptions.qvMin = std::stoi(min);
                    cliOptions.qvMax = std::stoi(max);
                } else {
                    throwErrorException("Quality value type not supported");
                }
            }
            LOG("Quality value range: [%d,%d]", cliOptions.qvMin, cliOptions.qvMax);
        }

        // Check command line options for decompression mode
        if (cliOptions.decompress == true) {
            LOG("Input file: %s", cliOptions.inFileName.c_str());
            if (fileNameExtension(cliOptions.inFileName) != std::string("cq")) {
                throwErrorException("Input file extension must be 'cq'");
            }

            if (cliOptions.outFileName.empty()) {
                cliOptions.outFileName.append(cliOptions.inFileName);
                cliOptions.outFileName.append(".qual");
            }
            LOG("Output file: %s", cliOptions.outFileName.c_str());
            if ((fileExists(cliOptions.outFileName) == true) && (cliOptions.force == false)) {
                throwErrorException("Output file already exists (use option 'f' to force overwriting)");
            }

            // Check if the SAM file exists and if it has the correct extension
           LOG("SAM file: %s", cliOptions.samFileName.c_str());
           if (fileExists(cliOptions.samFileName) == false) {
                throwErrorException("Cannot access SAM file");
            }
            if (fileNameExtension(cliOptions.samFileName) != std::string("sam")) {
                throwErrorException("SAM file extension must be 'sam'");
            }
        }

        // Compress or decompress
        if (cliOptions.decompress == false) {
            CalqEncoder calqEncoder(cliOptions);
            calqEncoder.encode();
        } else {
            CalqDecoder calqDecoder(cliOptions);
            calqDecoder.decode();
        }

    } catch (TCLAP::ArgException &tclapException) {
        ERROR("%s (argument: %s)", tclapException.error().c_str(), tclapException.argId().c_str());
        return EXIT_FAILURE;
    } catch (const ErrorException &errorException) {
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


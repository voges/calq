/** @file calq.cc
 *  @brief This file contains the main function of the calq compression tool.
 *
 *  This file contains the main function (and some help functions) of the calq
 *  compression tool.
 *  Encoding: The program takes a SAM file as input and encodes the quality
 *            values. The encoded bitstream is written to a CQ file.
 *  Decoding: The program reads in the compressed CQ file and reconstructs
 *            the quality scores. The quality scores are written to a QUAL
 *            file. The SAM file is -not- reconstructed.
 *
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  2016-05-30: added reference file(s) (voges)
 *  2016-05-29: added TCLAP for command line argument parsing; this btw
 *              is then compatible with Windows, too (as there is no
 *              getopt.h on Windows systems) (voges)
 */

#include "os_config.h"

#ifdef OS_WINDOWS
    #define TCLAP_NAMESTARTSTRING "~~"
    #define TCLAP_FLAGSTARTSTRING "/"
#else
    #define TCLAP_NAMESTARTSTRING "--"
    #define TCLAP_FLAGSTARTSTRING "-"
#endif

#include "CLIOptions.h"
#include "Codecs/CalqCodec.h"
#include "common.h"
#include "cmake_config.h"
#include "Exceptions.h"
#include <iostream>
#include <tclap/CmdLine.h>

CLIOptions cliOptions;

static void printVersion(void)
{
    std::cout << "Program: calq" << std::endl;
    std::cout << "Version: " << VERSION << std::endl;
    std::cout << "Build time: " << TIMESTAMP_UTC << std::endl;
    std::cout << "Git branch: " << GIT_BRANCH << std::endl;
    std::cout << "Git commit hash: " << GIT_COMMIT_HASH_SHORT << std::endl;
    std::cout << std::endl;
}

static void printCopyright(void)
{
    std::cout << "Copyright (c) 2015-" << BUILD_YEAR << std::endl;
    std::cout << "Leibniz Universitaet Hannover, Institut fuer ";
    std::cout << "Informationsverarbeitung (TNT)" << std::endl;
    std::cout << "Contact: Jan Voges <voges@tnt.uni-hannover.de>" << std::endl;
    std::cout << std::endl;
}

int main(int argc, char *argv[])
{
    printVersion();
    printCopyright();

    try {
        // TCLAP class
        TCLAP::CmdLine cmd("calq - Adaptive lossy compression of next-generation sequencing quality values", ' ', VERSION);

        // TCLAP arguments
        TCLAP::SwitchArg decompressSwitch("d", "decompress", "Decompress CQ file", cmd, false);
        TCLAP::SwitchArg forceSwitch("f", "force", "Force overwriting of output files", cmd, false);
        TCLAP::UnlabeledValueArg<std::string> infileArg("infile", "Input file", true, "", "string", cmd);
        TCLAP::ValueArg<std::string> outfileArg("o", "outfile", "Output file", false, "", "string", cmd);
        TCLAP::ValueArg<int> polyploidyArg("p", "polyploidy", "Polyploidy", false, 2, "unsigned int", cmd);
        TCLAP::MultiArg<std::string> referenceArg("r", "reference", "Reference file(s) (FASTA format)", true, "string", cmd);

        // let the TCLAP class parse the provided arguments
        cmd.parse(argc, argv);

        // get the value parsed by each arg
        cliOptions.decompress = decompressSwitch.getValue();
        cliOptions.force = forceSwitch.getValue();
        cliOptions.inFileName = infileArg.getValue();
        cliOptions.outFileName = outfileArg.getValue();
        cliOptions.polyploidy = polyploidyArg.getValue();
        cliOptions.refFileNames = referenceArg.getValue();

        // check polyploidy
        if (cliOptions.polyploidy < 1) {
            throwErrorException("Polyploidy must be greater than 0");
        }
        std::cout << ME << "Using polyploidy: " << cliOptions.polyploidy << std::endl;

        // check if the input file exists and if it has the correct extension
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

        // create correct output file name if it was not provided via the
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

        // check if the reference file(s) exist and check for FASTA file
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

        // compress or decompress
        if (cliOptions.decompress == false) {
            // invoke compressor
            CalqEncoder calqEncoder(cliOptions.inFileName, cliOptions.outFileName, cliOptions.refFileNames, cliOptions.polyploidy);
            calqEncoder.encode();
        } else {
            // invoke decompressor
            CalqDecoder calqDecoder(cliOptions.inFileName, cliOptions.outFileName, cliOptions.refFileNames);
            calqDecoder.decode();
        }
    std::cout << ME << "Finished" << std::endl;
    } catch (TCLAP::ArgException &tclapException) {
        std::cerr << "Argument error: " << tclapException.error() << " (argument: " << tclapException.argId() << ")" << std::endl;
        return EXIT_FAILURE;
    } catch (const ErrorException &errorException) {
        std::cerr << "Error: " << errorException.what() << std::endl;
        return EXIT_FAILURE;
    } catch (const std::exception &stdException) {
        std::cerr << "Fatal error: " << stdException.what() << std::endl;
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Unkown error occured" << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


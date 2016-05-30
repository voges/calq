/** @file calq.cc
 *  @brief This file contains the main function of the calq compression tool.
 *
 *  This file contains the main function (and some help functions) of the calq
 *  compression tool.
 *
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  2016-05-29: added TCLAP for command line argument parsing; this btw
 *              is then compatible with Windows, too (as there is no
 *              getopt.h)
 */

#include "os.h"

#ifdef OS_WINDOWS
    #define TCLAP_NAMESTARTSTRING "~~"
    #define TCLAP_FLAGSTARTSTRING "/"
#else
    //#define TCLAP_NAMESTARTSTRING "--"
    //#define TCLAP_FLAGSTARTSTRING "-"
#endif

#include "CalqCodec.h"
#include "common.h"
#include "config.h"
#include "Exceptions.h"
#include <iostream>
#include <tclap/CmdLine.h>

CLIOptions cliOptions;

static void printVersion(void)
{
    std::cout << "calq" << std::endl;
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

    // init CLI options
    cliOptions.decompress = false;
    cliOptions.force = false;
    cliOptions.inFileName = "";
    cliOptions.outFileName = "";
    cliOptions.referenceFileName = "";

    try {
        // TCLAP arguments definition and parsing
        TCLAP::CmdLine cmd("calq - Lossy compression of next-generation sequencing quality values", ' ', VERSION);

        TCLAP::SwitchArg decompressSwitch("d", "decompress", "Decompress CQ file", cmd, false);
        TCLAP::SwitchArg forceSwitch("f", "force", "Forces overwriting of output files", cmd, false);
        TCLAP::UnlabeledValueArg<std::string> infileArg("infile", "Input file", true, "", "string", cmd);
        TCLAP::ValueArg<std::string> outfileArg("o", "outfile", "Output file", false, "", "string", cmd);
        TCLAP::ValueArg<std::string> referenceArg("r", "reference", "Reference file (FASTA format)", true, "", "string", cmd);

        cmd.parse(argc, argv);

        // get the value parsed by each arg
        cliOptions.force = forceSwitch.getValue();
        cliOptions.inFileName = infileArg.getValue();
        cliOptions.outFileName = outfileArg.getValue();
        cliOptions.decompress = decompressSwitch.getValue();
        cliOptions.referenceFileName = referenceArg.getValue();

        // check the arguments for sanity
        if (!fileExists(cliOptions.inFileName)) {
            throwUserException("Cannot access input file");
        }
        if (!fileExists(cliOptions.referenceFileName)) {
            throwUserException("Cannot access reference file");
        }

        if (!cliOptions.decompress) {
            // check for correct infile extension
            if (filenameExtension(cliOptions.inFileName) != std::string("sam")) {
                throwUserException("Input file extension must be 'sam'");
            }

            // create correct output file name if it was not provided via the
            // command line options
            if (cliOptions.outFileName.empty()) {
                cliOptions.outFileName.append(cliOptions.inFileName);
                cliOptions.outFileName.append(".cq");
            }

            // check if output file is already there and if the user wants to
            // overwrite it in this case
            if (fileExists(cliOptions.outFileName) && cliOptions.force == false) {
                std::cout << "Output file already exists: " << cliOptions.outFileName << std::endl;
                std::cout << "Do you want to overwrite it? ";
                if (!yesno()) {
                    throwUserException("Exited because we do not overwrite the output file");
                }
            }

            // invoke compressor
            std::cout << "Compressing: " << cliOptions.inFileName << " > " << cliOptions.outFileName << std::endl;
            //CalqEncoder calqEncoder(cliOptions.inFileName, cliOptions.outFileName, cliOptions.referenceFileName);
            //calqEncoder.encode();
            std::cout << "Finished compression" << std::endl;
        } else {
            // check for correct infile extension
            if (filenameExtension(cliOptions.inFileName) != std::string("cq")) {
                throwUserException("Input file extension must be 'cq'");
            }

            // create correct output file name if it was not provided via the
            // command line options
            if (cliOptions.outFileName.empty()) {
                cliOptions.outFileName.append(cliOptions.inFileName);
                cliOptions.outFileName.append(".sam");
            }

            // check if output file is already there and if the user wants to
            // overwrite it in this case
            if (fileExists(cliOptions.outFileName) && cliOptions.force == false) {
                std::cout << "Output file already exists: " << cliOptions.outFileName << std::endl;
                std::cout << "Do you want to overwrite it? ";
            if (!yesno()) {
                    throwUserException("Exited because we do not overwrite the output file");
                }
            }

            // invoke decompressor
            std::cout << "Decompressing: " << cliOptions.inFileName << " > " << cliOptions.outFileName << std::endl;
            //CalqDecoder calqDecoder(cliOptions.inFileName, cliOptions.outFileName, cliOptions.referenceFileName);
            //calqDecoder.decode();
            std::cout << "Finished decompression" << std::endl;
        }
    } catch (TCLAP::ArgException &tclapException) {
        std::cerr << "Error: " << tclapException.error() << " for argument " << tclapException.argId() << std::endl;
        return EXIT_FAILURE;
    } catch (const UserException &userException) {
        std::cerr << userException.what() << std::endl << std::endl;
        return EXIT_FAILURE;
    } catch (const ErrorException &errorException) {
        std::cerr << "Error: " << errorException.what() << std::endl;
        return EXIT_FAILURE;
    } catch (const std::exception &stdException) {
        std::cerr << "Fatal error: " << stdException.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


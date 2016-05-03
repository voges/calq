/** @file calq.cc
 *  @brief This file contains the main function of the calq compression tool.
 *
 *  This file contains the main function (and some static helpers) of the calq
 *  compression tool. It furthermore defines the global object 'cliOptions'
 *  which contains the command line options that are therefore accessible to
 *  every linked module.
 *
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#include "CalqCodec.h"
#include "config.h"
#include "Exceptions.h"
#include <getopt.h>
#include <iostream>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

typedef struct {
    size_t blockSize;
    bool force;
    std::string infileName;
    enum {
        MODE_COMPRESS,
        MODE_DECOMPRESS,
        MODE_INFO 
    } mode;
    std::string outfileName;
} CLIOptions;

CLIOptions cliOptions;

static void printCopyright(void)
{
    std::cout << "Copyright (c) 2015-" << BUILD_YEAR << std::endl;
    std::cout << "Leibniz Universitaet Hannover, Institut fuer ";
    std::cout << "Informationsverarbeitung (TNT)" << std::endl;
    std::cout << "Contact: Jan Voges <voges@tnt.uni-hannover.de>" << std::endl;
}

static void printVersion(void)
{
    std::cout << "calq " << VERSION_MAJOR << "." << VERSION_MINOR << "." << VERSION_PATCH << std::endl;
    std::cout << "Build time: " << UTCTIMESTAMP << std::endl;
    std::cout << "Git revision: " << GITREVISION_LONG << std::endl;
    std::cout << std::endl;
    printCopyright();
}

static void printHelp(void)
{
    printVersion();
    std::cout << std::endl;
    std::cout << "Usage:" << std::endl;
    std::cout << "  Compress  : calq [-o FILE] [-b SIZE] [-f] file.sam" << std::endl;
    std::cout << "  Decompress: calq -d [-o FILE] [-f] file.cq" << std::endl;
    std::cout << "  Info      : calq -i file.cq" << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "  -b  --blocksz=SIZE Specify block SIZE" << std::endl;
    std::cout << "  -d  --decompress   Decompress" << std::endl;
    std::cout << "  -f, --force        Force overwriting of output file(s)" << std::endl;
    std::cout << "  -h, --help         Print this help" << std::endl;
    std::cout << "  -i, --info         Print information about CQ file" << std::endl;
    std::cout << "  -o, --output=FILE  Specify output FILE" << std::endl;
    std::cout << "  -v, --version      Display program version" << std::endl;
    std::cout << std::endl;
}

static void parseOptions(int argc, char *argv[])
{
    int opt;

    static struct option longOptions[] = {
        { "blocksz",    required_argument, NULL, 'b'},
        { "decompress", no_argument,       NULL, 'd'},
        { "force",      no_argument,       NULL, 'f'},
        { "help",       no_argument,       NULL, 'h'},
        { "info",       no_argument,       NULL, 'i'},
        { "output",     required_argument, NULL, 'o'},
        { "version",    no_argument,       NULL, 'v'},
        { NULL,         0,                 NULL,  0 }
    };

    const char *shortOptions = "b:dfhio:v";

    do {
        int optIdx = 0;
        opt = getopt_long(argc, argv, shortOptions, longOptions, &optIdx);
        switch (opt) {
        case -1:
            break;
        case 'b':
            if (atoi(optarg) <= 0) {
                throwUserException("Block size must be positive");
            } else {
                cliOptions.blockSize = (size_t)atoi(optarg);
            }
            break;
        case 'd':
            if (cliOptions.mode == CLIOptions::MODE_INFO) {
                throwUserException("Cannot decompress and get info at once");
            } else {
                cliOptions.mode = CLIOptions::MODE_DECOMPRESS;
            }
            break;
        case 'f':
            cliOptions.force = true;
            break;
        case 'h':
            printHelp();
            exit(EXIT_SUCCESS);
            break;
        case 'i':
            if (cliOptions.mode == CLIOptions::MODE_DECOMPRESS) {
                throwUserException("Cannot decompress and get info at once");
            } else {
                cliOptions.mode = CLIOptions::MODE_INFO;
            }
            break;
        case 'o':
            cliOptions.outfileName = optarg;
            break;
        case 'v':
            printVersion();
            exit(EXIT_SUCCESS);
            break;
        default:
            throwUserException("Unknown option(s)");
        }
    } while (opt != -1);

    // the input file must be the one remaining command line argument
    if (argc - optind > 1) {
        throwUserException("Only one input file allowed");
    } else if (argc - optind < 1) {
        throwUserException("Input file missing");
    } else {
        cliOptions.infileName = argv[optind];
    }

    // sanity checks
    if (cliOptions.mode == CLIOptions::MODE_COMPRESS) {
        // all possible options are legal in compression mode
        if (cliOptions.blockSize == 0) {
            std::cout << "Using default block size 10,000" << std::endl;
            cliOptions.blockSize = 10000; // default value
        }
    } else if (cliOptions.mode == CLIOptions::MODE_DECOMPRESS) {
        // option -b is illegal in decompression mode
        if (cliOptions.blockSize != 0) {
            throwUserException("Option -b is illegal in decompression mode");
        }
    } else { // CLIOptions::MODE_INFO
        // options -bf are illegal in info mode
        if (cliOptions.blockSize != 0) {
            throwUserException("Option -b is illegal in info mode");
        }
        if (cliOptions.force == true) {
            throwUserException("Option -f is illegal in info mode");
        }
    }
}

static const char * filenameExtension(const std::string &filename)
{
    const char *dot = strrchr(filename.c_str(), '.');
    if (!dot || dot == filename.c_str()) { return ""; }
    return (dot + 1);
}

static bool fileExists(const std::string &filename)
{
    return !access(filename.c_str(), F_OK | R_OK);
}

static bool yesno(void)
{
    int c = getchar();
    bool yes = c == 'y' || c == 'Y';
    while (c != '\n' && c != EOF) {
        c = getchar();
    }
    return yes;
}

static void handleSignal(const int sig)
{
    signal(sig, SIG_IGN); // ignore the signal
    std::cout << "Catched signal: " << sig << std::endl;
    signal(sig, SIG_DFL); // invoke default signal action
    raise(sig);
}

int main(int argc, char *argv[])
{
    // init CLI options
    cliOptions.blockSize = 0;
    cliOptions.infileName = "";
    cliOptions.outfileName = "";
    cliOptions.force = false;
    cliOptions.mode = CLIOptions::MODE_COMPRESS;

    // register custom signal handler(s)
    signal(SIGHUP,  handleSignal);
    signal(SIGQUIT, handleSignal);
    signal(SIGABRT, handleSignal);
    signal(SIGPIPE, handleSignal);
    signal(SIGTERM, handleSignal);
    signal(SIGXCPU, handleSignal);
    signal(SIGXFSZ, handleSignal);

    try {
        parseOptions(argc, argv);

        if (!fileExists(cliOptions.infileName)) {
            throwUserException("Cannot access input file");
        }

        switch (cliOptions.mode) {
        case CLIOptions::MODE_COMPRESS: {
            // check for correct infile extension
            if (filenameExtension(cliOptions.infileName) != std::string("sam")) {
                throwUserException("Input file extension must be 'sam'");
            }

            // create correct output file name
            if (cliOptions.outfileName.empty()) {
                cliOptions.outfileName.append(cliOptions.infileName);
                cliOptions.outfileName.append(".cq");
            }

            // check if output file is already there and if the user wants to
            // overwrite it in this case
            if (fileExists(cliOptions.outfileName) && cliOptions.force == false) {
                std::cout << "Output file already exists: " << cliOptions.outfileName << std::endl;
                std::cout << "Do you want to overwrite it? ";
                if (!yesno()) {
                    throwUserException("Exited because we do not overwrite the output file");
                }
            }

            // invoke compressor
            std::cout << "Compressing: " << cliOptions.infileName << " > " << cliOptions.outfileName << std::endl;
            CalqEncoder calqEncoder(cliOptions.infileName, cliOptions.outfileName, cliOptions.blockSize);
            calqEncoder.encode();
            std::cout << "Finished compression" << std::endl;

            break;
        }
        case CLIOptions::MODE_DECOMPRESS: {
            // check for correct infile extension
            if (filenameExtension(cliOptions.infileName) != std::string("cq")) {
                throwUserException("Input file extension must be 'cq'");
            }

            // create correct output file name
            if (cliOptions.outfileName.empty()) {
                cliOptions.outfileName.append(cliOptions.infileName);
                cliOptions.outfileName.append(".sam");
            }

            // check if output file is already there and if the user wants to
            // overwrite it in this case
            if (fileExists(cliOptions.outfileName) && cliOptions.force == false) {
                std::cout << "Output file already exists: " << cliOptions.outfileName << std::endl;
                std::cout << "Do you want to overwrite it? ";
                if (!yesno()) {
                    throwUserException("Exited because we do not overwrite the output file");
                }
            }

            // invoke decompressor
            std::cout << "Decompressing: " << cliOptions.infileName << " > " << cliOptions.outfileName << std::endl;
            CalqDecoder calqDecoder(cliOptions.infileName, cliOptions.outfileName);
            calqDecoder.decode();
            std::cout << "Finished decompression" << std::endl;

            break;
        }
        case CLIOptions::MODE_INFO: {
            // check for correct infile extension
            if (filenameExtension(cliOptions.infileName) != std::string("cq")) {
                throwUserException("Input file extension must be 'cq'");
            }

            // invoke info tool
            std::cout << "Reading information: " << cliOptions.infileName << std::endl;
            CalqInfoTool calqInfoTool(cliOptions.infileName);
            calqInfoTool.extractInfo();
            std::cout << "Finished reading information" << std::endl;

            break;
        }
        default:
            throwErrorException("Unknown mode");
        }
    } catch (ErrorException &errorException) {
        std::cerr << "Error: " << errorException.what() << std::endl;
    } catch (UserException &userException) {
        std::cerr << userException.what() << std::endl;
    }

    return EXIT_SUCCESS;
}


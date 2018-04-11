#if 0
/** @file calq.cc
 *  @brief This file contains the main function of CALQ.
 */

// Copyright 2015-2017 Leibniz Universitaet Hannover

#include "CalqEncoder.h"
#include "CalqDecoder.h"
#include "cmake.h"
#include "Common/Options.h"
#include "Common/Exceptions.h"
#include "Common/helpers.h"
#include "Common/log.h"
#include "tclap/CmdLine.h"

static void printVersionAndCopyright(void) {
    printf("-----------------------------------------------\n");
    printf("Program: CALQ\n");
    printf("Version: %s\n", CALQ_VERSION);
    printf("Build time: %s\n", CALQ_TIMESTAMP_UTC);
    printf("Git commit hash: %s\n", CALQ_GIT_COMMIT_HASH_SHORT);
    printf("-----------------------------------------------\n");
    printf("Copyright (c) 2015-%d\n", CALQ_BUILD_YEAR);
    printf("Leibniz Universitaet Hannover\n");
    printf("Institut fuer Informationsverarbeitung (TNT)\n");
    printf("Contact: Jan Voges <voges@tnt.uni-hannover.de>\n");
    printf("-----------------------------------------------\n");
    printf("Please report bugs to voges@tnt.uni-hannover.de\n");
    printf("-----------------------------------------------\n");
}

int main(int argc, char *argv[]) {
    printVersionAndCopyright();

    try {
        calq::Options options;

        // TCLAP class
        TCLAP::CmdLine cmd("CALQ", ' ', CALQ_VERSION);

        // TCLAP arguments (both compression and decompression)
        TCLAP::SwitchArg forceSwitch("f", "force", "Force overwriting of output files etc.", cmd, false);
        TCLAP::UnlabeledValueArg<std::string> inputFileNameArg("inputFileName", "Input file name", true, "", "string", cmd);
        TCLAP::ValueArg<std::string> outputFileNameArg("o", "outputFileName", "Output file name", false, "", "string", cmd);

        // TCLAP arguments (only compression)
        TCLAP::ValueArg<int> blockSizeArg("b", "blockSize", "Block size (in number of SAM records)", false, 10000, "int", cmd);
        TCLAP::ValueArg<int> polyploidyArg("p", "polyploidy", "Polyploidy", false, 2, "int", cmd);
        TCLAP::ValueArg<std::string> qualityValueTypeArg("q", "qualityValueType", "Quality value type (Sanger: Phred+33 [0,40]; Illumina-1.3+: Phred+64 [0,40]; Illumina-1.5+: Phred+64 [0,40]; Illumina-1.8+: Phred+33 [0,41]; Max33: Phred+33 [0,93]; Max64: Phred+64 [0,62])", false, "Illumina-1.8+", "string", cmd);
        TCLAP::MultiArg<std::string> referenceFileNamesArg("r", "referenceFileNames", "Reference file name(s) (FASTA format)", false, "string", cmd);

        // TCLAP arguments (only decompression)
        TCLAP::SwitchArg decompressSwitch("d", "decompress", "Decompress", cmd, false);
        TCLAP::ValueArg<std::string> sideInformationFileNameArg("s", "sideInformationFileName", "Side information file name", false, "", "string", cmd);

        // Let the TCLAP class parse the provided arguments
        cmd.parse(argc, argv);

        // Check for sanity in compression mode
        if (decompressSwitch.isSet() == false) {
//             if (referenceFileNamesArg.isSet() == false) {
//                throwErrorException("Argument 'r' required in compression mode");
//             }
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
            if (qualityValueTypeArg.isSet() == true) {
                throwErrorException("Argument 'q' forbidden in decompression mode");
            }
            if (referenceFileNamesArg.isSet() == true) {
                throwErrorException("Argument 'r' forbidden in decompression mode");
            }
            if (sideInformationFileNameArg.isSet() == false) {
                throwErrorException("Argument 's' required in decompression mode");
            }
        }

        // Get the value parsed by each arg
        options.force = forceSwitch.getValue();
        options.inputFileName = inputFileNameArg.getValue();
        options.outputFileName = outputFileNameArg.getValue();
        options.blockSize = blockSizeArg.getValue();
        options.polyploidy = polyploidyArg.getValue();
        options.qualityValueType = qualityValueTypeArg.getValue();
        options.referenceFileNames = referenceFileNamesArg.getValue();
        options.decompress = decompressSwitch.getValue();
        options.sideInformationFileName = sideInformationFileNameArg.getValue();

        // Check the options
        options.validate();

        // Compress or decompress
        if (options.decompress == false) {
            calq::CalqEncoder calqEncoder(options);
            calqEncoder.encode();
            CALQ_LOG("Finished encoding");
        } else {
            calq::CalqDecoder calqDecoder(options);
            calqDecoder.decode();
            CALQ_LOG("Finished decoding");
        }
    } catch (TCLAP::ArgException &tclapException) {
        CALQ_ERROR("%s (argument: %s)", tclapException.error().c_str(), tclapException.argId().c_str());
        return EXIT_FAILURE;
    } catch (const calq::ErrorException &errorException) {
        CALQ_ERROR("%s", errorException.what());
        return EXIT_FAILURE;
    } catch (const std::exception &stdException) {
        CALQ_ERROR("Fatal: %s", stdException.what());
        return EXIT_FAILURE;
    } catch (...) {
        CALQ_ERROR("Unkown error occured");
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
#else

#include "QualCodec/HaplotyperTest.h"

int main(int argc, char *argv[]) {
    haplotyperTest();
}

#endif


/** @file calq.cc
 *  @brief This file contains the main function of CALQ.
 */

// Copyright 2015-2017 Leibniz Universitaet Hannover

#include "CalqEncoder.h"
#include "CalqDecoder.h"
#include "Common/ErrorExceptionReporter.h"
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
        TCLAP::SwitchArg debugSwitch("", "debug", "Output verbose debug information regarding activity profile.", cmd, false);
        TCLAP::SwitchArg testSwitch("", "test", "Run test cases.", cmd, false);
        TCLAP::UnlabeledValueArg<std::string> inputFileNameArg("inputFileName", "Input file name", true, "", "string", cmd);
        TCLAP::ValueArg<std::string> outputFileNameArg("o", "outputFileName", "Output file name", false, "", "string", cmd);

        // TCLAP arguments (only compression)
        TCLAP::ValueArg<int> blockSizeArg("b", "blockSize", "Block size (in number of SAM records). Default 10000", false, 10000, "int", cmd);
        TCLAP::ValueArg<int> filterSizeArg("", "filterSize", "Haplotyper filter radius. Default 17. (v2 only)", false, 17, "int", cmd);
        TCLAP::ValueArg<int> quantizationMinArg("", "quantizationMin", "Minimum quantization steps. Default 2", false, 2, "int", cmd);
        TCLAP::ValueArg<int> quantizationMaxArg("", "quantizationMax", "Maximum quantization steps. Default 8", false, 8, "int", cmd);
        TCLAP::ValueArg<int> polyploidyArg("p", "polyploidy", "Polyploidy. Default 2", false, 2, "int", cmd);
        TCLAP::ValueArg<std::string> qualityValueTypeArg("q", "qualityValueType", "Quality value type (Sanger: Phred+33 [0,40]; Illumina-1.3+: Phred+64 [0,40]; Illumina-1.5+: Phred+64 [0,40]; Illumina-1.8+: Phred+33 [0,41]; Max33: Phred+33 [0,93]; Max64: Phred+64 [0,62])", false, "Illumina-1.8+", "string", cmd);
        TCLAP::ValueArg<std::string> referenceFileNamesArg("r", "referenceFileName", "Reference file (FASTA format) (v2 only)", false, "", "string", cmd);
        TCLAP::ValueArg<std::string> filterTypeArg("", "filterType", "Haplotyper Filter Type (Gauss; Rectangle). Default Gauss. (v2 only)", false, "Gauss", "string", cmd);
        TCLAP::ValueArg<std::string> quantizerTypeArg("", "quantizerType", "Quantizer type (Uniform; Lloyd). Default Uniform", false, "Uniform", "string", cmd);
        TCLAP::ValueArg<std::string> versionArg("", "CALQ-Version", "v1 or v2. Default v1", false, "v1", "string", cmd);
        TCLAP::SwitchArg noSquash("", "noSquash", "Do not squash activityscores between 0.0 and 1.0, which is done by default", cmd, false);

        // TCLAP arguments (only decompression)
        TCLAP::SwitchArg decompressSwitch("d", "decompress", "Decompress", cmd, false);
        TCLAP::ValueArg<std::string> sideInformationFileNameArg("s", "sideInformationFileName", "Side information file name", false, "", "string", cmd);

        // Let the TCLAP class parse the provided arguments
        cmd.parse(argc, argv);

        // Check for sanity in compression mode
        if (decompressSwitch.isSet() == false) {
            if (referenceFileNamesArg.isSet() == false) {
                throwErrorException("Argument 'r' required in compression mode");
            }
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
            if (filterSizeArg.isSet() == true) {
                throwErrorException("Argument 'filterSize' forbidden in decompression mode");
            }
            if (filterTypeArg.isSet() == true) {
                throwErrorException("Argument 'filterType' forbidden in decompression mode");
            }
            if (noSquash.isSet() == true) {
                throwErrorException("Argument 'noSquash' forbidden in decompression mode");
            }
        }

        // Get the value parsed by each arg
        options.force = forceSwitch.getValue();
        options.debug = debugSwitch.getValue();
        options.test = testSwitch.getValue();
        options.inputFileName = inputFileNameArg.getValue();
        options.outputFileName = outputFileNameArg.getValue();
        options.blockSize = blockSizeArg.getValue();
        options.filterSize = filterSizeArg.getValue();
        options.quantizationMin = quantizationMinArg.getValue();
        options.quantizationMax = quantizationMaxArg.getValue();
        options.polyploidy = polyploidyArg.getValue();
        options.qualityValueType = qualityValueTypeArg.getValue();
        options.filterTypeStr = filterTypeArg.getValue();
        options.quantizerTypeStr = quantizerTypeArg.getValue();
        options.referenceFileNames = referenceFileNamesArg.getValue();
        options.decompress = decompressSwitch.getValue();
        options.sideInformationFileName = sideInformationFileNameArg.getValue();
        options.squash = !noSquash.getValue();
        options.versionStr = versionArg.getValue();

        // Check the options
        options.validate();

        // Compress or decompress
        if (options.test) {
            calq::haplotyperTest();
            CALQ_LOG("Finished testing");
        } else if (options.decompress == false) {
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



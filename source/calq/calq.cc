#include "calq/calq_encoder.h"
#include "calq/calq_decoder.h"
#include "calq/error_exception_reporter.h"
#include "calq/tclap/CmdLine.h"
#include "calq/log.h"

#include "calq/SAMFileHandler.h"

int main(int argc, char* argv[]) {
    try {
        calq::Options options;

        // TCLAP class
        TCLAP::CmdLine cmd("CALQ", ' ', "");

        // TCLAP arguments (both compression and decompression)
        TCLAP::SwitchArg forceSwitch("f", "force", "Force overwriting of output files etc.", cmd, false);
        TCLAP::SwitchArg debugSwitch("", "debug", "Output verbose debug information regarding activity profile.", cmd, false);
        TCLAP::UnlabeledValueArg<std::string> inputFileNameArg("inputFileName", "Input file name", true, "", "string", cmd);
        TCLAP::ValueArg<std::string> outputFileNameArg("o", "outputFileName", "Output file name", false, "", "string", cmd);

        // TCLAP arguments (only compression)
        TCLAP::ValueArg<size_t> blockSizeArg("b", "blockSize", "Block size (in number of SAM records). Default 10000", false, 10000, "int", cmd);
        TCLAP::ValueArg<size_t> filterSizeArg("", "filterSize", "Haplotyper filter radius. Default 17. (v2 only)", false, 17, "int", cmd);
        TCLAP::ValueArg<size_t> quantizationMinArg("", "quantizationMin", "Minimum quantization steps. Default 2", false, 2, "int", cmd);
        TCLAP::ValueArg<size_t> quantizationMaxArg("", "quantizationMax", "Maximum quantization steps. Default 8", false, 8, "int", cmd);
        TCLAP::ValueArg<size_t> polyploidyArg("p", "polyploidy", "Polyploidy. Default 2", false, 2, "int", cmd);
        TCLAP::ValueArg<std::string> qualityValueTypeArg("q", "qualityValueType", "Quality value type (Sanger: Phred+33 [0,40]; "
                                                                                  "Illumina-1.3+: Phred+64 [0,40]; Illumina-1.5+: Phred+64 [0,40]; "
                                                                                  "Illumina-1.8+: Phred+33 [0,41]; Max33: Phred+33 [0,93]; Max64: "
                                                                                  "Phred+64 [0,62])", false, "Illumina-1.8+", "string", cmd);
        TCLAP::ValueArg<std::string> referenceFileNamesArg("r", "referenceFileName", "Reference file (FASTA format) (v2 only)", false, "", "string", cmd);
        TCLAP::ValueArg<std::string> filterTypeArg("", "filterType", "Haplotyper Filter Type (Gauss; Rectangle). Default Gauss. (v2 only)",
                                                   false, "Gauss", "string", cmd);
        TCLAP::ValueArg<std::string> quantizerTypeArg("", "quantizerType", "Quantizer type (Uniform; Lloyd). Default Uniform", false, "Uniform", "string", cmd);
        TCLAP::ValueArg<std::string> versionArg("", "CALQ-Version", "v1 or v2. Default v1", false, "v1", "string", cmd);
        TCLAP::SwitchArg noSquash("", "noSquash", "Do not squash activityscores between 0.0 and 1.0, which is done by default", cmd, false);

        // TCLAP arguments (only decompression)
        TCLAP::SwitchArg decompressSwitch("d", "decompress", "Decompress", cmd, false);
        TCLAP::ValueArg<std::string> sideInformationFileNameArg("s", "sideInformationFileName", "Side information file name", false, "", "string", cmd);

        // Let the TCLAP class parse the provided arguments
        cmd.parse(argc, argv);

        // Check for sanity in compression mode
        if (!decompressSwitch.isSet()) {
            if (!referenceFileNamesArg.isSet()) {
                throwErrorException("Argument 'r' required in compression mode");
            }
            if (sideInformationFileNameArg.isSet()) {
                throwErrorException("Argument 's' forbidden in compression mode");
            }
        }

        // Check for sanity in decompression mode
        if (decompressSwitch.isSet()) {
            if (blockSizeArg.isSet()) {
                throwErrorException("Argument 'b' forbidden in decompression mode");
            }
            if (polyploidyArg.isSet()) {
                throwErrorException("Argument 'p' forbidden in decompression mode");
            }
            if (qualityValueTypeArg.isSet()) {
                throwErrorException("Argument 'q' forbidden in decompression mode");
            }
            if (referenceFileNamesArg.isSet()) {
                throwErrorException("Argument 'r' forbidden in decompression mode");
            }
            if (!sideInformationFileNameArg.isSet()) {
                throwErrorException("Argument 's' required in decompression mode");
            }
            if (filterSizeArg.isSet()) {
                throwErrorException("Argument 'filterSize' forbidden in decompression mode");
            }
            if (filterTypeArg.isSet()) {
                throwErrorException("Argument 'filterType' forbidden in decompression mode");
            }
            if (noSquash.isSet()) {
                throwErrorException("Argument 'noSquash' forbidden in decompression mode");
            }
        }

        // Get the value parsed by each arg
        options.force = forceSwitch.getValue();
        options.debug = debugSwitch.getValue();
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
        if (!options.decompress) {
            calq::SAMFileHandler sH(options.inputFileName);
            while (sH.readBlock(options.blockSize) != 0) {

                // Container to be filled by lib
                std::vector<uint8_t> quantizerIndices;
                std::vector<std::vector<uint8_t>> stepIndices;
                std::vector<std::vector<uint8_t>> codeBooks;

                std::string test = sH.getSequences().at(0);
                std::cout << test << std::endl;

                // To-DO: Implement this:
                /* calqEncoder.encode( sH.getPositions(),
                                    sH.getSequences(),
                                    sH.getCigars(),
                                    sH.getQualityScores(),
                                    Other Params,
                                    &quantizerIndices,
                                    &codeBooks);
                */

            // calq::CalqEncoder calqEncoder(options);
            // calqEncoder.encode();
            }
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
        CALQ_ERROR("Unknown error occurred");
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

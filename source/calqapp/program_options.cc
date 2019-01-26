#include "calqapp/program_options.h"

// -----------------------------------------------------------------------------

#include "calqapp/error_reporter.h"
#include "calqapp/helpers.h"
#include "calqapp/logging.h"

// -----------------------------------------------------------------------------

#include <boost/program_options.hpp>

// -----------------------------------------------------------------------------

namespace calqapp {

// -----------------------------------------------------------------------------

ProgramOptions::ProgramOptions(
        int argc,
        char *argv[]
)
        : force(),
        test(false),
        help(false),
        inputFilePath(),
        outputFilePath(),
        blockSize(),
        qualityValueType(),
        referenceFilePath(),
        filterTypeStr(),
        quantizerTypeStr(),
        versionStr(),
        decompress(),
        sideInformationFilePath(){
    processCommandLine(argc, argv);
}

// -----------------------------------------------------------------------------

ProgramOptions::~ProgramOptions() = default;

// -----------------------------------------------------------------------------

void ProgramOptions::validateCompress(){
    if (calqapp::fileNameExtension(inputFilePath) != std::string("sam")) {
        throwErrorException("Input file name extension must be 'sam'");
    }

    CALQ_LOG("Block size: %d", static_cast<int>(blockSize));
    if (blockSize < 1) {
        throwErrorException("Block size must be greater than 0");
    }

    CALQ_LOG("Quantization min steps: %d",
             static_cast<int>(options.quantizationMin));
    if (options.quantizationMin < 2) {
        throwErrorException("Quantization must be greater than 1");
    }

    CALQ_LOG("Quantization max steps: %d",
             static_cast<int>(options.quantizationMax));
    if (options.quantizationMax < 2
        || options.quantizationMax < options.quantizationMin) {
        throwErrorException(
                "Quantization must be greater than"
                " 1 and quantizationMin"
        );
    }

    CALQ_LOG("Polyploidy: %d", static_cast<int>(options.polyploidy));
    if (options.polyploidy < 1) {
        throwErrorException("Polyploidy must be greater than 0");
    }

    CALQ_LOG("Quantizer type: %s", quantizerTypeStr.c_str());
    if (quantizerTypeStr == "Uniform") {
        options.quantizerType =
                calq::QuantizerType::UNIFORM;
    } else if (quantizerTypeStr == "Lloyd") {
        options.quantizerType =
                calq::QuantizerType::LLOYD_MAX;
    } else {
        throwErrorException("Quantizer type not supported");
    }

    CALQ_LOG("Quality value type: %s", qualityValueType.c_str());
    if (qualityValueType == "Sanger") {
        // Sanger: Phred+33 [0,40]
        options.qualityValueOffset = 33;
        options.qualityValueMin = 0;
        options.qualityValueMax = 40;
    } else if (qualityValueType == "Illumina-1.3+") {
        // Illumina 1.3+: Phred+64 [0,40]
        options.qualityValueOffset = 64;
        options.qualityValueMin = 0;
        options.qualityValueMax = 40;
    } else if (qualityValueType == "Illumina-1.5+") {
        // Illumina 1.5+: Phred+64 [0,40] with 0=unused,
        // 1=unused, 2=Read Segment Quality Control Indicator ('B')
        CALQ_LOG("Warning: Read Segment Quality Control Indicator"
                 " will not be treated specifically by CALQ"
        );
        options.qualityValueOffset = 64;
        options.qualityValueMin = 0;
        options.qualityValueMax = 40;
    } else if (qualityValueType == "Illumina-1.8+") {
        // Illumina 1.8+ Phred+33 [0,41]
        options.qualityValueOffset = 33;
        options.qualityValueMin = 0;
        options.qualityValueMax = 41;
    } else if (qualityValueType == "Max33") {
        // Max33 Phred+33 [0,93]
        options.qualityValueOffset = 33;
        options.qualityValueMin = 0;
        options.qualityValueMax = 93;
    } else if (qualityValueType == "Max64") {
        // Max64 Phred+64 [0,62]
        options.qualityValueOffset = 64;
        options.qualityValueMin = 0;
        options.qualityValueMax = 62;
    } else {
        throwErrorException("Quality value type not supported");
    }
    CALQ_LOG("Quality value offset: %d",
             static_cast<int>(options.qualityValueOffset)
    );
    CALQ_LOG("Quality value range: [%d,%d]",
             static_cast<int>(options.qualityValueMin),
             static_cast<int>(options.qualityValueMax));

    if (options.debugPileup) {
        CALQ_LOG("Outputting full pileup as debug information");
    }

    if (debugStreams) {
        CALQ_LOG("Outputting full streams as debug information");
    }
}

// -----------------------------------------------------------------------------

void ProgramOptions::validateV1(){
    options.version = calq::Version::V1;
    CALQ_LOG("Compressing using CALQ version 1");
}

// -----------------------------------------------------------------------------

void ProgramOptions::validateV2(){
    options.version = calq::Version::V2;
    CALQ_LOG("Compressing using CALQ version 2");

    if (options.squash) {
        CALQ_LOG("Acitivity scores are squashed between 0.0 and 1.0");
    } else {
        CALQ_LOG("Acitivity scores are !NOT! squashed between 0.0 and 1.0");
    }

    CALQ_LOG("Filter size: %d", static_cast<int>(options.filterSize));
    if (options.filterSize < 1) {
        throwErrorException("Filter size must be greater than 0");
    }

    CALQ_LOG("hqSoftClipPropagation: %d",
             static_cast<int>(options.hqSoftClipPropagation)
    );

    CALQ_LOG("hqSoftClipStreak: %d",
             static_cast<int>(options.hqSoftClipStreak)
    );

    CALQ_LOG("hqSoftClipThreshold: %d",
             static_cast<int>(options.hqSoftClipThreshold)
    );

    CALQ_LOG("filterCutOff: %d",
             static_cast<int>(options.filterCutOff)
    );

    CALQ_LOG("Filter type: %s", filterTypeStr.c_str());
    if (filterTypeStr == "Gauss") {
        options.filterType = calq::FilterType::GAUSS;
    } else if (filterTypeStr == "Rectangle") {
        options.filterType = calq::FilterType::RECTANGLE;
    } else {
        throwErrorException("Filter type not supported");
    }

    CALQ_LOG("  %s", referenceFilePath.c_str());
    if (referenceFilePath.empty()) {
        throwErrorException("Reference file name not provided");
    }
    if (!calqapp::fileExists(referenceFilePath)) {
        throwErrorException("Cannot access reference file");
    }
    if (calqapp::fileNameExtension(referenceFilePath)
        != std::string("fa")
        && calqapp::fileNameExtension(referenceFilePath)
           != std::string("fasta")) {
        throwErrorException(
                "Reference file name extension must "
                "be 'fa' or 'fasta'"
        );
    }
}

// -----------------------------------------------------------------------------

void ProgramOptions::validateDecompress(){
    CALQ_LOG("Decompressing");

    if (calqapp::fileNameExtension(inputFilePath) != std::string("cq")) {
        CALQ_LOG("Warning: Input file name extension is not 'cq'");
    }

    CALQ_LOG("Side information file name: %s",
             sideInformationFilePath.c_str()
    );
    if (sideInformationFilePath.empty()) {
        throwErrorException("No side information file name provided");
    }
    if (calqapp::fileNameExtension(sideInformationFilePath)
        != std::string("sam")) {
        throwErrorException(
                "Side information file name "
                "extension must be 'sam'"
        );
    }
    if (!calqapp::fileExists(sideInformationFilePath)) {
        throwErrorException("Cannot access side information file");
    }
}

// -----------------------------------------------------------------------------

void ProgramOptions::validateCommon(){
    this->options.squash = !this->options.squash;
    if (this->quantizationMin > 255) {
        throwErrorException("Option quantizationMin too big");
    }
    this->options.quantizationMin = uint8_t(this->quantizationMin);

    if (this->quantizationMax > 255) {
        throwErrorException("Option quantizationMax too big");
    }
    this->options.quantizationMax = uint8_t(this->quantizationMax);

    if (this->polyploidy > 255) {
        throwErrorException("Option polyploidy too big");
    }
    this->options.polyploidy = uint8_t(this->polyploidy);

    if (this->hqSoftClipThreshold > 255) {
        throwErrorException("Option hqSoftClipThreshold too big");
    }
    this->options.hqSoftClipThreshold = uint8_t(this->hqSoftClipThreshold);



    // inputFileName
    CALQ_LOG("Input file name: %s", inputFilePath.c_str());
    if (inputFilePath.empty()) {
        throwErrorException("No input file name provided");
    }

    if (!calqapp::fileExists(inputFilePath)) {
        throwErrorException("Cannot access input file");
    }

    CALQ_LOG("Output file name: %s", outputFilePath.c_str());
    if (calqapp::fileExists(outputFilePath)) {
        if (!force) {
            throwErrorException(
                    "Not overwriting output file "
                    "(use option 'f' to force overwriting)"
            );
        } else {
            CALQ_LOG("Force switch set - overwriting output file(s)");
        }
    }
}

// -----------------------------------------------------------------------------

void ProgramOptions::validate(){
    validateCommon();

    if (decompress) {
        validateDecompress();
    } else {
        if (versionStr == "v1") {
            validateV1();
        } else if (versionStr == "v2") {
            validateV2();
        } else {
            throwErrorException("Unknown CALQ version");
        }
        validateCompress();
    }
}

// -----------------------------------------------------------------------------

void ProgramOptions::processCommandLine(
        int argc,
        char *argv[]
){
    namespace po = boost::program_options;

    // Declare the supported options
    po::options_description options("Options");

    options.add_options()
            (
                    "help,h",
                    "[D, V1, V2] Print help"
            )
            (
                    "force,f",
                    po::bool_switch(&(this->force))->default_value(false),
                    "[D, V1, V2] Force overwriting of output files\n"
            )
            (
                    "input_file_path,i",
                    po::value<std::string>(&(this->inputFilePath))->required(),
                    "[D, V1, V2] Input file path"
            )
            (
                    "output_file_path,o",
                    po::value<std::string>(&(this->outputFilePath))->required(),
                    "[D, V1, V2] Output file path"
            )
            (
                    "pileup",
                    po::bool_switch(&(this->options.debugPileup)),
                    "[-, V1, V2] Be verbose and print pileup"
            )
            (
                    "streams",
                    po::bool_switch(&(this->debugStreams))
                            ->default_value(false),
                    "[-, V1, V2] Be verbose and print raw streams\n"
            )
            (
                    "blocksize,b",
                    po::value<size_t>(&(this->blockSize))->default_value(10000),
                    "[-, V1, V2] Block size (in number of SAM records)\n"
            )
            (
                    "quantization_min",
                    po::value<size_t>(&(this->quantizationMin))->
                            default_value(2),
                    "[-, V1, V2] Minimum quantization steps"
            )
            (
                    "quantization_max",
                    po::value<size_t>(&(this->quantizationMax))->
                            default_value(8),
                    "[-, V1, V2] Maximum quantization steps"
            )
            (
                    "polyploidy,p",
                    po::value<size_t>(&(this->polyploidy))->
                            default_value(2),
                    "[-, V1, V2] Polyploidy\n"
            )
            (
                    "qual_type,q",
                    po::value<std::string>(&(this->qualityValueType))->
                            default_value("Illumina-1.8+"),
                    "[-, V1, V2] Quality value type \nSanger: Phred+33 "
                    "[0,40];\nIllumina-1.3+: Phred+64 [0,40];\nIllumina-1.5+:"
                    " Phred+64 [0,40];\n"
                    "Illumina-1.8+: Phred+33 [0,41];\nMax33: Phred+33 [0,93];\n"
                    "Max64: Phred+64 [0,62];\n\n"
            )
            (
                    "quantizer_type",
                    po::value<std::string>(&(this->quantizerTypeStr))->
                            default_value("Uniform"),
                    "[-, V1, V2] Quantizer type\n(Uniform; Lloyd)\n"
            )
            (
                    "calq_version",
                    po::value<std::string>(&(this->versionStr))->
                            default_value("v1"),
                    "[-, V1, V2] v1 or v2"
            )
            (
                    "decompress,d",
                    po::bool_switch(&(this->decompress))->default_value(false),
                    "[D, --, --] Decompress"
            )
            (
                    "side_info_file_path,s",
                    po::value<std::string>(&(this->sideInformationFilePath)),
                    "[D, --, --] Side information file path"
            )
            (
                    "filter_size",
                    po::value<size_t>(&(this->options.filterSize))->
                            default_value(17),
                    "[-, --, V2] Haplotyper filter radius"
            )
            (
                    "filter_type",
                    po::value<std::string>(&(this->filterTypeStr))->
                            default_value("Gauss"),
                    "[-, --, V2] Haplotyper Filter Type (Gauss; Rectangle)\n"
            )
            (
                    "filter_cutoff",
                    po::value<size_t>(&(this->options.filterCutOff))->
                            default_value(50),
                    "[-, --, V2] Haplotyper filter cutoff radius\n"
            )
            (
                    "hq_softclip_threshold",
                    po::value<size_t>(&(this->hqSoftClipThreshold))->
                            default_value(29),
                    "[-, --, V2] Quality (no offset) at which a "
                    "softclip is considered HQ\n"
            )
            (
                    "hq_softclip_streak",
                    po::value<size_t>(&(this->options.hqSoftClipStreak))->
                            default_value(7),
                    "[-, --, V2] Number of hq-softclips in one "
                    "position which triggers spreading of activity values\n"
            )
            (
                    "hq_softclip_propagation",
                    po::value<size_t>(&(this->options.hqSoftClipPropagation))->
                            default_value(50),
                    "[-, --, V2] Distance at which hq-softclips impact "
                    "the activity score\n"
            )
            (
                    "reference_file_path,r",
                    po::value<std::string>(&(this->referenceFilePath)),
                    "[-, --, V2] Reference file (FASTA format)\n"
            )
            (
                    "no_squash",
                    po::bool_switch(&(this->options.squash))
                            ->default_value(false),
                    "[-, --, V2] Don't squash activity values "
                    "between 0.0 and 1.0"
            );

    // Parse the command line
    po::variables_map optionsMap;
    po::store(po::parse_command_line(argc, argv, options), optionsMap);

    // First thing to do is to print the help
    if (optionsMap.count("help") || optionsMap.count("h")) {
        std::cout << "Usage calq V1:\n $ ./calqapp -i [input_file] "
                     "-o [output_file] [optional parameters]\n\n"
                     "Usage calq V2:\n $ ./calqapp -i [input_file] "
                     "-o [output_file] -r [reference file] "
                     "--calq_version v2 [optional parameters]\n\n"
                     "Usage calq Decoder:\n $ ./calqapp -i [input_file] "
                     "-o [output_file] -s [sam side file] "
                     "-d [optional parameters]\n" << std::endl;
        std::cout << options;
        help = true;
        return;
    }

    // po::notify() will // throw on erroneous program options, that's why we
    // call it after printing the help
    po::notify(optionsMap);
}

// -----------------------------------------------------------------------------

}  // namespace calqapp

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

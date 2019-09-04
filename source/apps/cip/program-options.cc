#include "program-options.h"

#include <sstream>
#include "errors.h"
#include "util/helpers.h"
// #include "logging.h"
#include <cli11/cli11.h>

namespace cip {

ProgramOptions::ProgramOptions(int argc, char *argv[])
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
      sideInformationFilePath() {
    processCommandLine(argc, argv);
}

void ProgramOptions::validateCompress() {
    //     if (cip::fileNameExtension(inputFilePath) != std::string("sam")) {
    //         throwErrorException("Input file name extension must be 'sam'");
    //     }
    //
    //     CALQ_LOG("Block size: %d", static_cast<int>(blockSize));
    //     if (blockSize < 1) {
    //         throwErrorException("Block size must be greater than 0");
    //     }
    //
    //     CALQ_LOG("Quantization min steps: %d",
    //              static_cast<int>(options.quantizationMin));
    //     if (options.quantizationMin < 2) {
    //         throwErrorException("Quantization must be greater than 1");
    //     }
    //
    //     CALQ_LOG("Quantization max steps: %d",
    //              static_cast<int>(options.quantizationMax));
    //     if (options.quantizationMax < 2 ||
    //         options.quantizationMax < options.quantizationMin) {
    //         throwErrorException(
    //             "Quantization must be greater than"
    //             " 1 and quantizationMin");
    //     }
    //
    //     CALQ_LOG("Polyploidy: %d", static_cast<int>(options.polyploidy));
    //     if (options.polyploidy < 1) {
    //         throwErrorException("Polyploidy must be greater than 0");
    //     }
    //
    //     CALQ_LOG("Quantizer type: %s", quantizerTypeStr.c_str());
    //     if (quantizerTypeStr == "Uniform") {
    //         options.quantizerType = calq::QuantizerType::UNIFORM;
    //     } else if (quantizerTypeStr == "Lloyd") {
    //         options.quantizerType = calq::QuantizerType::LLOYD_MAX;
    //     } else {
    //         throwErrorException("Quantizer type not supported");
    //     }
    //
    //     CALQ_LOG("Quality value type: %s", qualityValueType.c_str());
    //     if (qualityValueType == "Sanger") {
    //         // Sanger: Phred+33 [0,40]
    //         options.qualityValueOffset = 33;
    //         options.qualityValueMin = 0;
    //         options.qualityValueMax = 40;
    //     } else if (qualityValueType == "Illumina-1.3+") {
    //         // Illumina 1.3+: Phred+64 [0,40]
    //         options.qualityValueOffset = 64;
    //         options.qualityValueMin = 0;
    //         options.qualityValueMax = 40;
    //     } else if (qualityValueType == "Illumina-1.5+") {
    //         // Illumina 1.5+: Phred+64 [0,40] with 0=unused,
    //         // 1=unused, 2=Read Segment Quality Control Indicator ('B')
    //         CALQ_LOG(
    //             "Warning: Read Segment Quality Control Indicator"
    //             " will not be treated specifically by CALQ");
    //         options.qualityValueOffset = 64;
    //         options.qualityValueMin = 0;
    //         options.qualityValueMax = 40;
    //     } else if (qualityValueType == "Illumina-1.8+") {
    //         // Illumina 1.8+ Phred+33 [0,41]
    //         options.qualityValueOffset = 33;
    //         options.qualityValueMin = 0;
    //         options.qualityValueMax = 41;
    //     } else if (qualityValueType == "Max33") {
    //         // Max33 Phred+33 [0,93]
    //         options.qualityValueOffset = 33;
    //         options.qualityValueMin = 0;
    //         options.qualityValueMax = 93;
    //     } else if (qualityValueType == "Max64") {
    //         // Max64 Phred+64 [0,62]
    //         options.qualityValueOffset = 64;
    //         options.qualityValueMin = 0;
    //         options.qualityValueMax = 62;
    //     } else {
    //         throwErrorException("Quality value type not supported");
    //     }
    //     CALQ_LOG("Quality value offset: %d",
    //              static_cast<int>(options.qualityValueOffset));
    //     CALQ_LOG("Quality value range: [%d,%d]",
    //              static_cast<int>(options.qualityValueMin),
    //              static_cast<int>(options.qualityValueMax));
    //
    //     if (options.debugPileup) {
    //         CALQ_LOG("Outputting full pileup as debug information");
    //     }
    //
    //     if (debugStreams) {
    //         CALQ_LOG("Outputting full streams as debug information");
    //     }
}
//
// void ProgramOptions::validateV2() {
//     options.version = calq::Version::V2;
//     CALQ_LOG("Compressing using CALQ version 2");
//
//     if (options.squash) {
//         CALQ_LOG("Acitivity scores are squashed between 0.0 and 1.0");
//     } else {
//         CALQ_LOG("Acitivity scores are !NOT! squashed between 0.0 and 1.0");
//     }
//
//     CALQ_LOG("Filter size: %d", static_cast<int>(options.filterSize));
//     if (options.filterSize < 1) {
//         throwErrorException("Filter size must be greater than 0");
//     }
//
//     CALQ_LOG("hqSoftClipPropagation: %d",
//              static_cast<int>(options.hqSoftClipPropagation));
//
//     CALQ_LOG("hqSoftClipStreak: %d",
//              static_cast<int>(options.hqSoftClipStreak));
//
//     CALQ_LOG("hqSoftClipThreshold: %d",
//              static_cast<int>(options.hqSoftClipThreshold));
//
//     CALQ_LOG("filterCutOff: %d", static_cast<int>(options.filterCutOff));
//
//     CALQ_LOG("Filter type: %s", filterTypeStr.c_str());
//     if (filterTypeStr == "Gauss") {
//         options.filterType = calq::FilterType::GAUSS;
//     } else if (filterTypeStr == "Rectangle") {
//         options.filterType = calq::FilterType::RECTANGLE;
//     } else {
//         throwErrorException("Filter type not supported");
//     }
//
//     CALQ_LOG("  %s", referenceFilePath.c_str());
//     if (referenceFilePath.empty()) {
//         throwErrorException("Reference file name not provided");
//     }
//     if (!cip::fileExists(referenceFilePath)) {
//         throwErrorException("Cannot access reference file");
//     }
//     if (cip::fileNameExtension(referenceFilePath) != std::string("fa") &&
//         cip::fileNameExtension(referenceFilePath) != std::string("fasta")) {
//         throwErrorException(
//             "Reference file name extension must "
//             "be 'fa' or 'fasta'");
//     }
// }

void ProgramOptions::validateDecompress() {
    if (util::fileNameExtension(sideInformationFilePath) != std::string("sam")) {
        throwErrorException("Side information file name extension must be '.sam'");
    }
    if (!util::fileExists(sideInformationFilePath)) {
        throwErrorException("Cannot access side information file");
    }
}

void ProgramOptions::validateCommon() {
    //     this->options.squash = !this->options.squash;
    //     if (this->quantizationMin > 255) {
    //         throwErrorException("Option quantizationMin too big");
    //     }
    //     this->options.quantizationMin = uint8_t(this->quantizationMin);
    //
    //     if (this->quantizationMax > 255) {
    //         throwErrorException("Option quantizationMax too big");
    //     }
    //     this->options.quantizationMax = uint8_t(this->quantizationMax);
    //
    //     if (this->polyploidy > 255) {
    //         throwErrorException("Option polyploidy too big");
    //     }
    //     this->options.polyploidy = uint8_t(this->polyploidy);
    //
    //     if (this->hqSoftClipThreshold > 255) {
    //         throwErrorException("Option hqSoftClipThreshold too big");
    //     }
    //     this->options.hqSoftClipThreshold = uint8_t(this->hqSoftClipThreshold);

    // force
    if (force) {
        std::cout << "Force switch set - overwriting output file(s)" << std::endl;
    }

    // inputFilePath
    if (inputFilePath.empty()) {
        throwErrorException("No input file name provided");
    }
    if (!util::fileExists(inputFilePath)) {
        throwErrorException("Cannot access input file");
    }

    // polyploidy
    if (polyploidy < 1 || polyploidy > 6) {
    }

    // outputFilePath
    if (util::fileExists(outputFilePath) && !force) {
        throwErrorException("Not overwriting output file (use option '-f' to force overwriting)");
    }
}

void ProgramOptions::validate() {
    validateCommon();

    if (decompress) {
        validateDecompress();
    } else {
        validateCompress();
    }
}

void ProgramOptions::processCommandLine(int argc, char *argv[]) {
    CLI::App app("CALQ");

    std::stringstream qualityValueTypeDesc;
    qualityValueTypeDesc << "Quality value type" << std::endl;
    qualityValueTypeDesc << "  Sanger:        Phred+33 [0,40]" << std::endl;
    qualityValueTypeDesc << "  Illumina-1.3+: Phred+64 [0,40]" << std::endl;
    qualityValueTypeDesc << "  Illumina-1.5+: Phred+64 [0,40]" << std::endl;
    qualityValueTypeDesc << "  Illumina-1.8+: Phred+33 [0,41]" << std::endl;
    qualityValueTypeDesc << "  Max33:         Phred+33 [0,93]" << std::endl;
    qualityValueTypeDesc << "  Max64:         Phred+64 [0,62]" << std::endl;

    std::stringstream quantizerTypeDesc;
    quantizerTypeDesc << "Quantizer type" << std::endl;
    quantizerTypeDesc << "  Uniform" << std::endl;
    quantizerTypeDesc << "  UniformMinMax" << std::endl;
    quantizerTypeDesc << "  MaxLloyd" << std::endl;

    app.add_option("-b,--block-size", blockSize, "Block size (in number of SAM records)")->default_val("10000");
    app.add_flag("-d,--decompress", decompress, "Decompress");
    app.add_flag("-f,--force", force, "Force overwriting of output file(s)");
    app.add_option("-i,--input-file", inputFilePath, "Input file")->required();
    app.add_option("--max-q-steps", manNumQuantSteps, "Maximum number of quantization steps")->default_val("8");
    app.add_option("--min-q-steps", minNumQuantSteps, "Minimum number of quantization steps")->default_val("2");
    app.add_option("-o,--output-file", outputFilePath, "Output file")->required();
    app.add_option("-p,--polyploidy", polyploidy, "Polyploidy")->default_val("2");
    app.add_option("--qual-type", qualityValueType, qualityValueTypeDesc.str())->default_str("Illumina-1.8+");
    // app.add_option("--quant-type", quantizerType, quantizerTypeDesc.str())->default_str("Uniform");

    //     options.add_options()(
    //         "calq_version",
    //         po::value<std::string>(&(this->versionStr))->default_value("v1"),
    //         "[-, V1, V2] v1 or v2")(
    //         "side_info_file_path,s",
    //         po::value<std::string>(&(this->sideInformationFilePath)),
    //         "[D, --, --] Side information file path")(
    //         "filter_size",
    //         po::value<size_t>(&(this->options.filterSize))->default_value(17),
    //         "[-, --, V2] Haplotyper filter radius")(
    //         "filter_type",
    //         po::value<std::string>(&(this->filterTypeStr))->default_value("Gauss"),
    //         "[-, --, V2] Haplotyper Filter Type (Gauss; Rectangle)\n")(
    //         "filter_cutoff",
    //         po::value<size_t>(&(this->options.filterCutOff))->default_value(50),
    //         "[-, --, V2] Haplotyper filter cutoff radius\n")(
    //         "hq_softclip_threshold",
    //         po::value<size_t>(&(this->hqSoftClipThreshold))->default_value(29),
    //         "[-, --, V2] Quality (no offset) at which a "
    //         "softclip is considered HQ\n")(
    //         "hq_softclip_streak",
    //         po::value<size_t>(&(this->options.hqSoftClipStreak))->default_value(7),
    //         "[-, --, V2] Number of hq-softclips in one "
    //         "position which triggers spreading of activity values\n")(
    //         "hq_softclip_propagation",
    //         po::value<size_t>(&(this->options.hqSoftClipPropagation))
    //             ->default_value(50),
    //         "[-, --, V2] Distance at which hq-softclips impact "
    //         "the activity score\n")(
    //         "reference_file_path,r",
    //         po::value<std::string>(&(this->referenceFilePath)),
    //         "[-, --, V2] Reference file (FASTA format)\n")(
    //         "no_squash",
    //         po::bool_switch(&(this->options.squash))->default_value(false),
    //         "[-, --, V2] Don't squash activity values "
    //         "between 0.0 and 1.0");

    try {
        app.parse(argc, argv);
    } catch (const CLI::ParseError &e) {
        throwErrorException("Command line parsing failed: " + std::to_string(app.exit(e)));
    }

    validate();
}

}  // namespace cip

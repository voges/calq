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
        verbose(),
        test(false),
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

void ProgramOptions::validate(void){
    if (versionStr == "v1")
    {
        options.version = calq::Version::V1;
        CALQ_LOG("Using CALQ version 1");
    }
    else if (versionStr == "v2")
    {
        options.version = calq::Version::V2;
        CALQ_LOG("Using CALQ version 2");
    }
    else
    {
        throwErrorException("Unknown CALQ version");
    }

    // force
    if (force)
    {
        CALQ_LOG("Force switch set - overwriting output file(s)");
    }

    if (verbose)
    {
        CALQ_LOG("Debug switch set - verbose output");
    }

    if (test)
    {
        CALQ_LOG("Test switch set - running test cases instead of compression");
    }

    if (options.squash)
    {
        CALQ_LOG("Acitivity scores are squashed between 0.0 and 1.0");
    }
    else
    {
        CALQ_LOG("Acitivity scores are !NOT! squashed between 0.0 and 1.0");
    }

    // inputFileName
    CALQ_LOG("Input file name: %s", inputFilePath.c_str());
    if (inputFilePath.empty())
    {
        throwErrorException("No input file name provided");
    }
    if (!decompress)
    {
        if (calqapp::fileNameExtension(inputFilePath) != std::string("sam"))
        {
            throwErrorException("Input file name extension must be 'sam'");
        }
    }
    else
    {
        if (calqapp::fileNameExtension(inputFilePath) != std::string("cq"))
        {
            CALQ_LOG("Warning: Input file name extension is not 'cq'");
//             throwErrorException("Input file name extension must be 'cq'");
        }
    }
    if (!calqapp::fileExists(inputFilePath))
    {
        throwErrorException("Cannot access input file");
    }

    // outputFileName
    if (!decompress)
    {
        if (outputFilePath.empty())
        {
            CALQ_LOG("No output file name provided - constructing output"
                     " file name from input file name"
            );
            outputFilePath += inputFilePath + ".cq";
        }
    }
    else
    {
        if (outputFilePath.empty())
        {
            CALQ_LOG("No output file name provided - constructing output"
                     " file name from input file name"
            );
            outputFilePath += inputFilePath + ".qual";
        }
    }
    CALQ_LOG("Output file name: %s", outputFilePath.c_str());
    if (calqapp::fileExists(outputFilePath))
    {
        if (!force)
        {
            throwErrorException(
                    "Not overwriting output file "
                    "(use option 'f' to force overwriting)"
            );
        }
    }

    // blockSize
    if (!decompress)
    {
        CALQ_LOG("Block size: %d", static_cast<int>(blockSize));
        if (blockSize < 1)
        {
            throwErrorException("Block size must be greater than 0");
        }
    }

    // Haplotyper filter size
    if (!decompress)
    {
        CALQ_LOG("Filter size: %d", static_cast<int>(options.filterSize));
        if (options.filterSize < 1)
        {
            throwErrorException("Filter size must be greater than 0");
        }
    }

    // Quantization
    if (!decompress)
    {
        CALQ_LOG("Quantization min steps: %d",
                 static_cast<int>(options.quantizationMin));
        if (options.quantizationMin < 2)
        {
            throwErrorException("Quantization must be greater than 1");
        }
    }

    // Quantization
    if (!decompress)
    {
        CALQ_LOG("Quantization max steps: %d",
                 static_cast<int>(options.quantizationMax));
        if (options.quantizationMax < 2
            || options.quantizationMax < options.quantizationMin)
        {
            throwErrorException(
                    "Quantization must be greater than"
                    " 1 and quantizationMin"
            );
        }
    }

    // polyploidy
    if (!decompress)
    {
        CALQ_LOG("Polyploidy: %d", static_cast<int>(options.polyploidy));
        if (options.polyploidy < 1)
        {
            throwErrorException("Polyploidy must be greater than 0");
        }
    }

    // qualityValueType
    if (!decompress)
    {
        CALQ_LOG("Filter type: %s", filterTypeStr.c_str());
        if (filterTypeStr == "Gauss")
        {
            options.filterType = calq::FilterType::GAUSS;
        }
        else if (filterTypeStr == "Rectangle")
        {
            options.filterType = calq::FilterType::RECTANGLE;
        }
        else
        {
            throwErrorException("Filter type not supported");
        }
    }

    // qualityValueType
    if (!decompress)
    {
        CALQ_LOG("Quantizer type: %s", quantizerTypeStr.c_str());
        if (quantizerTypeStr == "Uniform")
        {
            options.quantizerType =
                    calq::QuantizerType::UNIFORM;
        }
        else if (quantizerTypeStr == "Lloyd")
        {
            options.quantizerType =
                    calq::QuantizerType::LLOYD_MAX;
        }
        else
        {
            throwErrorException("Quantizer type not supported");
        }
    }

    // qualityValueType
    if (!decompress)
    {
        CALQ_LOG("Quality value type: %s", qualityValueType.c_str());
        if (qualityValueType == "Sanger")
        {
            // Sanger: Phred+33 [0,40]
            options.qualityValueOffset = 33;
            options.qualityValueMin = 0;
            options.qualityValueMax = 40;
        }
        else if (qualityValueType == "Illumina-1.3+")
        {
            // Illumina 1.3+: Phred+64 [0,40]
            options.qualityValueOffset = 64;
            options.qualityValueMin = 0;
            options.qualityValueMax = 40;
        }
        else if (qualityValueType == "Illumina-1.5+")
        {
            // Illumina 1.5+: Phred+64 [0,40] with 0=unused,
            // 1=unused, 2=Read Segment Quality Control Indicator ('B')
            CALQ_LOG("Warning: Read Segment Quality Control Indicator"
                     " will not be treated specifically by CALQ"
            );
            options.qualityValueOffset = 64;
            options.qualityValueMin = 0;
            options.qualityValueMax = 40;
        }
        else if (qualityValueType == "Illumina-1.8+")
        {
            // Illumina 1.8+ Phred+33 [0,41]
            options.qualityValueOffset = 33;
            options.qualityValueMin = 0;
            options.qualityValueMax = 41;
        }
        else if (qualityValueType == "Max33")
        {
            // Max33 Phred+33 [0,93]
            options.qualityValueOffset = 33;
            options.qualityValueMin = 0;
            options.qualityValueMax = 93;
        }
        else if (qualityValueType == "Max64")
        {
            // Max64 Phred+64 [0,62]
            options.qualityValueOffset = 64;
            options.qualityValueMin = 0;
            options.qualityValueMax = 62;
        }
        else
        {
            throwErrorException("Quality value type not supported");
        }
        CALQ_LOG("Quality value offset: %d",
                 static_cast<int>(options.qualityValueOffset)
        );
        CALQ_LOG("Quality value range: [%d,%d]",
                 static_cast<int>(options.qualityValueMin),
                 static_cast<int>(options.qualityValueMax));

        // referenceFiles
        if (!decompress)
        {
            if (referenceFilePath.empty())
            {
                CALQ_LOG("Operating without reference file");
            }
            else
            {
                CALQ_LOG("Operating with reference file:");
                CALQ_LOG("  %s", referenceFilePath.c_str());
                if (referenceFilePath.empty())
                {
                    throwErrorException("Reference file name not provided");
                }
                if (!calqapp::fileExists(referenceFilePath))
                {
                    throwErrorException("Cannot access reference file");
                }
                if (calqapp::fileNameExtension(referenceFilePath)
                    != std::string("fa")
                    && calqapp::fileNameExtension(referenceFilePath)
                       != std::string("fasta"))
                {
                    throwErrorException(
                            "Reference file name extension must "
                            "be 'fa' or 'fasta'"
                    );
                }
            }
        }
    }

    // decompress
    if (!decompress)
    {
        CALQ_LOG("Compressing");
    }
    else
    {
        CALQ_LOG("Decompressing");
    }

    // sideInformationFilePath
    if (decompress)
    {
        CALQ_LOG("Side information file name: %s",
                 sideInformationFilePath.c_str()
        );
        if (sideInformationFilePath.empty())
        {
            throwErrorException("No side information file name provided");
        }
        if (calqapp::fileNameExtension(sideInformationFilePath)
            != std::string("sam"))
        {
            throwErrorException(
                    "Side information file name "
                    "extension must be 'sam'"
            );
        }
        if (!calqapp::fileExists(sideInformationFilePath))
        {
            throwErrorException("Cannot access side information file");
        }
    }
}

// -----------------------------------------------------------------------------

static void writeUsage(std::ostream&){

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
                    "Print help"
            )
            (
                    "force,f",
                    po::bool_switch(&(this->force))->default_value(false),
                    "Force overwriting of output files."
            )
            (
                    "verbose,v",
                    po::bool_switch(&(this->verbose))->default_value(false),
                    "Be verbose."
            )
            (
                    "input_file_path,i",
                    po::value<std::string>(&(this->inputFilePath))->required(),
                    "Input file path"
            )
            (
                    "output_file_path,o",
                    po::value<std::string>(&(this->outputFilePath))->required(),
                    "Output file path"
            )
            (
                    "blocksize,b",
                    po::value<size_t>(&(this->blockSize))->default_value(10000),
                    "Block size (in number of SAM records). Default 10000."
            )
            (
                    "filtersize",
                    po::value<size_t>(&(this->options.filterSize))->
                            default_value(17),
                    "Haplotyper filter radius. Default 17. (v2 only)"
            )
            (
                    "quantization_min",
                    po::value<uint8_t>(&(this->options.quantizationMin))->
                            default_value(2),
                    "Minimum quantization steps. Default 2."
            )
            (
                    "quantization_max",
                    po::value<uint8_t>(&(this->options.quantizationMax))->
                            default_value(8),
                    "Maximum quantization steps. Default 8."
            )
            (
                    "polyploidy,p",
                    po::value<uint8_t>(&(this->options.polyploidy))->
                            default_value(2),
                    "Polyploidy. Default 2."
            )
            (
                    "quality_value_type,q",
                    po::value<std::string>(&(this->qualityValueType))->
                            default_value("Illumina-1.8+"),
                    "Quality value type (Sanger: Phred+33 [0,40]; "
                    "Illumina-1.3+: Phred+64 [0,40]; Illumina-1.5+:"
                    " Phred+64 [0,40]; "
                    "Illumina-1.8+: Phred+33 [0,41]; Max33: Phred+33 [0,93];"
                    " Max64: Phred+64 [0,62])"
            )
            (
                    "reference_file_path,r",
                    po::value<std::string>(&(this->referenceFilePath)),
                    "Reference file (FASTA format) (v2 only)"
            )
            (
                    "filter_type",
                    po::value<std::string>(&(this->filterTypeStr))->
                            default_value("Gauss"),
                    "Haplotyper Filter Type (Gauss; Rectangle). "
                    "Default Gauss. (v2 only)"
            )
            (
                    "quantizer_type",
                    po::value<std::string>(&(this->quantizerTypeStr))->
                            default_value("Uniform"),
                    "Quantizer type (Uniform; Lloyd). Default Uniform."
            )
            (
                    "CALQ-version",
                    po::value<std::string>(&(this->versionStr))->
                            default_value("v1"),
                    "v1 or v2. Default v1"
            )
            (
                    "decompress,d",
                    po::bool_switch(&(this->decompress))->default_value(false),
                    "Decompress."
            )
            (
                    "side_information_file_path,s",
                    po::value<std::string>(&(this->sideInformationFilePath)),
                    "Side information file path"
            );

    this->options.debug = false;
    this->options.squash = true;



    // Parse the command line
    po::variables_map optionsMap;
    po::store(po::parse_command_line(argc, argv, options), optionsMap);

    // First thing to do is to print the help
    if (optionsMap.count("help") || optionsMap.count("h"))
    {
        writeUsage(std::cout);
        std::cout << options;
        // TO-DO: Change this to Error instead of log.
        // Program should terminate here.
        CALQ_LOG("aborting program after printing help");
    }

    // po::notify() will // throw on erroneous program options, that's why we
    // call it after printing the help
    po::notify(optionsMap);

    // Yes, we are verbose
    if (this->verbose)
    {
        std::cout << "calqapp: yes, we are verbose" << std::endl;
    }
}

// -----------------------------------------------------------------------------

}  // namespace calqapp

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
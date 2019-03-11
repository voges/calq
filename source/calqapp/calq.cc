#include <chrono>
#include <numeric>

// -----------------------------------------------------------------------------

#include <boost/program_options.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

// -----------------------------------------------------------------------------

#include "calq/calq_coder.h"
#include "calq/exceptions.h"

// -----------------------------------------------------------------------------

#include "calqapp/cq_file.h"
#include "calqapp/fasta_file.h"
#include "calqapp/logging.h"
#include "calqapp/program_options.h"
#include "calqapp/SAMFileHandler.h"

// -----------------------------------------------------------------------------

#include "gabac/configuration.h"

// -----------------------------------------------------------------------------

void json2config(const std::string& json, gabac::EncodingConfiguration *enc) {
    try
    {
        // Read the stringstream JSON data to a property tree
        std::stringstream tmp(json);
        boost::property_tree::ptree propertyTree;
        boost::property_tree::read_json(tmp, propertyTree);

        // Convert the property tree contents to our internal structure
        enc->wordSize
                = propertyTree.get<unsigned int>("word_size");
        enc->sequenceTransformationId
                =
                static_cast<gabac::SequenceTransformationId>(propertyTree.get<unsigned int>("sequence_transformation_id"));
        enc->sequenceTransformationParameter
                = propertyTree.get<unsigned int>("sequence_transformation_parameter");
        for (const auto& child : propertyTree.get_child("transformed_sequences")) {
            // Declare a transformed sequence configuration
            gabac::TransformedSequenceConfiguration transformedSequenceConfiguration;

            // Fill the transformed sequence configuration
            transformedSequenceConfiguration.lutTransformationEnabled
                    = static_cast<bool>(child.second.get<unsigned int>("lut_transformation_enabled"));
            transformedSequenceConfiguration.lutBits = enc->wordSize*8;
            transformedSequenceConfiguration.lutOrder = 0;
            if (transformedSequenceConfiguration.lutTransformationEnabled) {
                if(child.second.count("lut_transformation_bits") > 0) {
                    transformedSequenceConfiguration.lutBits
                            = child.second.get<unsigned int>("lut_transformation_bits");
                }
                if(child.second.count("lut_transformation_order") > 0) {
                    transformedSequenceConfiguration.lutOrder
                            = child.second.get<unsigned int>("lut_transformation_order");
                }
            } else {
                transformedSequenceConfiguration.lutBits = 0;
                transformedSequenceConfiguration.lutOrder = 0;
            }
            transformedSequenceConfiguration.diffCodingEnabled
                    = child.second.get<bool>("diff_coding_enabled");
            transformedSequenceConfiguration.binarizationId
                    = static_cast<gabac::BinarizationId>(child.second.get<unsigned int>("binarization_id"));
            for (const auto& grandchild : child.second.get_child("binarization_parameters"))
            {
                transformedSequenceConfiguration.binarizationParameters
                        .push_back(grandchild.second.get_value<unsigned int>());
            }
            transformedSequenceConfiguration.contextSelectionId
                    = static_cast<gabac::ContextSelectionId>(child.second.get<unsigned int>("context_selection_id"));

            // Append the filled transformed sequence configuration to our
            // list of transformed sequence configurations
            enc->transformedSequenceConfigurations.push_back(transformedSequenceConfiguration);
        }
    }
    catch (const boost::property_tree::ptree_error& e)
    {
        throwErrorException("Boost ptree failed.");  
    }
} 

// -----------------------------------------------------------------------------


static int encode(const calqapp::ProgramOptions& programOptions){
    auto startTime = std::chrono::steady_clock::now();

    size_t compressedMappedQualSize = 0;
    size_t compressedUnmappedQualSize = 0;
    size_t uncompressedMappedQualSize = 0;
    size_t uncompressedUnmappedQualSize = 0;

    calqapp::SAMFileHandler sH(
            programOptions.inputFilePath,
            programOptions.referenceFilePath
    );

    calqapp::CQFile file(
            programOptions.outputFilePath,
            calqapp::File::Mode::MODE_WRITE
    );
    file.writeHeader(programOptions.blockSize);
        
    gabac::EncodingConfiguration config;
    calqapp::File configurationFile(
            "config.json",
            calqapp::File::Mode::MODE_READ
            );
    std::string jsonInput("\0", configurationFile.size());
    configurationFile.read(&jsonInput[0], jsonInput.size());
    json2config(jsonInput, &config);


    while (sH.readBlock(programOptions.blockSize) != 0) {
        calq::SideInformation encSide;
        sH.getSideInformation(&encSide);

        calqapp::UnmappedInformation unmappedInfo;
        sH.getUnmappedBlock(&unmappedInfo);

        calq::EncodingBlock encBlock;
        sH.getMappedBlock(&encBlock);

        calq::DecodingBlock decBlock;


        calq::encode(
                programOptions.options,
                encSide,
                encBlock,
                &decBlock
        );

        // Build single string out of unmapped q-values
        std::string unmappedString;
        unmappedString = std::accumulate(
                unmappedInfo.unmappedQualityScores.begin(),
                unmappedInfo.unmappedQualityScores.end(),
                unmappedString
        );

        size_t mappedSize;
        size_t unmappedSize;

        file.writeBlock(
                programOptions.options,
                config,
                decBlock,
                encSide,
                unmappedString,
                programOptions.debugStreams,
                &mappedSize,
                &unmappedSize
        );

        compressedMappedQualSize += mappedSize;
        compressedUnmappedQualSize += unmappedSize;

        uncompressedMappedQualSize = std::accumulate(
                encBlock.qvalues.begin(),
                encBlock.qvalues.end(),
                uncompressedMappedQualSize,
                [](const size_t& a, const std::string& b)->size_t{
                    return a + b.size();
                }
        );
        uncompressedUnmappedQualSize += unmappedString.length();
    }

    file.close();

    auto stopTime = std::chrono::steady_clock::now();
    auto diffTime = stopTime - startTime;
    auto diffTimeMs = std::chrono::duration_cast<std::chrono::milliseconds>(
            diffTime
    ).count();
    auto diffTimeS = std::chrono::duration_cast<std::chrono::seconds>(
            diffTime
    ).count();
    auto diffTimeM = std::chrono::duration_cast<std::chrono::minutes>(
            diffTime
    ).count();
    auto diffTimeH = std::chrono::duration_cast<std::chrono::hours>(
            diffTime
    ).count();

    CALQ_LOG("COMPRESSION STATISTICS");
    CALQ_LOG("  Took %d ms ~= %d s ~= %d m ~= %d h",
             (int) diffTimeMs,
             (int) diffTimeS,
             (int) diffTimeM,
             (int) diffTimeH);
    const unsigned MB = 1 * 1000 * 1000;
    CALQ_LOG("  Speed (uncompressed size/time): %.2f MB/s",
             ((static_cast<double>(uncompressedMappedQualSize
                                   + uncompressedUnmappedQualSize) /
               static_cast<double>(MB))) / (static_cast<double>(diffTimeS)));
    CALQ_LOG("  Wrote %zu block(s)", sH.nrBlocksRead());
    CALQ_LOG("  Record(s):  %12zu", sH.nrRecordsRead());
    CALQ_LOG("    Mapped:   %12zu", sH.nrMappedRecordsRead());
    CALQ_LOG("    Unmapped: %12zu", sH.nrUnmappedRecordsRead());
    CALQ_LOG("  Uncompressed size: %12zu",
             uncompressedMappedQualSize + uncompressedUnmappedQualSize);
    CALQ_LOG("    Mapped:          %12zu", uncompressedMappedQualSize);
    CALQ_LOG("    Unmapped:        %12zu", uncompressedUnmappedQualSize);
    CALQ_LOG("  Compressed size: %12zu", file.nrWrittenBytes());
    CALQ_LOG("    File format:   %12zu", file.nrWrittenFileFormatBytes());
    CALQ_LOG("    Mapped:        %12zu", compressedMappedQualSize);
    CALQ_LOG("    Unmapped:      %12zu", compressedUnmappedQualSize);
    CALQ_LOG("  Compression ratio: %4.2f%%",
             (double) file.nrWrittenBytes() *
             100 /
             (double) (uncompressedMappedQualSize
                       + uncompressedUnmappedQualSize));
    CALQ_LOG("    Mapped:          %4.2f%%",
             (double) compressedMappedQualSize * 100
             / (double) (uncompressedMappedQualSize));
    CALQ_LOG("    Unmapped:        %4.2f%%",
             (double) compressedUnmappedQualSize * 100
             / (double) (uncompressedUnmappedQualSize));
    CALQ_LOG("  Compression factor: %4.2f",
             (double) (uncompressedMappedQualSize
                       + uncompressedUnmappedQualSize) /
             (double) file.nrWrittenBytes());
    CALQ_LOG("    Mapped:           %4.2f",
             (double) (uncompressedMappedQualSize)
             / (double) compressedMappedQualSize);
    CALQ_LOG("    Unmapped:         %4.2f",
             (double) (uncompressedUnmappedQualSize)
             / (double) compressedUnmappedQualSize);
    CALQ_LOG("  Bits per quality value: %2.4f",
             (static_cast<double>(file.nrWrittenBytes()) * 8) /
             static_cast<double>(uncompressedMappedQualSize
                                 + uncompressedUnmappedQualSize));
    CALQ_LOG("    Mapped:               %2.4f",
             (static_cast<double>(compressedMappedQualSize) * 8) /
             static_cast<double>(uncompressedMappedQualSize));
    CALQ_LOG("    Unmapped:             %2.4f",
             (static_cast<double>(compressedUnmappedQualSize) * 8) /
             static_cast<double>(uncompressedUnmappedQualSize));

    return EXIT_SUCCESS;
}

// -----------------------------------------------------------------------------

static int decode(const calqapp::ProgramOptions& po){
    calqapp::ProgramOptions programOptions = po;

    auto startTime = std::chrono::steady_clock::now();

    calqapp::CQFile file(
            programOptions.inputFilePath,
            calqapp::File::Mode::MODE_READ
    );

    calqapp::File qualFile(
            programOptions.outputFilePath,
            calqapp::File::Mode::MODE_WRITE
    );

    file.readHeader(&programOptions.blockSize);

    calqapp::SAMFileHandler sH(
            programOptions.sideInformationFilePath,
            programOptions.referenceFilePath
    );

    gabac::EncodingConfiguration configuration;
    calqapp::File configurationFile(
            "config.json",
            calqapp::File::Mode::MODE_READ
            );
    std::string jsonInput("\0", configurationFile.size());
    configurationFile.read(&jsonInput[0], jsonInput.size());
    json2config(jsonInput, &configuration);

    while (sH.readBlock(programOptions.blockSize) != 0) {
        calq::DecodingBlock input;
        calq::EncodingBlock output;
        calq::SideInformation side;
        calqapp::UnmappedInformation unmappedInfo;

        sH.getSideInformation(&side);
        sH.getUnmappedBlock(&unmappedInfo);

        side.qualOffset = programOptions.options.qualityValueOffset;
        input.quantizerIndices.clear();

        std::string unmappedValues;
        calq::DecodingOptions opts;
        file.readBlock(&input, &side, &unmappedValues, configuration);
        calq::decode(opts, side, input, &output);


        auto mappedIt = output.qvalues.begin();
        auto unmappedPos = 0u;
        auto unmappedSideIt = unmappedInfo.unmappedQualityScores.begin();
        for (const auto& b : unmappedInfo.mappedFlags) {
            if (b) {
                qualFile.write(mappedIt->c_str(), mappedIt->length());
                qualFile.writeByte('\n');
                ++mappedIt;
            } else {
                std::string read = unmappedValues.substr(
                        unmappedPos,
                        unmappedSideIt->length()
                );
                qualFile.write(read.c_str(), read.length());
                qualFile.writeByte('\n');

                unmappedPos += read.length();
                ++unmappedSideIt;
            }
        }
    }

    auto stopTime = std::chrono::steady_clock::now();
    auto diffTime = stopTime - startTime;
    auto diffTimeMs =
            std::chrono::duration_cast<std::chrono::milliseconds>(diffTime)
                    .count();
    auto diffTimeS =
            std::chrono::duration_cast<std::chrono::seconds>(diffTime).count();
    auto diffTimeM =
            std::chrono::duration_cast<std::chrono::minutes>(diffTime).count();
    auto diffTimeH =
            std::chrono::duration_cast<std::chrono::hours>(diffTime).count();

    CALQ_LOG("DECOMPRESSION STATISTICS");
    CALQ_LOG("  Took %d ms ~= %d s ~= %d m ~= %d h",
             (int) diffTimeMs, (int) diffTimeS, (int) diffTimeM, (int) diffTimeH
    );
    const unsigned MB = 1 * 1000 * 1000;
    CALQ_LOG("  Speed (compressed size/time): %.2f MB/s",
             ((static_cast<double>(file.nrReadBytes()) /
               static_cast<double>(MB))) /
             (static_cast<double>(diffTimeS)));
    CALQ_LOG("  Decoded %zu block(s)", sH.nrBlocksRead());

    return EXIT_SUCCESS;
}

// -----------------------------------------------------------------------------

int main(int argc,
         char *argv[]
){
    try {
        calqapp::ProgramOptions programOptions(argc, argv);
        if (programOptions.help) {
            return EXIT_SUCCESS;
        }
        programOptions.validate();

        if (!programOptions.decompress) {
            return encode(programOptions);
        } else {
            return decode(programOptions);
        }
    }
    catch (const calq::ErrorException& errorException) {
        CALQ_ERROR("%s", errorException.what());
        return EXIT_FAILURE;
    }
    catch (const std::exception& stdException) {
        CALQ_ERROR("Fatal: %s", stdException.what());
        return EXIT_FAILURE;
    }
    catch (...) {
        CALQ_ERROR("Unknown error occurred");
        return EXIT_FAILURE;
    }
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

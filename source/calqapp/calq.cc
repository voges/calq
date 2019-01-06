#include <boost/program_options.hpp>

#include "calq/calq_encoder.h"
#include "calq/calq_decoder.h"
#include "calq/structs.h"
#include "calq/error_exception_reporter.h"
#include "calq/options.h"
#include "calqapp/program_options.h"
/*
#include "tclap/CmdLine.h"
*/
#include "calq/log.h"

#include "calq/cq_file.h"


#include "calqapp/SAMFileHandler.h"


size_t writeBlock(const calq::EncodingOptions& opts, const calq::DecodingBlock& block, const calq::EncodingSideInformation& side, const std::string& unmappedQualityValues_, calq::CQFile* cqFile) {

    // Write block parameters
    cqFile->writeUint32(side.positionStart);
    cqFile->writeUint32((uint32_t) opts.qualityValueOffset);

    // Write inverse quantization LUTs
    // compressedMappedQualSize_ += cqFile->writeQuantizers(quantizers_);

    if (opts.debug) {
        std::cerr << "New block. Quantizers:" << std::endl;
        for (auto &q : block.codeBooks) {
            std::cerr << "Quantizer " << ", " << q.size() << " steps:" << std::endl;
            size_t ctr = 0;
            for (auto &lut : q) {
                std::cerr << ctr << ": " << lut << std::endl;
                ctr++;
            }
        }
    }

    // Write unmapped quality values
    auto* uqv = (unsigned char*) unmappedQualityValues_.c_str();
    size_t uqvSize = unmappedQualityValues_.length();
    if (uqvSize > 0) {
        cqFile->writeUint8(0x01);
        cqFile->writeQualBlock(uqv, uqvSize);
    } else {
        cqFile->writeUint8(0x00);
    }

    // Write mapped quantizer indices
    std::string mqiString;
    for (auto const &mappedQuantizerIndex : block.quantizerIndices) {
        mqiString += std::to_string(mappedQuantizerIndex);
    }
    auto* mqi = (unsigned char*) mqiString.c_str();
    size_t mqiSize = mqiString.length();
    if (mqiSize > 0) {
        cqFile->writeUint8(0x01);
        cqFile->writeQualBlock(mqi, mqiSize);
    } else {
        cqFile->writeUint8(0x00);
    }

    // Write mapped quality value indices
    for (int i = 0; i < opts.quantizationMax - opts.quantizationMin + 1; ++i) {
        std::vector<uint8_t> mqviStream = block.stepindices[i];
        std::string mqviString;
        for (auto const &mqviInt : mqviStream) {
            mqviString += std::to_string(mqviInt);
        }
        auto* mqvi = (unsigned char*) mqviString.c_str();
        size_t mqviSize = mqviString.length();
        if (mqviSize > 0) {
            cqFile->writeUint8(0x01);
            cqFile->writeQualBlock(mqvi, mqviSize);
            // compressedMappedQualSize_ += cqFile->write(mqvi, mqviSize);
        } else {
            cqFile->writeUint8(0x00);
        }
    }

    return 0; // TODO: Add statistics
}


size_t readBlock(calq::CQFile* cqFile, calq::DecodingBlock* out, calq::DecodingSideInformation* side, std::string* unmapped) {
    size_t ret = 0;

    std::string buffer;

    // Read block parameters
    ret += cqFile->readUint32(&side->posOffset);
    ret += cqFile->readUint32(reinterpret_cast<uint32_t*>(&side->qualOffset));

    std::map<int, calq::Quantizer> quantizers_;

    for (size_t i = 0; i < quantizers_.size(); ++i) {
        out->codeBooks.emplace_back();
        for(size_t j = 0; j < quantizers_[i].inverseLut().size(); ++j) {
            out->codeBooks[i].push_back(quantizers_[i].inverseLut().at(j));
        }
    }

    // Read inverse quantization LUTs
    cqFile->readQuantizers(&quantizers_);

    // Read unmapped quality values
    uint8_t uqvFlags = 0;
    ret += cqFile->readUint8(&uqvFlags);
    if (uqvFlags & 0x01) { //NOLINT
        ret += cqFile->readQualBlock(unmapped);
    }

    // Read mapped quantizer indices
    uint8_t mqiFlags = 0;
    ret += cqFile->readUint8(&mqiFlags);
    if (mqiFlags & 0x1) { //NOLINT
        ret += cqFile->readQualBlock(&buffer);
        std::copy(buffer.begin(), buffer.end(), out->quantizerIndices.begin());
        buffer.clear();
    }

    // Read mapped quality value indices
    for (int i = 0; i < static_cast<int>(quantizers_.size()); ++i) {
        out->stepindices.emplace_back();
        uint8_t mqviFlags = 0;
        ret += cqFile->readUint8(&mqviFlags);
        if (mqviFlags & 0x1) { //NOLINT
            ret += cqFile->readQualBlock(&buffer);
            std::copy(buffer.begin(), buffer.end(), out->stepindices[i].begin());
            buffer.clear();
        }
    }

    return ret;
}

int main(int argc, char *argv[]){
    try
    {
        /*
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
        */
        // Compress or decompress

        calqapp::ProgramOptions ProgramOptions(argc, argv);
        ProgramOptions.validate();
        /*
        if (!options.decompress)
        {
            calq::SAMFileHandler sH(options.inputFileName);
            while (sH.readBlock(options.blockSize) != 0)
            {

                // Container to be filled by lib
                std::string reference; // TODO: read from file and cut for block
                std::string unmappedQualityScores; // TODO: Read
                uint64_t blockStartPosition = 0; // TODO: Calculate
                calq::EncodingOptions opts; //TODO: fill
                calq::EncodingBlock encBlock{sH.getQualityScores()};
                calq::DecodingBlock decBlock;
                calq::EncodingSideInformation encSide{
                        sH.getPositions(),
                        sH.getSequences(),
                        sH.getCigars(),
                        reference,
                        blockStartPosition
                };

                calq::encode(opts, encSide, encBlock, &decBlock);

                calq::CQFile file(options.outputFileName, calq::File::Mode::MODE_WRITE);

                writeBlock(opts, decBlock, encSide, unmappedQualityScores, &file);

            }
            CALQ_LOG("Finished encoding");
        }
        else
        {
            calq::DecodingBlock input;
            calq::EncodingBlock output;
            calq::DecodingSideInformation side; // TODO: Load from file
            std::string unmappedValues;

            calq::CQFile file(options.inputFileName, calq::File::Mode::MODE_READ);

            readBlock(&file, &input, &side, &unmappedValues);

            calq::decode(side, input, &output);

            //TODO: Write output & unmappedValues to file
            CALQ_LOG("Finished decoding");
        }
        /* } catch (TCLAP::ArgException &tclapException) {
             CALQ_ERROR("%s (argument: %s)", tclapException.error().c_str(), tclapException.argId().c_str());
             return EXIT_FAILURE;
         */
    }
    catch (const calq::ErrorException& errorException)
    {
        CALQ_ERROR("%s", errorException.what());
        return EXIT_FAILURE;
    }
    catch (const std::exception& stdException)
    {
        CALQ_ERROR("Fatal: %s", stdException.what());
        return EXIT_FAILURE;
    }
    catch (...)
    {
        CALQ_ERROR("Unknown error occurred");
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

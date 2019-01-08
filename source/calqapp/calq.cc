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
    cqFile->writeUint32(side.positions[0]);
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
        // Compress or decompress

        calqapp::ProgramOptions ProgramOptions(argc, argv);
        ProgramOptions.validate();

        // TO-Do: Fill structs below with information from SAMFileHandler

        if (!ProgramOptions.decompress)
        {
            calq::SAMFileHandler sH(ProgramOptions.inputFilePath);
            while (sH.readBlock(ProgramOptions.blockSize) != 0)
            {

                // Container to be filled by lib
                std::string reference; // TODO: read from file and cut for block
                std::string unmappedQualityScores; // TODO: Read
                calq::EncodingBlock encBlock{sH.getQualityScores()};
                calq::DecodingBlock decBlock;
                calq::EncodingSideInformation encSide{
                        sH.getPositions(),
                        sH.getSequences(),
                        sH.getCigars(),
                        reference
                };

                calq::encode(ProgramOptions.options, encSide, encBlock, &decBlock);

                calq::CQFile file(ProgramOptions.outputFilePath, calq::File::Mode::MODE_WRITE);

                writeBlock(ProgramOptions.options, decBlock, encSide, unmappedQualityScores, &file);

            }
            CALQ_LOG("Finished encoding");
        }
        else
        {
            calq::DecodingBlock input;
            calq::EncodingBlock output;
            calq::DecodingSideInformation side; // TODO: Load from file
            std::string unmappedValues;

            calq::CQFile file(ProgramOptions.inputFilePath, calq::File::Mode::MODE_READ);

            readBlock(&file, &input, &side, &unmappedValues);

            calq::decode(side, input, &output);

            //TODO: Write output & unmappedValues to file
            CALQ_LOG("Finished decoding");
        }
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

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
#include "calq/fasta_file.h"


#include "calqapp/SAMFileHandler.h"

#define STREAMOUT 0

size_t writeBlock(const calq::EncodingOptions& opts,
                  const calq::DecodingBlock& block,
                  const calq::EncodingSideInformation& side,
                  const std::string& unmappedQualityValues_,
                  calq::CQFile *cqFile
){

    size_t fileSize = 0;
    // Write block parameters
    fileSize += cqFile->writeUint32(side.positions[0] - 1);
    fileSize += cqFile->writeUint32((uint32_t) opts.qualityValueOffset);

    // Write inverse quantization LUTs
    fileSize += cqFile->writeQuantizers(block.codeBooks);
    
    // Write unmapped quality values
    auto *uqv = (unsigned char *) unmappedQualityValues_.c_str();
    size_t uqvSize = unmappedQualityValues_.length();
    if (uqvSize > 0)
    {
        fileSize += cqFile->writeUint8(0x01);

        if (STREAMOUT)
        {
            std::cerr << "unmapped qvalues:" << std::endl;

            std::cerr << unmappedQualityValues_;

            std::cerr << std::endl;
        }

        fileSize += cqFile->writeQualBlock(uqv, uqvSize);
    }
    else
    {
        fileSize += cqFile->writeUint8(0x00);
    }

    // Write mapped quantizer indices
    std::string mqiString;

    for (auto const& mappedQuantizerIndex : block.quantizerIndices)
    {
        mqiString += std::to_string(mappedQuantizerIndex);
    }

    if (STREAMOUT)
    {
        std::cerr << "quantizer indices:" << std::endl;

        std::cerr << mqiString;

        std::cerr << std::endl;
    }

    auto *mqi = (unsigned char *) mqiString.c_str();
    size_t mqiSize = mqiString.length();
    if (mqiSize > 0)
    {
        fileSize += cqFile->writeUint8(0x01);
        fileSize += cqFile->writeQualBlock(mqi, mqiSize);
    }
    else
    {
        fileSize += cqFile->writeUint8(0x00);
    }

    // Write mapped quality value indices
    for (int i = 0; i < opts.quantizationMax - opts.quantizationMin + 1; ++i)
    {
        std::vector<uint8_t> mqviStream = block.stepindices[i];
        std::string mqviString;


        for (auto const& mqviInt : mqviStream)
        {
            mqviString += std::to_string(mqviInt);
        }
        auto *mqvi = (unsigned char *) mqviString.c_str();
        size_t mqviSize = mqviString.length();

        if (STREAMOUT)
        {
            std::cerr << "Step indices" << i << ":" << std::endl;


            std::cerr << mqviString;

            std::cerr << std::endl;
        }
        if (mqviSize > 0)
        {
            fileSize += cqFile->writeUint8(0x01);
            fileSize += cqFile->writeQualBlock(mqvi, mqviSize);
        }
        else
        {
            fileSize += cqFile->writeUint8(0x00);
        }
    }

    return fileSize;
}


size_t readBlock(calq::CQFile *cqFile,
                 calq::DecodingBlock *out,
                 calq::DecodingSideInformation *side,
                 std::string *unmapped
){
    size_t ret = 0;

    std::string buffer;

    // Read block parameters
    ret += cqFile->readUint32(&side->posOffset);
    ret += cqFile->readUint32(reinterpret_cast<uint32_t *>(&side->qualOffset));

    std::map<int, calq::Quantizer> quantizers_;

    for (size_t i = 0; i < quantizers_.size(); ++i)
    {
        out->codeBooks.emplace_back();
        for (size_t j = 0; j < quantizers_[i].inverseLut().size(); ++j)
        {
            out->codeBooks[i].push_back(quantizers_[i].inverseLut().at(j));
        }
    }

    // Read inverse quantization LUTs
    cqFile->readQuantizers(&quantizers_);

    // Read unmapped quality values
    uint8_t uqvFlags = 0;
    ret += cqFile->readUint8(&uqvFlags);
    if (uqvFlags & 0x01)
    { //NOLINT
        ret += cqFile->readQualBlock(unmapped);
    }

    // Read mapped quantizer indices
    uint8_t mqiFlags = 0;
    ret += cqFile->readUint8(&mqiFlags);
    if (mqiFlags & 0x1)
    { //NOLINT
        ret += cqFile->readQualBlock(&buffer);
        std::copy(buffer.begin(), buffer.end(), out->quantizerIndices.begin());
        buffer.clear();
    }

    // Read mapped quality value indices
    for (int i = 0; i < static_cast<int>(quantizers_.size()); ++i)
    {
        out->stepindices.emplace_back();
        uint8_t mqviFlags = 0;
        ret += cqFile->readUint8(&mqviFlags);
        if (mqviFlags & 0x1)
        { //NOLINT
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

            std::unique_ptr<calq::FASTAFile> fastaFile;

            if (ProgramOptions.options.version == calq::EncodingOptions::Version::V2)
            {
                fastaFile = std::unique_ptr<calq::FASTAFile>(new calq::FASTAFile(ProgramOptions.referenceFilePath));
            }

            calq::CQFile file(ProgramOptions.outputFilePath, calq::File::Mode::MODE_WRITE);
            file.writeHeader(ProgramOptions.blockSize);

            while (sH.readBlock(ProgramOptions.blockSize) != 0)
            {

                // Container to be filled by lib
                std::string reference;
                std::string rname;
                if (ProgramOptions.options.version == calq::EncodingOptions::Version::V2)
                {
                    sH.getRname(&rname);
                    reference = fastaFile->getReferencesInRange(rname, sH.getRefStart(), sH.getRefEnd());
                    std::transform(reference.begin(), reference.end(),reference.begin(), ::toupper);
                }
                std::string unmappedQualityScores;
                sH.getUnmappedQualityScores(&unmappedQualityScores);
                calq::EncodingBlock encBlock;
                sH.getMappedQualityScores(&encBlock.qvalues);
                calq::DecodingBlock decBlock;

                calq::EncodingSideInformation encSide;
                sH.getPositions(&encSide.positions);
                sH.getSequences(&encSide.sequences);
                sH.getCigars(&encSide.cigars);
                encSide.reference.swap(reference);

                calq::encode(ProgramOptions.options, encSide, encBlock, &decBlock);

                writeBlock(ProgramOptions.options, decBlock, encSide, unmappedQualityScores, &file);

            }

            file.close();

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

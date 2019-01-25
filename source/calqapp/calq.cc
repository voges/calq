#include <chrono>

// -----------------------------------------------------------------------------

#include <boost/program_options.hpp>

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

size_t writeBlock(const calq::EncodingOptions& opts,
                  const calq::DecodingBlock& block,
                  const calq::SideInformation& side,
                  const std::string& unmappedQualityValues_,
                  bool STREAMOUT,
                  calqapp::CQFile *cqFile,
                  size_t *compressedSizeMapped,
                  size_t *compressedSizeUnmapped
){

    *compressedSizeMapped = 0;
    *compressedSizeUnmapped = 0;
    // Write block parameters
    *compressedSizeMapped += cqFile->writeUint32(side.positions[0] - 1);
    *compressedSizeMapped += cqFile->writeUint32((uint32_t) opts.qualityValueOffset);

    // Write inverse quantization LUTs
    *compressedSizeMapped += cqFile->writeQuantizers(block.codeBooks);

    // Write unmapped quality values
    auto *uqv = (unsigned char *) unmappedQualityValues_.c_str();
    size_t uqvSize = unmappedQualityValues_.length();
    if (uqvSize > 0) {
        *compressedSizeUnmapped += cqFile->writeUint8(0x01);

        if (STREAMOUT) {
            std::cerr << "unmapped qvalues:" << std::endl;

            std::cerr << unmappedQualityValues_;

            std::cerr << std::endl;
        }

        *compressedSizeUnmapped += cqFile->writeQualBlock(uqv, uqvSize);
    } else {
        *compressedSizeUnmapped += cqFile->writeUint8(0x00);
    }

    // Write mapped quantizer indices
    std::string mqiString;

    for (auto const& mappedQuantizerIndex : block.quantizerIndices) {
        mqiString += std::to_string(mappedQuantizerIndex);
    }

    if (STREAMOUT) {
        std::cerr << "quantizer indices:" << std::endl;

        std::cerr << mqiString;

        std::cerr << std::endl;
    }

    auto *mqi = (unsigned char *) mqiString.c_str();
    size_t mqiSize = mqiString.length();
    if (mqiSize > 0) {
        *compressedSizeMapped += cqFile->writeUint8(0x01);
        *compressedSizeMapped += cqFile->writeQualBlock(mqi, mqiSize);
    } else {
        *compressedSizeMapped += cqFile->writeUint8(0x00);
    }

    // Write mapped quality value indices
    for (int i = 0; i < opts.quantizationMax - opts.quantizationMin + 1; ++i) {
        std::vector<uint8_t> mqviStream = block.stepindices[i];
        std::string mqviString;


        for (auto const& mqviInt : mqviStream) {
            mqviString += std::to_string(mqviInt);
        }
        auto *mqvi = (unsigned char *) mqviString.c_str();
        size_t mqviSize = mqviString.length();

        if (STREAMOUT) {
            std::cerr << "Step indices" << i << ":" << std::endl;


            std::cerr << mqviString;

            std::cerr << std::endl;
        }
        if (mqviSize > 0) {
            *compressedSizeMapped += cqFile->writeUint8(0x01);
            *compressedSizeMapped += cqFile->writeQualBlock(mqvi, mqviSize);
        } else {
            *compressedSizeMapped += cqFile->writeUint8(0x00);
        }
    }

    return *compressedSizeUnmapped + *compressedSizeMapped;
}

// -----------------------------------------------------------------------------

size_t readBlock(calqapp::CQFile *cqFile,
                 calq::DecodingBlock *out,
                 calq::SideInformation *side,
                 std::string *unmapped
){
    out->codeBooks.clear();
    out->stepindices.clear();
    out->quantizerIndices.clear();

    unmapped->clear();

    size_t ret = 0;

    std::string buffer;

    // Read block parameters
    ret += cqFile->readUint32(&side->posOffset);
    ret += cqFile->readUint32(reinterpret_cast<uint32_t *>(&side->qualOffset));

    // Read inverse quantization LUTs
    cqFile->readQuantizers(&out->codeBooks);
    /*    
    for (size_t i = 0; i < quantizers_.size(); ++i)
    {
        out->codeBooks.emplace_back();
        for (size_t j = 0; j < quantizers_[i].inverseLut().size(); ++j)
        {
            out->codeBooks[i].push_back(quantizers_[i].inverseLut().at(j));
        }
    }
    */

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
        std::copy(buffer.begin(), buffer.end(), std::back_inserter(out->quantizerIndices));
        buffer.clear();
    }

    // Read mapped quality value indices
    for (int i = 0; i < static_cast<int>(out->codeBooks.size()); ++i) {
        out->stepindices.emplace_back();
        uint8_t mqviFlags = 0;
        ret += cqFile->readUint8(&mqviFlags);
        if (mqviFlags & 0x1) { //NOLINT
            ret += cqFile->readQualBlock(&buffer);
            std::copy(
                    buffer.begin(),
                    buffer.end(),
                    std::back_inserter(out->stepindices[i]));
            buffer.clear();
        }
    }

    return ret;
}

// -----------------------------------------------------------------------------

int main(int argc,
         char *argv[]
){
    try {
        // calq::setLogging(calq::loggingPresets::getSilent());
        // Compress or decompress

        calqapp::ProgramOptions ProgramOptions(argc, argv);
        ProgramOptions.validate();

        if (!ProgramOptions.decompress) {
            auto startTime = std::chrono::steady_clock::now();

            size_t compressedMappedQualSize = 0;
            size_t compressedUnmappedQualSize = 0;
            size_t uncompressedMappedQualSize = 0;
            size_t uncompressedUnmappedQualSize = 0;

            calqapp::SAMFileHandler sH(ProgramOptions.inputFilePath, ProgramOptions.referenceFilePath);

            calqapp::CQFile file(
                    ProgramOptions.outputFilePath,
                    calqapp::File::Mode::MODE_WRITE
            );
            file.writeHeader(ProgramOptions.blockSize);

            while (sH.readBlock(ProgramOptions.blockSize) != 0) {
                calq::SideInformation encSide;
                sH.getSideInformation(&encSide);

                calqapp::UnmappedInformation unmappedInfo;
                sH.getUnmappedBlock(&unmappedInfo);

                calq::EncodingBlock encBlock;
                sH.getMappedBlock(&encBlock);

                calq::DecodingBlock decBlock;


                calq::encode(
                        ProgramOptions.options,
                        encSide,
                        encBlock,
                        &decBlock
                );

                // Build single string out of unmapped q-values
                std::string unmappedString;
                for (const auto& s : unmappedInfo.unmappedQualityScores) {
                    unmappedString += s;
                }

                size_t mappedSize;
                size_t unmappedSize;
                writeBlock(
                        ProgramOptions.options,
                        decBlock,
                        encSide,
                        unmappedString,
                        ProgramOptions.debugStreams,
                        &file,
                        &mappedSize,
                        &unmappedSize
                );


                compressedMappedQualSize += mappedSize;//qualEncoder.compressedMappedQualSize();
                compressedUnmappedQualSize += unmappedSize;//qualEncoder.compressedUnmappedQualSize();
                for (const auto& q : encBlock.qvalues) {
                    uncompressedMappedQualSize += q.size();
                }
                uncompressedUnmappedQualSize += unmappedString.length();

            }

            file.close();

            auto stopTime = std::chrono::steady_clock::now();
            auto diffTime = stopTime - startTime;
            auto diffTimeMs = std::chrono::duration_cast<std::chrono::milliseconds>(diffTime).count();
            auto diffTimeS = std::chrono::duration_cast<std::chrono::seconds>(diffTime).count();
            auto diffTimeM = std::chrono::duration_cast<std::chrono::minutes>(diffTime).count();
            auto diffTimeH = std::chrono::duration_cast<std::chrono::hours>(diffTime).count();

            CALQ_LOG("COMPRESSION STATISTICS");
            CALQ_LOG("  Took %d ms ~= %d s ~= %d m ~= %d h",
                     (int) diffTimeMs,
                     (int) diffTimeS,
                     (int) diffTimeM,
                     (int) diffTimeH);
            const unsigned MB = 1 * 1000 * 1000;
            CALQ_LOG("  Speed (uncompressed size/time): %.2f MB/s",
                     ((static_cast<double>(uncompressedMappedQualSize + uncompressedUnmappedQualSize) /
                       static_cast<double>(MB))) / (static_cast<double>(diffTimeS)));
            CALQ_LOG("  Wrote %zu block(s)", sH.nrBlocksRead());
            CALQ_LOG("  Record(s):  %12zu", sH.nrRecordsRead());
            CALQ_LOG("    Mapped:   %12zu", sH.nrMappedRecordsRead());
            CALQ_LOG("    Unmapped: %12zu", sH.nrUnmappedRecordsRead());
            CALQ_LOG("  Uncompressed size: %12zu", uncompressedMappedQualSize + uncompressedUnmappedQualSize);
            CALQ_LOG("    Mapped:          %12zu", uncompressedMappedQualSize);
            CALQ_LOG("    Unmapped:        %12zu", uncompressedUnmappedQualSize);
            CALQ_LOG("  Compressed size: %12zu", file.nrWrittenBytes());
            CALQ_LOG("    File format:   %12zu", file.nrWrittenFileFormatBytes());
            CALQ_LOG("    Mapped:        %12zu", compressedMappedQualSize);
            CALQ_LOG("    Unmapped:      %12zu", compressedUnmappedQualSize);
            CALQ_LOG("  Compression ratio: %4.2f%%",
                     (double) file.nrWrittenBytes() *
                     100 /
                     (double) (uncompressedMappedQualSize + uncompressedUnmappedQualSize));
            CALQ_LOG("    Mapped:          %4.2f%%",
                     (double) compressedMappedQualSize * 100 / (double) (uncompressedMappedQualSize));
            CALQ_LOG("    Unmapped:        %4.2f%%",
                     (double) compressedUnmappedQualSize * 100 / (double) (uncompressedUnmappedQualSize));
            CALQ_LOG("  Compression factor: %4.2f",
                     (double) (uncompressedMappedQualSize + uncompressedUnmappedQualSize) /
                     (double) file.nrWrittenBytes());
            CALQ_LOG("    Mapped:           %4.2f",
                     (double) (uncompressedMappedQualSize) / (double) compressedMappedQualSize);
            CALQ_LOG("    Unmapped:         %4.2f",
                     (double) (uncompressedUnmappedQualSize) / (double) compressedUnmappedQualSize);
            CALQ_LOG("  Bits per quality value: %2.4f",
                     (static_cast<double>(file.nrWrittenBytes()) * 8) /
                     static_cast<double>(uncompressedMappedQualSize + uncompressedUnmappedQualSize));
            CALQ_LOG("    Mapped:               %2.4f",
                     (static_cast<double>(compressedMappedQualSize) * 8) /
                     static_cast<double>(uncompressedMappedQualSize));
            CALQ_LOG("    Unmapped:             %2.4f",
                     (static_cast<double>(compressedUnmappedQualSize) * 8) /
                     static_cast<double>(uncompressedUnmappedQualSize));
        } else {
            auto startTime = std::chrono::steady_clock::now();

            calqapp::CQFile file(
                    ProgramOptions.inputFilePath,
                    calqapp::File::Mode::MODE_READ
            );

            calqapp::File qualFile(
                    ProgramOptions.outputFilePath,
                    calqapp::File::Mode::MODE_WRITE
            );

            file.readHeader(&ProgramOptions.blockSize);

            calqapp::SAMFileHandler sH(
                    ProgramOptions.sideInformationFilePath,
                    ProgramOptions.referenceFilePath
            );

            while (sH.readBlock(ProgramOptions.blockSize) != 0) {
                //file side
                calq::DecodingBlock input;
                calq::EncodingBlock output;
                calq::SideInformation side;
                calqapp::UnmappedInformation unmappedInfo;

                sH.getSideInformation(&side);
                sH.getUnmappedBlock(&unmappedInfo);

                side.qualOffset = ProgramOptions.options.qualityValueOffset;
                input.quantizerIndices.clear();

                std::string unmappedValues;
                calq::DecodingOptions opts;
                readBlock(&file, &input, &side, &unmappedValues);
                calq::decode(opts, side, input, &output);


                auto mappedIt = output.qvalues.begin();
                auto unmappedPos = 0;
                auto unmappedSideIt = unmappedInfo.unmappedQualityScores.begin();
                for (const auto& b : unmappedInfo.mappedFlags) {
                    if (b) {
                        qualFile.write(mappedIt->c_str(), mappedIt->length());
                        qualFile.writeByte('\n');
                        ++mappedIt;
                    } else {
                        std::string read = unmappedValues.substr(unmappedPos, unmappedSideIt->length());
                        qualFile.write(read.c_str(), read.length());
                        qualFile.writeByte('\n');

                        unmappedPos += read.length();
                        ++unmappedSideIt;
                    }
                }
            }

            auto stopTime = std::chrono::steady_clock::now();
            auto diffTime = stopTime - startTime;
            auto diffTimeMs = std::chrono::duration_cast<std::chrono::milliseconds>(diffTime).count();
            auto diffTimeS = std::chrono::duration_cast<std::chrono::seconds>(diffTime).count();
            auto diffTimeM = std::chrono::duration_cast<std::chrono::minutes>(diffTime).count();
            auto diffTimeH = std::chrono::duration_cast<std::chrono::hours>(diffTime).count();

            CALQ_LOG("DECOMPRESSION STATISTICS");
            CALQ_LOG("  Took %d ms ~= %d s ~= %d m ~= %d h",
                     (int) diffTimeMs, (int) diffTimeS, (int) diffTimeM, (int) diffTimeH
            );
            const unsigned MB = 1 * 1000 * 1000;
            CALQ_LOG("  Speed (compressed size/time): %.2f MB/s",
                     ((static_cast<double>(file.nrReadBytes()) / static_cast<double>(MB))) /
                     (static_cast<double>(diffTimeS)));
            CALQ_LOG("  Decoded %zu block(s)", sH.nrBlocksRead());
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

    return EXIT_SUCCESS;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

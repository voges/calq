/**
 * @file calq-codec.cc
 */

#include "calq-codec.h"
#include <map>
#include <utility>
#include "errors.h"
#include "lloyd-max-quantizer.h"
#include "qual-decoder.h"
#include "qual-encoder.h"
#include "uniform-min-max-quantizer.h"

namespace calq {

// #include "calq/calq_coder.h"
//
// // -----------------------------------------------------------------------------
//
// #include "calqapp/SAMFileHandler.h"
// #include "calqapp/cq_file.h"
// #include "calqapp/fasta_file.h"
// #include "calqapp/logging.h"
// #include "calqapp/program_options.h"
//
// static int encode(const calqapp::ProgramOptions& programOptions) {
//     auto startTime = std::chrono::steady_clock::now();
//
//     size_t compressedMappedQualSize = 0;
//     size_t compressedUnmappedQualSize = 0;
//     size_t uncompressedMappedQualSize = 0;
//     size_t uncompressedUnmappedQualSize = 0;
//
//     calqapp::SAMFileHandler sH(programOptions.inputFilePath, programOptions.referenceFilePath);
//
//     calqapp::CQFile file(programOptions.outputFilePath, calqapp::File::Mode::MODE_WRITE);
//     file.writeHeader(programOptions.blockSize);
//
//     gabac::EncodingConfiguration config;
//     calqapp::File configurationFile("config.json", calqapp::File::Mode::MODE_READ);
//     std::string jsonInput("\0", configurationFile.size());
//     configurationFile.read(&jsonInput[0], jsonInput.size());
//     json2config(jsonInput, &config);
//
//     while (sH.readBlock(programOptions.blockSize) != 0) {
//         calq::SideInformation encSide;
//         sH.getSideInformation(&encSide);
//
//         calqapp::UnmappedInformation unmappedInfo;
//         sH.getUnmappedBlock(&unmappedInfo);
//
//         calq::EncodingBlock encBlock;
//         sH.getMappedBlock(&encBlock);
//
//         calq::DecodingBlock decBlock;
//
//         calq::encode(programOptions.options, encSide, encBlock, &decBlock);
//
//         // Build single string out of unmapped q-values
//         std::string unmappedString;
//         unmappedString = std::accumulate(unmappedInfo.unmappedQualityScores.begin(),
//                                          unmappedInfo.unmappedQualityScores.end(), unmappedString);
//
//         size_t mappedSize;
//         size_t unmappedSize;
//
//         file.writeBlock(programOptions.options, config, decBlock, encSide, unmappedString,
//         programOptions.debugStreams,
//                         &mappedSize, &unmappedSize);
//
//         compressedMappedQualSize += mappedSize;
//         compressedUnmappedQualSize += unmappedSize;
//
//         uncompressedMappedQualSize =
//             std::accumulate(encBlock.qvalues.begin(), encBlock.qvalues.end(), uncompressedMappedQualSize,
//                             [](const size_t& a, const std::string& b) -> size_t { return a + b.size(); });
//         uncompressedUnmappedQualSize += unmappedString.length();
//     }
//
//     file.close();
//
//     auto stopTime = std::chrono::steady_clock::now();
//     auto diffTime = stopTime - startTime;
//     auto diffTimeMs = std::chrono::duration_cast<std::chrono::milliseconds>(diffTime).count();
//     auto diffTimeS = std::chrono::duration_cast<std::chrono::seconds>(diffTime).count();
//     auto diffTimeM = std::chrono::duration_cast<std::chrono::minutes>(diffTime).count();
//     auto diffTimeH = std::chrono::duration_cast<std::chrono::hours>(diffTime).count();
//
//     CALQ_LOG("COMPRESSION STATISTICS");
//     CALQ_LOG("  Took %d ms ~= %d s ~= %d m ~= %d h", (int)diffTimeMs, (int)diffTimeS, (int)diffTimeM,
//     (int)diffTimeH); const unsigned MB = 1 * 1000 * 1000; CALQ_LOG(
//         "  Speed (uncompressed size/time): %.2f MB/s",
//         ((static_cast<double>(uncompressedMappedQualSize + uncompressedUnmappedQualSize) / static_cast<double>(MB)))
//         /
//             (static_cast<double>(diffTimeS)));
//     CALQ_LOG("  Wrote %zu block(s)", sH.nrBlocksRead());
//     CALQ_LOG("  Record(s):  %12zu", sH.nrRecordsRead());
//     CALQ_LOG("    Mapped:   %12zu", sH.nrMappedRecordsRead());
//     CALQ_LOG("    Unmapped: %12zu", sH.nrUnmappedRecordsRead());
//     CALQ_LOG("  Uncompressed size: %12zu", uncompressedMappedQualSize + uncompressedUnmappedQualSize);
//     CALQ_LOG("    Mapped:          %12zu", uncompressedMappedQualSize);
//     CALQ_LOG("    Unmapped:        %12zu", uncompressedUnmappedQualSize);
//     CALQ_LOG("  Compressed size: %12zu", file.nrWrittenBytes());
//     CALQ_LOG("    File format:   %12zu", file.nrWrittenFileFormatBytes());
//     CALQ_LOG("    Mapped:        %12zu", compressedMappedQualSize);
//     CALQ_LOG("    Unmapped:      %12zu", compressedUnmappedQualSize);
//     CALQ_LOG("  Compression ratio: %4.2f%%",
//              (double)file.nrWrittenBytes() * 100 / (double)(uncompressedMappedQualSize +
//              uncompressedUnmappedQualSize));
//     CALQ_LOG("    Mapped:          %4.2f%%",
//              (double)compressedMappedQualSize * 100 / (double)(uncompressedMappedQualSize));
//     CALQ_LOG("    Unmapped:        %4.2f%%",
//              (double)compressedUnmappedQualSize * 100 / (double)(uncompressedUnmappedQualSize));
//     CALQ_LOG("  Compression factor: %4.2f",
//              (double)(uncompressedMappedQualSize + uncompressedUnmappedQualSize) / (double)file.nrWrittenBytes());
//     CALQ_LOG("    Mapped:           %4.2f", (double)(uncompressedMappedQualSize) / (double)compressedMappedQualSize);
//     CALQ_LOG("    Unmapped:         %4.2f",
//              (double)(uncompressedUnmappedQualSize) / (double)compressedUnmappedQualSize);
//     CALQ_LOG("  Bits per quality value: %2.4f",
//              (static_cast<double>(file.nrWrittenBytes()) * 8) /
//                  static_cast<double>(uncompressedMappedQualSize + uncompressedUnmappedQualSize));
//     CALQ_LOG("    Mapped:               %2.4f",
//              (static_cast<double>(compressedMappedQualSize) * 8) / static_cast<double>(uncompressedMappedQualSize));
//     CALQ_LOG("    Unmapped:             %2.4f",
//              (static_cast<double>(compressedUnmappedQualSize) * 8) /
//              static_cast<double>(uncompressedUnmappedQualSize));
//
//     return EXIT_SUCCESS;
// }
//
// // -----------------------------------------------------------------------------
//
// static int decode(const calqapp::ProgramOptions& po) {
//     calqapp::ProgramOptions programOptions = po;
//
//     auto startTime = std::chrono::steady_clock::now();
//
//     calqapp::CQFile file(programOptions.inputFilePath, calqapp::File::Mode::MODE_READ);
//
//     calqapp::File qualFile(programOptions.outputFilePath, calqapp::File::Mode::MODE_WRITE);
//
//     file.readHeader(&programOptions.blockSize);
//
//     calqapp::SAMFileHandler sH(programOptions.sideInformationFilePath, programOptions.referenceFilePath);
//
//     gabac::EncodingConfiguration configuration;
//     calqapp::File configurationFile("config.json", calqapp::File::Mode::MODE_READ);
//     std::string jsonInput("\0", configurationFile.size());
//     configurationFile.read(&jsonInput[0], jsonInput.size());
//     json2config(jsonInput, &configuration);
//
//     while (sH.readBlock(programOptions.blockSize) != 0) {
//         calq::DecodingBlock input;
//         calq::EncodingBlock output;
//         calq::SideInformation side;
//         calqapp::UnmappedInformation unmappedInfo;
//
//         sH.getSideInformation(&side);
//         sH.getUnmappedBlock(&unmappedInfo);
//
//         side.qualOffset = programOptions.options.qualityValueOffset;
//         input.quantizerIndices.clear();
//
//         std::string unmappedValues;
//         calq::DecodingOptions opts;
//         file.readBlock(&input, &side, &unmappedValues, configuration);
//         calq::decode(opts, side, input, &output);
//
//         auto mappedIt = output.qvalues.begin();
//         auto unmappedPos = 0u;
//         auto unmappedSideIt = unmappedInfo.unmappedQualityScores.begin();
//         for (const auto& b : unmappedInfo.mappedFlags) {
//             if (b) {
//                 qualFile.write(mappedIt->c_str(), mappedIt->length());
//                 qualFile.writeByte('\n');
//                 ++mappedIt;
//             } else {
//                 std::string read = unmappedValues.substr(unmappedPos, unmappedSideIt->length());
//                 qualFile.write(read.c_str(), read.length());
//                 qualFile.writeByte('\n');
//
//                 unmappedPos += read.length();
//                 ++unmappedSideIt;
//             }
//         }
//     }
//
//     auto stopTime = std::chrono::steady_clock::now();
//     auto diffTime = stopTime - startTime;
//     auto diffTimeMs = std::chrono::duration_cast<std::chrono::milliseconds>(diffTime).count();
//     auto diffTimeS = std::chrono::duration_cast<std::chrono::seconds>(diffTime).count();
//     auto diffTimeM = std::chrono::duration_cast<std::chrono::minutes>(diffTime).count();
//     auto diffTimeH = std::chrono::duration_cast<std::chrono::hours>(diffTime).count();
//
//     CALQ_LOG("DECOMPRESSION STATISTICS");
//     CALQ_LOG("  Took %d ms ~= %d s ~= %d m ~= %d h", (int)diffTimeMs, (int)diffTimeS, (int)diffTimeM,
//     (int)diffTimeH); const unsigned MB = 1 * 1000 * 1000; CALQ_LOG("  Speed (compressed size/time): %.2f MB/s",
//              ((static_cast<double>(file.nrReadBytes()) / static_cast<double>(MB))) /
//              (static_cast<double>(diffTimeS)));
//     std::cout << "  Decoded %zu block(s)" << sH.nrBlocksRead() << std::endl;
//
//     return EXIT_SUCCESS;
// }

static uint32_t computeLength(const std::string& cigar) {
    // Compute 0-based first position and 0-based last position this record
    // is mapped to on the reference used for alignment
    uint32_t posMax = 0;

    size_t cigarIdx = 0;
    size_t cigarLen = cigar.length();
    uint32_t opLen = 0;  // length of current CIGAR operation

    for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++) {
        if (isdigit(cigar[cigarIdx])) {
            opLen = opLen * 10 + (uint32_t)cigar[cigarIdx] - (uint32_t)'0';
            continue;
        }
        switch (cigar[cigarIdx]) {
            case 'M':
            case '=':
            case 'X':
                posMax += opLen;
                break;
            case 'I':
            case 'S':
                break;
            case 'D':
            case 'N':
                posMax += opLen;
                break;
            case 'H':
            case 'P':
                break;  // these have been clipped
            default:
                throwErrorException("Bad CIGAR string");
        }
        opLen = 0;
    }
    return posMax;
}

void encode(const EncodingOptions& opt, const SideInformation& sideInformation, const EncodingBlock& input,
            DecodingBlock* output) {
    ProbabilityDistribution pdf(opt.qualityValueMin, opt.qualityValueMax);

    // Check quality value range
    for (auto const& samRecord : input.qvalues) {
        for (auto const& q : samRecord) {
            if ((static_cast<int>(q) - opt.qualityValueOffset) < opt.qualityValueMin) {
                throwErrorException("Quality value too small");
            }
            if ((static_cast<int>(q) - opt.qualityValueOffset) > opt.qualityValueMax) {
                throwErrorException("Quality value too large");
            }
            pdf.add((static_cast<size_t>(q) - opt.qualityValueOffset));
        }
    }

    std::map<int, Quantizer> quantizers;

    for (auto i = static_cast<int>(opt.quantizationMin); i <= static_cast<int>(opt.quantizationMax); ++i) {
        if (opt.quantizerType == QuantizerType::UNIFORM) {
            UniformMinMaxQuantizer quantizer(static_cast<const int&>(opt.qualityValueMin),
                                             static_cast<const int&>(opt.qualityValueMax), i);
            quantizers.insert(std::pair<int, Quantizer>(static_cast<const int&>(i - opt.quantizationMin), quantizer));
        } else if (opt.quantizerType == QuantizerType::LLOYD_MAX) {
            LloydMaxQuantizer quantizer(static_cast<size_t>(i));
            quantizer.build(pdf);
            quantizers.insert(std::pair<int, Quantizer>(static_cast<const int&>(i - opt.quantizationMin), quantizer));
        } else {
            throwErrorException("Quantization Type not supported");
        }
    }

    // Encode the quality values
    QualEncoder qualEncoder(opt, quantizers, output);
    for (size_t i = 0; i < sideInformation.positions.size(); ++i) {
        std::string ref;
        uint32_t len = computeLength(sideInformation.cigars[i]);
        if (opt.version == Version::V2) {
            ref = sideInformation.reference.substr(sideInformation.positions[i] - sideInformation.positions[0], len);
        }
        MinSamRecord r = {sideInformation.positions[i], sideInformation.positions[i] + len, input.qvalues[i],
                          sideInformation.cigars[i],    sideInformation.sequences[i],       ref};
        qualEncoder.addMappedRecordToBlock(r);
    }

    qualEncoder.finishBlock();
}

void decode(const DecodingOptions&, const SideInformation& sideInformation, const DecodingBlock& input,
            EncodingBlock* output) {
    // Decode the quality values
    QualDecoder qualDecoder(input, sideInformation.positions[0], sideInformation.qualOffset, output);
    output->qvalues.clear();
    for (size_t i = 0; i < sideInformation.positions.size(); ++i) {
        DecodingRead r = {sideInformation.positions[i], sideInformation.cigars[i]};
        qualDecoder.decodeMappedRecordFromBlock(r);
    }
}

}  // namespace calq

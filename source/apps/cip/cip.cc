#include <iostream>
#include "program-options.h"
#include "calq/errors.h"

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

int calq_main(int argc, char* argv[]) {
    try {
        cip::ProgramOptions programOptions(argc, argv);

        if (programOptions.help) {
            return 0;
        }
        // programOptions.validate();
        //
        // if (!programOptions.decompress) {
        //     return encode(programOptions);
        // } else {
        //     return decode(programOptions);
        // }
    } catch (const calq::ErrorException& errorException) {
        std::cerr << "CALQ: Error: " << errorException.what() << std::endl;
        return EXIT_FAILURE;
    } catch (const std::exception& stdException) {
        std::cerr << "CALQ: Fatal: " << stdException.what() << std::endl;
        return -1;
    } catch (...) {
        std::cerr << "CALQ: Fatal: Unknown error occurred" << std::endl;
        return -1;
    }

    return 0;
}

int main(int argc, char* argv[]) {
    int mainRc = calq_main(argc, argv);
    if (mainRc != 0) {
        std::cerr << "CALQ: Error: Failed to run CALQ" << std::endl;
    }

    // The C standard makes no guarantees as to when output to stdout or stderr
    // (standard error) is actually flushed. If e.g. stdout is directed to a
    // file and an error occurs while flushing the data (after program
    // termination), then the output may be lost.
    // Thus we explicitly flush stdout and stderr. On failure, we notify the OS
    // by returning with EXIT_FAILURE.
    if (fflush(stdout) == EOF) {
        return EXIT_FAILURE;
    }
    if (fflush(stderr) == EOF) {
        return EXIT_FAILURE;
    }

    // Return to the caller
    int callerRc = ((mainRc == 0) ? EXIT_SUCCESS : EXIT_FAILURE);
    std::cout << "Exiting with return code: " << callerRc << std::endl;
    return callerRc;
}

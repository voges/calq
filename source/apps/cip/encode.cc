#include "encode.h"
#include <chrono>
#include <iostream>
#include <thread>
#include <util/log.h>

namespace cip {

void encode(const ProgramOptions &programOptions) {
    (void)programOptions;  // Silence warning about unused variable

    auto startTime = std::chrono::steady_clock::now();

    auto sleepTimeS = std::chrono::seconds(2);
    LOG_INFO("waiting for " + std::to_string(sleepTimeS.count()) + "s");
    std::this_thread::sleep_for(sleepTimeS);

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

    auto stopTime = std::chrono::steady_clock::now();
    auto diffTime = stopTime - startTime;
    auto diffTimeMs = std::chrono::duration_cast<std::chrono::milliseconds>(diffTime).count();
    auto diffTimeS = std::chrono::duration_cast<std::chrono::seconds>(diffTime).count();
    auto diffTimeMin = std::chrono::duration_cast<std::chrono::minutes>(diffTime).count();
    auto diffTimeH = std::chrono::duration_cast<std::chrono::hours>(diffTime).count();

    LOG_INFO("took " + std::to_string(diffTimeMs) + " ms ~= " + std::to_string(diffTimeS) + " s ~= " + std::to_string(diffTimeMin) + " min ~= " + std::to_string(diffTimeH) + " h");

    //
    //     CALQ_LOG("COMPRESSION STATISTICS");
    //     CALQ_LOG("  Took %d ms ~= %d s ~= %d m ~= %d h", (int)diffTimeMs, (int)diffTimeS, (int)diffTimeM,
    //     (int)diffTimeH); const unsigned MB = 1 * 1000 * 1000; CALQ_LOG(
    //         "  Speed (uncompressed size/time): %.2f MB/s",
    //         ((static_cast<double>(uncompressedMappedQualSize + uncompressedUnmappedQualSize) /
    //         static_cast<double>(MB)))
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
    //              (double)(uncompressedMappedQualSize + uncompressedUnmappedQualSize) /
    //              (double)file.nrWrittenBytes());
    //     CALQ_LOG("    Mapped:           %4.2f", (double)(uncompressedMappedQualSize) /
    //     (double)compressedMappedQualSize); CALQ_LOG("    Unmapped:         %4.2f",
    //              (double)(uncompressedUnmappedQualSize) / (double)compressedUnmappedQualSize);
    //     CALQ_LOG("  Bits per quality value: %2.4f",
    //              (static_cast<double>(file.nrWrittenBytes()) * 8) /
    //                  static_cast<double>(uncompressedMappedQualSize + uncompressedUnmappedQualSize));
    //     CALQ_LOG("    Mapped:               %2.4f",
    //              (static_cast<double>(compressedMappedQualSize) * 8) /
    //              static_cast<double>(uncompressedMappedQualSize));
    //     CALQ_LOG("    Unmapped:             %2.4f",
    //              (static_cast<double>(compressedUnmappedQualSize) * 8) /
    //              static_cast<double>(uncompressedUnmappedQualSize));
    //
}

}  // namespace cip

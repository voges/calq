#include "decode.h"
#include <chrono>
#include <iostream>
#include <thread>
#include <util/log.h>
#include "common.h"

namespace cip {

void decode(const ProgramOptions& programOptions) {
    (void)programOptions;  // Silence warning about unused variable

    auto t1 = std::chrono::high_resolution_clock::now();

    auto sleepTimeS = std::chrono::seconds(2);
    std::cout << "Waiting for " << sleepTimeS.count() << "s" << std::endl;
    std::this_thread::sleep_for(sleepTimeS);

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

    auto now = std::chrono::high_resolution_clock::now();
    logProgress(now - t1);

    //
    //     CALQ_LOG("DECOMPRESSION STATISTICS");
    //     CALQ_LOG("  Took %d ms ~= %d s ~= %d m ~= %d h", (int)diffTimeMs, (int)diffTimeS, (int)diffTimeM,
    //     (int)diffTimeH); const unsigned MB = 1 * 1000 * 1000; CALQ_LOG("  Speed (compressed size/time): %.2f MB/s",
    //              ((static_cast<double>(file.nrReadBytes()) / static_cast<double>(MB))) /
    //              (static_cast<double>(diffTimeS)));
    //     std::cout << "  Decoded %zu block(s)" << sH.nrBlocksRead() << std::endl;
    //
    //     return EXIT_SUCCESS;
}

}  // namespace cip

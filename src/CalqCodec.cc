/** @file CalqCodec.cc
 *  @brief This file contains the implementations of the CalqEncoder and
 *         CalqDecoder classes.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "CalqCodec.h"
#include "Common/constants.h"
#include "Common/Exceptions.h"
#include "Common/helpers.h"
#include "IO/FASTA/FASTAFile.h"
#include <chrono>
#include <limits>

cq::CalqEncoder::CalqEncoder(const CLIOptions &cliOptions)
    : m_blockSize(cliOptions.blockSize)
    , m_cqFile(cliOptions.outputFileName, CQFile::MODE_WRITE)
    , m_force(cliOptions.force)
    , m_polyploidy(cliOptions.polyploidy)
    //, m_qualEncoder(cqFile, cliOptions)
    , m_referenceFileNames(cliOptions.referenceFileNames)
    , m_samFile(cliOptions.inputFileName)
{
    // Check arguments
    if (cliOptions.blockSize < 1) {
        throwErrorException("blockSize must be greater than zero");
    }
    if (cliOptions.inputFileName.empty() == true) {
        throwErrorException("inputFileName is empty");
    }
    if (cliOptions.outputFileName.empty() == true) {
        throwErrorException("outputFileName is empty");
    }
    if (cliOptions.polyploidy < 1) {
        throwErrorException("polyploidy must be greater than zero");
    }
    //if (referenceFileNames.empty() == true) {
    //    throwErrorException("referenceFileNames is empty");
    //}

    // Check and, in case they are provided, get reference sequences
    if (m_referenceFileNames.empty() == true) {
        CQ_LOG("No reference file name(s) given - operating without reference sequence(s)");
    } else {
        CQ_LOG("Looking in %zu reference file(s) for reference sequence(s)", m_referenceFileNames.size());
        for (auto const &referenceFileName : m_referenceFileNames) {
            CQ_LOG("  Parsing reference file: %s", referenceFileName.c_str());
            FASTAFile fastaFile(referenceFileName);
            CQ_LOG("  Found %zu reference(s):", fastaFile.references.size());
            for (auto const &reference : fastaFile.references) {
                CQ_LOG("    %s (length: %zu)", reference.first.c_str(), reference.second.length());
            }
        }
    }
}

cq::CalqEncoder::~CalqEncoder(void)
{
    // empty
}

void cq::CalqEncoder::encode(void)
{
    auto startTime = std::chrono::steady_clock::now();

    size_t uncompressedSize = 0;
    size_t compressedSize = 0;

    compressedSize += m_cqFile.writeHeader();

    while (m_samFile.readBlock(m_blockSize) != 0) {
        CQ_LOG("Training Lloyd-Max quantizer(s)");

        int qMin = std::numeric_limits<int>::max();
        int qMax = std::numeric_limits<int>::min();

        for (auto const &samRecord : m_samFile.currentBlock.records) {
            uncompressedSize += samRecord.qual.length();

            for (auto const &q : samRecord.qual) {
                if ((int)q < qMin) qMin = q;
                if ((int)q > qMax) qMax = q;
            }
        }

        CQ_LOG("qMin: %d, qMax: %d", qMin, qMax);

        CQ_LOG("Encoding quality values");
        //qualEncoder.startBlock();
        for (auto const &samRecord : m_samFile.currentBlock.records) {
            if (samRecord.isMapped() == true) {
                //qualEncoder.addMappedRecordToBlock(samRecord);
            } else {
                //qualEncoder.addUnmappedRecordToBlock(samRecord);
            }
        }
        //qualEncoder.finishBlock();
    }

    auto stopTime = std::chrono::steady_clock::now();
    auto diffTime = stopTime - startTime;
    auto diffTimeMs = std::chrono::duration_cast<std::chrono::milliseconds>(diffTime).count();
    auto diffTimeS = std::chrono::duration_cast<std::chrono::seconds>(diffTime).count();
    auto diffTimeM = std::chrono::duration_cast<std::chrono::minutes>(diffTime).count();
    auto diffTimeH = std::chrono::duration_cast<std::chrono::hours>(diffTime).count();

    //qualEncoder.printStats();

    CQ_LOG("COMPRESSION STATISTICS");
    CQ_LOG("  Took %d ms ~= %d s ~= %d m ~= %d h", (int)diffTimeMs, (int)diffTimeS, (int)diffTimeM, (int)diffTimeH);
    CQ_LOG("  Compressed %zu mapped + %zu unmapped = %zu record(s) in %zu block(s)", m_samFile.numMappedRecordsRead(), m_samFile.numUnmappedRecordsRead(), m_samFile.numRecordsRead(), m_samFile.numBlocksRead());
    CQ_LOG("  Uncompressed size: %zu", uncompressedSize);
    CQ_LOG("  Compressed size: %zu", compressedSize);
    CQ_LOG("  Compression ratio: %.2f%%", (double)compressedSize*100/(double)uncompressedSize);
    CQ_LOG("  Compression factor: %.2f", (double)uncompressedSize/(double)compressedSize);
    CQ_LOG("  Bits per quality value: %.4f", ((double)compressedSize * 8)/(double)uncompressedSize);
    CQ_LOG("  Speed (uncompressed size/time): %.2f MB/s", ((double)(uncompressedSize/MB))/(double)((double)diffTimeMs/1000));
}

cq::CalqDecoder::CalqDecoder(const CLIOptions &cliOptions)
    : m_cqFile(cliOptions.inputFileName, CQFile::MODE_READ)
    , m_qualFile(cliOptions.outputFileName, File::MODE_WRITE)
    //, m_qualDecoder(cqFile, qualFile, cliOptions)
    , m_sideInformationFile(cliOptions.sideInformationFileName)
{
    // Check arguments
    if (cliOptions.inputFileName.empty() == true) {
        throwErrorException("cliOptions.inputFileName is empty");
    }
    if (cliOptions.outputFileName.empty() == true) {
        throwErrorException("cliOptions.outputFileName is empty");
    }
    if (cliOptions.sideInformationFileName.empty() == true) {
        throwErrorException("cliOptions.sideInformationFileName is empty");
    }
}

cq::CalqDecoder::~CalqDecoder(void)
{
    // empty
}

void cq::CalqDecoder::decode(void)
{
    auto startTime = std::chrono::steady_clock::now();

    size_t compressedSize = 0;

    compressedSize += m_cqFile.readHeader();

    size_t numBlocks = 0;

    //while (cqFile.eof() == false) {
        //
        //numBlocks++;
    //}

    auto stopTime = std::chrono::steady_clock::now();
    auto diffTime = stopTime - startTime;
    auto diffTimeMs = std::chrono::duration_cast<std::chrono::milliseconds>(diffTime).count();
    auto diffTimeS = std::chrono::duration_cast<std::chrono::seconds>(diffTime).count();
    auto diffTimeM = std::chrono::duration_cast<std::chrono::minutes>(diffTime).count();
    auto diffTimeH = std::chrono::duration_cast<std::chrono::hours>(diffTime).count();

    //qualDecoder.printStats();

    CQ_LOG("DECOMPRESSION STATISTICS");
    CQ_LOG("  Took %d ms ~= %d s ~= %d m ~= %d h", (int)diffTimeMs, (int)diffTimeS, (int)diffTimeM, (int)diffTimeH);
    CQ_LOG("  Decoded %zu block(s)", numBlocks);
    CQ_LOG("  Compressed size: %zu", compressedSize);
    CQ_LOG("  Speed (compressed size/time): %.2f MB/s", ((double)(compressedSize/MB))/(double)((double)diffTimeMs/1000));
}

// void MappedRecord::extractObservations(const uint32_t &observedPosMin,
//                                        std::deque<std::string> &observedNucleotides,
//                                        std::deque<std::string> &observedQualityValues)
// {
//     size_t cigarIdx = 0;
//     size_t cigarLen = cigar.length();
//     size_t opLen = 0; // length of current CIGAR operation
//     size_t idx = 0;
//     size_t observedIdx = posMin - observedPosMin;
// 
//     for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++) {
//         if (isdigit(cigar[cigarIdx])) {
//             opLen = opLen*10 + (size_t)cigar[cigarIdx] - (size_t)'0';
//             continue;
//         }
// 
//         switch (cigar[cigarIdx]) {
//         case 'M':
//         case '=':
//         case 'X':
//             // Add matching parts
//             for (size_t i = 0; i < opLen; i++) {
//                 observedNucleotides[observedIdx] += nucleotides[idx];
//                 observedQualityValues[observedIdx] += qualityValues[idx];
//                 idx++;
//                 observedIdx++;
//             }
//             break;
//         case 'I':
//         case 'S':
//             idx += opLen;
//             break;
//         case 'D':
//         case 'N':
//             observedIdx += opLen;
//             break;
//         case 'H':
//         case 'P':
//             break; // these have been clipped
//         default:
//             LOG("CIGAR string: %s", cigar.c_str());
//             throwErrorException("Bad CIGAR string");
//         }
// 
//         opLen = 0;
//     }
// }


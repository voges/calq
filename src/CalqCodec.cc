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
#include "IO/FASTA.h"
#include <chrono>

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
        LOG("No reference file name(s) given - operating without reference sequence(s)");
    } else {
        LOG("Looking in %zu reference file(s) for reference sequence(s)", m_referenceFileNames.size());
        for (auto const &referenceFileName : m_referenceFileNames) {
            LOG("  Parsing reference file: %s", referenceFileName.c_str());
            FASTAFile fastaFile(referenceFileName);
            LOG("  Found %zu reference(s):", fastaFile.references.size());
            for (auto const &reference : fastaFile.references) {
                LOG("    %s (length: %zu)", reference.first.c_str(), reference.second.length());
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
    size_t numBlocks = 0;
    size_t numRecords = 0;
    size_t numMappedRecords = 0;
    size_t numUnmappedRecords = 0;

    compressedSize += m_cqFile.writeHeader();

    while (m_samFile.eof() == false) {
        m_samFile.readBlock(m_blockSize);

        // 1st iteration
        LOG("Training Lloyd-Max quantizer(s)");
        for (auto const &samRecord : m_samFile.currentBlock.records) {
            uncompressedSize += samRecord.qual.length();
        }

        // 2nd iteration
        LOG("Encoding quality values");
        //qualEncoder.startBlock();
        for (auto const &samRecord : m_samFile.currentBlock.records) {
            if (samRecord.isMapped() == true) {
                //qualEncoder.addMappedRecordToBlock(samRecord);
            } else {
                //qualEncoder.addUnmappedRecordToBlock(samRecord);
            }
        }
        //qualEncoder.finishBlock();

        numMappedRecords += m_samFile.currentBlock.numUnmappedRecords();
        numUnmappedRecords += m_samFile.currentBlock.numMappedRecords();
        numRecords += m_samFile.currentBlock.numRecords();
        numBlocks++;
    }

    auto stopTime = std::chrono::steady_clock::now();
    auto diffTime = stopTime - startTime;
    auto diffTimeMs = std::chrono::duration_cast<std::chrono::milliseconds>(diffTime).count();
    auto diffTimeS = std::chrono::duration_cast<std::chrono::seconds>(diffTime).count();
    auto diffTimeM = std::chrono::duration_cast<std::chrono::minutes>(diffTime).count();
    auto diffTimeH = std::chrono::duration_cast<std::chrono::hours>(diffTime).count();

    //qualEncoder.printStats();

    LOG("COMPRESSION STATISTICS");
    LOG("  Took %ld ms ~= %ld s ~= %d m ~= %d h", diffTimeMs, diffTimeS, diffTimeM, diffTimeH);
    LOG("  Compressed %zu mapped + %zu unmapped = %zu record(s) in %zu block(s)", numMappedRecords, numUnmappedRecords, numRecords, numBlocks);
    LOG("  Uncompressed size: %zu", uncompressedSize);
    LOG("  Compressed size: %zu", compressedSize);
    LOG("  Compression ratio: %.2f%%", (double)compressedSize*100/(double)uncompressedSize);
    LOG("  Compression factor: %.2f", (double)uncompressedSize/(double)compressedSize);
    LOG("  Bits per quality value: %.4f", ((double)compressedSize * 8)/(double)uncompressedSize);
    LOG("  Speed (uncompressed size/time): %.2f MB/s", ((double)(uncompressedSize/MB))/(double)((double)diffTimeMs/1000));
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

    LOG("DECOMPRESSION STATISTICS");
    LOG("  Took %ld ms ~= %ld s ~= %d m ~= %d h", diffTimeMs, diffTimeS, diffTimeM, diffTimeH);
    LOG("  Decoded %zu block(s)", numBlocks);
    LOG("  Compressed size: %zu", compressedSize);
    LOG("  Speed (compressed size/time): %.2f MB/s", ((double)(compressedSize/MB))/(double)((double)diffTimeMs/1000));
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


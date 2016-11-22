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
#include "cmake_config.h"
#include "Common/constants.h"
#include "Common/Exceptions.h"
#include "Common/helpers.h"
#include <chrono>
#include <iostream>
#include <string.h>

#include "IO/FASTAFile.h"

// static bool samRecordIsMapped(const SAMRecord &samRecord)
// {
//     if (   ((samRecord.flag & 0x4) != 0)
//         || (strlen(samRecord.rname) == 0)
//         || ((samRecord.rname[0] == '*') && (strlen(samRecord.rname) == 1))
//         || (samRecord.pos == 0)
//         || (strlen(samRecord.cigar) == 0)
//         || ((samRecord.cigar[0] == '*') && (strlen(samRecord.rname) == 1))
//         || (strlen(samRecord.seq) == 0)
//         || ((samRecord.seq[0] == '*') && (strlen(samRecord.rname) == 1))
//         || (strlen(samRecord.qual) == 0)
//         || ((samRecord.qual[0] == '*') && (strlen(samRecord.rname) == 1))) {
//         return false;
//     }
// 
//     return true;
// }

CalqEncoder::CalqEncoder(const CLIOptions &cliOptions)
    : force(cliOptions.force)
    , samFile(cliOptions.inputFileName, "rb", cliOptions.blockSize)
    , cqFile(cliOptions.outputFileName, "w")
    , blockSize(cliOptions.blockSize)
    , polyploidy(cliOptions.polyploidy)
    , qualityValueMax(cliOptions.qualityValueMax)
    , qualityValueMin(cliOptions.qualityValueMin)
    , referenceFileNames(cliOptions.referenceFileNames)


    //, qualEncoder(cqFile, cliOptions)
{
    // Check class initialization
    if (cliOptions.inputFileName.empty() == true) {
        throwErrorException("inputFileName is empty");
    }
    if (cliOptions.outputFileName.empty() == true) {
        throwErrorException("outputFileName is empty");
    }
    if (blockSize < 1) {
        throwErrorException("blockSize must be greater than zero");
    }
    if (polyploidy < 1) {
        throwErrorException("polyploidy must be greater than zero");
    }
    if (referenceFileNames.empty() == true) {
        throwErrorException("referenceFileNames is empty");
    }

    // Check and, in case they are provided, get reference sequences
    if (referenceFileNames.empty() == true) {
        LOG("No reference file name(s) given - operating without reference sequence(s)");
    } else {
        LOG("Looking in %zu reference file(s) for reference sequence(s)", referenceFileNames.size());
        for (std::string const &referenceFileName : referenceFileNames) {
            LOG("  Parsing reference file: %s", referenceFileName.c_str());
            //FASTAFile fastaFile(referenceFileName, "r");
            //LOG("  Found %zu reference(s):", fastaFile.references.size());
            //for (auto const &reference : fastaFile.references) {
            //    LOG("    %s (length: %zu)", reference.first.c_str(), reference.second.length());
            //}
        }
    }
    
    // TODO: check first read block
    LOG("1st block contains %d SAM records", samFile.block.size());
    for (SAMRecord const &samRecord : samFile.block) {
        //samRecord.print();
    }
    
    while (samFile.eof() == false) {
        samFile.readBlock();
        LOG("\nREAD BLOCK\n");
        for (SAMRecord const &samRecord : samFile.block) {
            //samRecord.print();
        }
    }
}

CalqEncoder::~CalqEncoder(void)
{
    // empty
}

void CalqEncoder::encode(void)
{
    auto startTime = std::chrono::steady_clock::now();

    size_t uncompressedSize = 0;
    size_t compressedSize = 0;
    size_t numBlocks = 0;
    size_t numRecords = 0;
    size_t numMappedRecords = 0;
    size_t numUnmappedRecords = 0;
    size_t numRecordsInBlock = 0;
// 
//     size_t fileHeaderSize = writeFileHeader();
// 
//     std::string rnamePrev("");
//     uint32_t posPrev = 0;
//     //qualEncoder.startBlock();
// 
//     do {
//         uncompressedSize += strlen(samParser.curr.qual);
// 
//         if (numRecordsInBlock == blockSize) {
//             // Start a new block
//             //compressedSize += qualEncoder.finishBlock();
//             numBlocks++;
//             //qualEncoder.startBlock();
//             numRecordsInBlock = 0;
//             rnamePrev = "";
//             posPrev = 0;
//         }
// 
//         if (samRecordIsMapped(samParser.curr) == true) {
//             if (rnamePrev.empty() == false) {
//                 if (samParser.curr.rname != rnamePrev) {
//                     // Start a new block
//                     //compressedSize += qualEncoder.finishBlock();
//                     numBlocks++;
//                     //qualEncoder.startBlock();
//                     numRecordsInBlock = 0;
//                     rnamePrev = "";
//                     posPrev = 0;
//                 } else {
//                     if (samParser.curr.pos < posPrev) {
//                         throwErrorException("SAM file is not sorted");
//                     }
//                 }
//             }
//             posPrev = samParser.curr.pos;
//             rnamePrev = samParser.curr.rname;
//             //qualEncoder.addMappedRecordToBlock(samParser.curr);
//             numMappedRecords++;
//         } else {
//             //qualEncoder.addUnmappedRecordToBlock(samParser.curr);
//             numUnmappedRecords++;
//         }
// 
//         numRecords++;
//         numRecordsInBlock++;
//     } while (samParser.next());
// 
//     compressedSize += qualEncoder.finishBlock();
//     numBlocks++;
// 
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
    //LOG("  Compressed size: %zu (+ file header size: %zu)", compressedSize, fileHeaderSize);
    LOG("  Compression ratio: %.2f%%", (double)compressedSize*100/(double)uncompressedSize);
    LOG("  Compression factor: %.2f", (double)uncompressedSize/(double)compressedSize);
    LOG("  Bits per quality value: %.4f", ((double)compressedSize * 8)/(double)uncompressedSize);
    LOG("  Speed (uncompressed size/time): %.2f MB/s", ((double)(uncompressedSize/MB))/(double)((double)diffTimeMs/1000));
}

size_t CalqEncoder::writeFileHeader(void)
{
    size_t ret = 0;

    const size_t magicSize = 5;
    char magic[magicSize] = "calq";
    const size_t versionSize = 6;
    char version[versionSize] = VERSION;
    ret += cqFile.write(magic, magicSize);
    ret += cqFile.write(version, versionSize);
    ret += cqFile.writeUint32(polyploidy);

    return ret;
}

CalqDecoder::CalqDecoder(const CLIOptions &cliOptions)
//     : CalqCodec(cliOptions.inFileName, cliOptions.outFileName)
//     , cqFile(cliOptions.inFileName, "r")
//     , qualFile(cliOptions.outFileName, "w")
//     , qualDecoder(cqFile, qualFile, cliOptions)
//     , samParser(cliOptions.alignmentsFileName)
{
//     if (cliOptions.inFileName.empty() == true) {
//         throwErrorException("No input file name given");
//     }
//     if (cliOptions.outFileName.empty() == true) {
//         throwErrorException("No output file name given");
//     }
//     if (cliOptions.alignmentsFileName.empty() == true) {
//         throwErrorException("No alignments (i.e. SAM) file name given");
//     }
}

CalqDecoder::~CalqDecoder(void)
{
    // empty
}

void CalqDecoder::decode(void)
{
//     auto startTime = std::chrono::steady_clock::now();
// 
//     size_t compressedSize = 0;
//     size_t fileHeaderSize = readFileHeader();
// 
//     size_t numBlocks = 0;
//     while (cqFile.tell() < cqFile.size()) {
//         //std::cout << ME << "Decoding block " << numBlocks << std::endl;
//         compressedSize += qualDecoder.decodeBlock();
// 
//         // Get numRecordsInBlock from SAM file
// //         do {
// //             if (samRecordIsMapped(samParser.curr) == true) {
// //
// //             } else {
// //
// //             }
// //         } while (samParser.next());
// 
//         numBlocks++;
//     }
// 
//     auto stopTime = std::chrono::steady_clock::now();
//     auto diffTime = stopTime - startTime;
//     auto diffTimeMs = std::chrono::duration_cast<std::chrono::milliseconds>(diffTime).count();
//     auto diffTimeS = std::chrono::duration_cast<std::chrono::seconds>(diffTime).count();
//     auto diffTimeM = std::chrono::duration_cast<std::chrono::minutes>(diffTime).count();
//     auto diffTimeH = std::chrono::duration_cast<std::chrono::hours>(diffTime).count();
// 
//     qualDecoder.printStats();
// 
//     LOG("DECOMPRESSION STATISTICS");
//     LOG("  Took %ld ms ~= %ld s ~= %d m ~= %d h", diffTimeMs, diffTimeS, diffTimeM, diffTimeH);
//     LOG("  Decoded %zu block(s)", numBlocks);
//     LOG("  Compressed size: %zu (+ file header size: %zu)", compressedSize, fileHeaderSize);
//     LOG("  Speed (compressed size/time): %.2f MB/s", ((double)(compressedSize/MB))/(double)((double)diffTimeMs/1000));
}
// 
// size_t CalqDecoder::readFileHeader(void)
// {
//     size_t ret = 0;
// 
//     const size_t magicSize = 5;
//     char magic[magicSize];
//     const size_t versionSize = 6;
//     char version[versionSize];
//     unsigned int polyploidy;
//     ret += cqFile.read(magic, magicSize);
//     ret += cqFile.read(version, versionSize);
//     ret += cqFile.readUint32(&polyploidy);
// 
//     if (strncmp(version, VERSION, 5) != 0) {
//         LOG("Program version: %s", VERSION);
//         LOG("File version: %s", version);
//         throwErrorException("CQ file was compressed with another version");
//     }
// 
//     LOG("Magic: %s", magic);
//     LOG("Version: %s",  version);
//     LOG("Polyploidy: %u",  polyploidy);
// 
//     return ret;
// }


/** @file CalqCodec.cc
 *  @brief This file contains the implementations of the CalqEncoder and
 *         CalqDecoder class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "CalqCodec.h"
#include "cmake_config.h"
#include "Common/Exceptions.h"
#include "Common/log.h"
#include <chrono>
#include <iostream>
#include <string.h>

static bool samRecordIsMapped(const SAMRecord &samRecord)
{
    if (   ((samRecord.flag & 0x4) != 0)
        || (strlen(samRecord.rname) == 0 || samRecord.rname[0] == '*')
        || (samRecord.pos == 0)
        || (strlen(samRecord.cigar) == 0 || samRecord.cigar[0] == '*')
        || (strlen(samRecord.seq) == 0 || samRecord.seq[0] == '*')
        || (strlen(samRecord.qual) == 0 || samRecord.qual[0] == '*')) {
        return false;
    }

    return true;
}

CalqCodec::CalqCodec(const std::string &inFileName,
                     const std::string &outFileName)
    : inFileName(inFileName)
    , outFileName(outFileName)
{
    if (inFileName.empty()) {
        throwErrorException("No input file name given");
    }
    if (outFileName.empty()) {
        throwErrorException("No output file name given");
    }
}

CalqCodec::~CalqCodec(void)
{
    // empty
}

CalqEncoder::CalqEncoder(const CLIOptions &cliOptions)
    : CalqCodec(cliOptions.inFileName, cliOptions.outFileName)
    , blockSize(cliOptions.blockSize)
    , cqFile(cliOptions.outFileName, "w")
    , fastaReferences()
    , polyploidy(cliOptions.polyploidy)
    , qualEncoder(cqFile, cliOptions)
    , samParser(cliOptions.inFileName)
{
    if (cliOptions.inFileName.empty()) {
        throwErrorException("No input file name given");
    }
    if (cliOptions.outFileName.empty()) {
        throwErrorException("No output file name given");
    }
    if (cliOptions.blockSize < 1) {
        throwErrorException("Block size must be greater than zero");
    }
    if (cliOptions.polyploidy < 1) {
        throwErrorException("Polyploidy must be greater than zero");
    }
    if (cliOptions.refFileNames.size() == 0) {
        throwErrorException("No reference file names given");
    }

    // Get reference sequences
    FASTAParser fastaParser;
    for (auto const &fastaFileName : cliOptions.refFileNames) {
        LOG("Parsing reference file: %s", fastaFileName.c_str());
        fastaParser.parseFile(fastaFileName, fastaReferences);
    }

    // Print info about found reference sequences
    LOG("Found the following reference sequence(s):");
    for (auto const &fastaReference : fastaReferences) {
        LOG("  %s (sequence length: %lu)", fastaReference.header.c_str(), fastaReference.sequence.size());
    }

    // Now pass the FASTA references to the qualEncoder
    qualEncoder.fastaReferences = fastaReferences;
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

    size_t fileHeaderSize = writeFileHeader();

    // Send all records to the QualEncoder
    std::string rnamePrev("");
    uint32_t posPrev = 0;
    qualEncoder.startBlock();

    do {
        uncompressedSize += strlen(samParser.curr.qual);

        if (numRecordsInBlock == blockSize) {
            // Start a new block
            compressedSize += qualEncoder.finishBlock();
            numBlocks++;
            qualEncoder.startBlock();
            numRecordsInBlock = 0;
            rnamePrev = "";
            posPrev = 0;
        }

        if (samRecordIsMapped(samParser.curr) == true) {
            if (rnamePrev.empty() == false) {
                if (samParser.curr.rname != rnamePrev) {
                    // Start a new block
                    compressedSize += qualEncoder.finishBlock();
                    numBlocks++;
                    qualEncoder.startBlock();
                    numRecordsInBlock = 0;
                    rnamePrev = "";
                    posPrev = 0;
                } else {
                    if (samParser.curr.pos < posPrev) {
                        throwErrorException("SAM file is not sorted");
                    }
                }
            }
            posPrev = samParser.curr.pos;
            rnamePrev = samParser.curr.rname;
            qualEncoder.addMappedRecordToBlock(samParser.curr);
            numMappedRecords++;
        } else {
            qualEncoder.addUnmappedRecordToBlock(samParser.curr);
            numUnmappedRecords++;
        }

        numRecords++;
        numRecordsInBlock++;
    } while (samParser.next());

    compressedSize += qualEncoder.finishBlock();
    numBlocks++;

    // Print summary
    auto stopTime = std::chrono::steady_clock::now();
    auto diffTime = stopTime - startTime;
    auto diffTimeMs = std::chrono::duration_cast<std::chrono::milliseconds>(diffTime).count();
    auto diffTimeS = std::chrono::duration_cast<std::chrono::seconds>(diffTime).count();
    std::cout << ME << "COMPRESSION STATISTICS" << std::endl;
    std::cout << ME << "  Took " << diffTimeMs << " ms" << " ~= " << diffTimeS << " s" << std::endl;
    std::cout << ME << "  Compressed " << numMappedRecords << " mapped + " << numUnmappedRecords << " unmapped = " << numRecords << " record(s) in " << numBlocks << " block(s)" << std::endl;
    std::cout << ME << "  Uncompressed size: " << uncompressedSize << std::endl;
    std::cout << ME << "  Compressed size: " << compressedSize << " (+ file header size: " << fileHeaderSize << ")" << std::endl;
    std::cout << ME << "  Compression ratio: " << (double)compressedSize*100/(double)uncompressedSize << "%" << std::endl;
    std::cout << ME << "  Compression factor: " << (double)uncompressedSize/(double)compressedSize << std::endl;
    std::cout << ME << "  Bits per quality value: " << ((double)compressedSize * 8)/(double)uncompressedSize << std::endl;
    std::cout << ME << "  Speed (uncompressed size/time): " << ((double)(uncompressedSize/MB))/(double)((double)diffTimeMs/1000) << " MB/s" << std::endl;
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
    : CalqCodec(cliOptions.inFileName, cliOptions.outFileName)
    , cqFile(cliOptions.inFileName, "r")
    , qualFile(cliOptions.outFileName, "w")
    , qualDecoder(cqFile, qualFile, cliOptions)
    , samParser(cliOptions.samFileName)
{
    // empty
}

CalqDecoder::~CalqDecoder(void)
{
    // empty
}

void CalqDecoder::decode(void)
{
    auto startTime = std::chrono::steady_clock::now();

    size_t compressedSize = 0;
    size_t fileHeaderSize = readFileHeader();

    size_t numBlocks = 0;
    while (cqFile.tell() < cqFile.size()) {
        //std::cout << ME << "Decoding block " << numBlocks << std::endl;
        compressedSize += qualDecoder.decodeBlock();

        // Get numRecordsInBlock from SAM file
//         do {
//             if (samRecordIsMapped(samParser.curr) == true) {
// 
//             } else {
// 
//             }
//         } while (samParser.next());

        numBlocks++;
    }

    // Print summary
    auto stopTime = std::chrono::steady_clock::now();
    auto diffTime = stopTime - startTime;
    auto diffTimeMs = std::chrono::duration_cast<std::chrono::milliseconds>(diffTime).count();
    auto diffTimeS = std::chrono::duration_cast<std::chrono::seconds>(diffTime).count();
    std::cout << ME << "DECOMPRESSION STATISTICS" << std::endl;
    std::cout << ME << "  Took " << diffTimeMs << " ms" << " ~= " << diffTimeS << " s" << std::endl;
    std::cout << ME << "  Decoded " << numBlocks << " block(s)" << std::endl;
    std::cout << ME << "  Compressed size: " << compressedSize << " (+ file header size: " << fileHeaderSize << ")" << std::endl;
    std::cout << ME << "  Speed (compressed size/time): " << ((double)(compressedSize/MB))/(double)((double)diffTimeMs/1000) << " MB/s" << std::endl;
}

size_t CalqDecoder::readFileHeader(void)
{
    std::cout << ME << "Reading file header" << std::endl;
    size_t ret = 0;

    const size_t magicSize = 5;
    char magic[magicSize];
    const size_t versionSize = 6;
    char version[versionSize];
    unsigned int polyploidy;
    ret += cqFile.read(magic, magicSize);
    ret += cqFile.read(version, versionSize);
    ret += cqFile.readUint32(&polyploidy);

    if (strncmp(version, VERSION, 5) != 0) {
        std::cout << ME << "Program version: " << VERSION << std::endl;
        std::cout << ME << "File version: " << version << std::endl;
        throwErrorException("CQ file was compressed with another version");
    }

    std::cout << ME << "  Magic: " << magic << std::endl;
    std::cout << ME << "  Version: " << version << std::endl;
    std::cout << ME << "  Polyploidy: " << polyploidy << std::endl;

    return ret;
}


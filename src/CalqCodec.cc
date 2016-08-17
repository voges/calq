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
#include "Common/debug.h"
#include <chrono>
#include <iostream>
#include <string.h>

static bool samRecordIsMapped(SAMRecord &samRecord)
{
    if (   ((samRecord.flag & 0x4) == 1)
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
                     const std::string &outFileName,
                     const std::vector<std::string> &fastaFileNames)
    : fastaReferences()
    , inFileName(inFileName)
    , outFileName(outFileName)
{
    // Get reference sequences
    FASTAParser fastaParser;
    for (auto const &fastaFileName : fastaFileNames) {
        std::cout << ME << "Parsing FASTA file: " << fastaFileName << std::endl;
        fastaParser.parseFile(fastaFileName, fastaReferences);
        std::cout << ME << "OK" << std::endl;
    }

    // Print info about found reference sequences
    std::cout << ME << "Found the following reference sequence(s):" << std::endl;
    for (auto const &fastaReference : fastaReferences) {
        std::cout << ME << "  " << fastaReference.header << " (sequence length: " << fastaReference.sequence.size() << ")" << std::endl;
    }
}

CalqCodec::~CalqCodec(void)
{
    // empty
}

CalqEncoder::CalqEncoder(const std::string &samFileName,
                         const std::string &cqFileName,
                         const std::vector<std::string> &fastaFileNames,
                         const unsigned int &blockSize,
                         const unsigned int &polyploidy,
                         const int &qvMin,
                         const int &qvMax)
    : CalqCodec(samFileName, cqFileName, fastaFileNames)
    , blockSize(blockSize)
    , cqFile(cqFileName, "w")
    , polyploidy(polyploidy)
    , qualEncoder(cqFile, fastaReferences, polyploidy, qvMin, qvMax)
    , samParser(samFileName)
{
    // empty
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
    std::cout << ME << "STATISTICS" << std::endl;
    std::cout << ME << "  Took " << std::chrono::duration_cast<std::chrono::milliseconds>(diffTime).count() << " ms" << " ~= " << std::chrono::duration_cast<std::chrono::seconds>(diffTime).count() << " s" << std::endl;
    std::cout << ME << "  Compressed " << numRecords << " record(s) (" << numMappedRecords << " mapped record(s) and " << numUnmappedRecords << " unmapped record(s)) in " << numBlocks << " block(s)" << std::endl;
    std::cout << ME << "  Uncompressed size: " << uncompressedSize << std::endl;
    std::cout << ME << "  Compressed size: " << compressedSize << " (+ file header size: " << fileHeaderSize << ")" << std::endl;
    std::cout << ME << "  Compression ratio: " << (double)compressedSize*100/(double)uncompressedSize << "%" << std::endl;
    std::cout << ME << "  Compression factor: " << (double)uncompressedSize/(double)compressedSize << std::endl;
    std::cout << ME << "  Bits per quality value: " << ((double)compressedSize * 8)/(double)uncompressedSize << std::endl;
    std::cout << ME << "  Speed (uncompressed size/time): " << ((double)(uncompressedSize/MB))/(double)std::chrono::duration_cast<std::chrono::seconds>(diffTime).count() << " MB/s" << std::endl;
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

CalqDecoder::CalqDecoder(const std::string &cqFileName,
                         const std::string &qualFileName,
                         const std::string &samFileName,
                         const std::vector<std::string> &fastaFileNames)
    : CalqCodec(cqFileName, qualFileName, fastaFileNames)
    , cqFile(cqFileName, "r")
    , qualFile(qualFileName, "w")
    , qualDecoder(cqFile, qualFile, fastaReferences)
    , samParser(samFileName)
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
    /*size_t fileHeaderSize = */readFileHeader();

    size_t numBlocks = 0;
    while (cqFile.tell() < cqFile.size()) {
        std::cout << ME << "Decoding block " << numBlocks << std::endl;
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
    std::cout << ME << "STATISTICS" << std::endl;
    std::cout << ME << "  Took " << std::chrono::duration_cast<std::chrono::milliseconds>(diffTime).count() << " ms"
              << " ~= " << std::chrono::duration_cast<std::chrono::seconds>(diffTime).count() << " s" << std::endl;
    std::cout << ME << "  Decoded " << numBlocks << " block(s)" << std::endl;
    std::cout << ME << "  Speed (compressed size/time): " << ((double)(compressedSize/MB))/(double)std::chrono::duration_cast<std::chrono::seconds>(diffTime).count() << " MB/s" << std::endl;
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
        throwErrorException("CQ file was compressed with another version");
    }

    std::cout << ME << "Magic: " << magic << std::endl;
    std::cout << ME << "Version: " << version << std::endl;
    std::cout << ME << "Polyploidy: " << polyploidy << std::endl;

    return ret;
}


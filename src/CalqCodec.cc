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
        std::cout << "  " << fastaReference.header;
        std::cout << " (sequence length: " << fastaReference.sequence.size();
        std::cout << ")" << std::endl;
    }
}

CalqCodec::~CalqCodec(void)
{
    // empty
}

CalqEncoder::CalqEncoder(const std::string &samFileName,
                         const std::string &cqFileName,
                         const std::vector<std::string> &fastaFileNames,
                         const unsigned int &polyploidy)
    : CalqCodec(samFileName, cqFileName, fastaFileNames)
    , cqFile(cqFileName, "w")
    , polyploidy(polyploidy)
    , qualEncoder(cqFile, fastaReferences, polyploidy)
    , samParser(samFileName)
{
    writeFileHeader();
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
    size_t numRecords = 0;

    // Send all records to the QualEncoder
    std::string rnamePrev("");
    uint32_t posPrev = 0;
    bool first = true;
    qualEncoder.startBlock();

    do {
        uncompressedSize += strlen(samParser.curr.qual);
        numRecords++;
        if (samRecordIsMapped(samParser.curr) == true) {
            if (first == false) {
                if (samParser.curr.rname != rnamePrev) {
                    // Start a new block
                    compressedSize += qualEncoder.finishBlock();
                    qualEncoder.startBlock();
                } else {
                    if (samParser.curr.pos < posPrev) {
                        throwErrorException("SAM file is not sorted");
                    }
                }
            }
            posPrev = samParser.curr.pos;
            rnamePrev = samParser.curr.rname;
            qualEncoder.addMappedRecordToBlock(samParser.curr);
            first = false;
        } else {
            qualEncoder.addUnmappedRecordToBlock(samParser.curr);
        }
    } while (samParser.next());

    compressedSize += qualEncoder.finishBlock();

    // Print summary
    auto stopTime = std::chrono::steady_clock::now();
    auto diffTime = stopTime - startTime;
    std::cout << ME << "Took " << std::chrono::duration_cast<std::chrono::milliseconds>(diffTime).count() << " ms"
              << " ~= " << std::chrono::duration_cast<std::chrono::seconds>(diffTime).count() << " s" << std::endl;
    std::cout << ME << "Compressed " << numRecords << " record(s)" << std::endl;
    std::cout << ME << "Uncompressed size: " << uncompressedSize << std::endl;
    std::cout << ME << "Compressed size: " << compressedSize << std::endl;
}

void CalqEncoder::writeFileHeader(void)
{
    size_t fileHeaderSize = 0;

    const size_t magicSize = 5;
    char magic[magicSize] = "calq";
    const size_t versionSize = 6;
    char version[versionSize] = VERSION;
    fileHeaderSize += cqFile.write(magic, magicSize);
    fileHeaderSize += cqFile.write(version, versionSize);
    fileHeaderSize += cqFile.writeUint32(polyploidy);

    std::cout << ME << "File header size: " << fileHeaderSize << std::endl;
}

CalqDecoder::CalqDecoder(const std::string &cqFileName,
                         const std::string &qualFileName,
                         const std::vector<std::string> &fastaFileNames)
    : CalqCodec(cqFileName, qualFileName, fastaFileNames)
    , cqFile(cqFileName, "r")
    , qualFile(qualFileName, "w")
    , qualDecoder(cqFile, qualFile, fastaReferences)
{
    readFileHeader();
}

CalqDecoder::~CalqDecoder(void)
{
    // empty
}

void CalqDecoder::decode(void)
{
    size_t numBlocks = 0;
    qualDecoder.decodeBlock();
    std::cout << ME << "Decoded " << numBlocks << " block(s)" << std::endl;
}

void CalqDecoder::readFileHeader(void)
{
    std::cout << ME << "Reading file header" << std::endl;

    const size_t magicSize = 5;
    char magic[magicSize];
    const size_t versionSize = 6;
    char version[versionSize];
    unsigned int polyploidy;
    cqFile.read(magic, magicSize);
    cqFile.read(version, versionSize);
    cqFile.readUint32(&polyploidy);

    if (strncmp(version, VERSION, 5) != 0) {
        throwErrorException("CQ file was compressed with another version");
    }

    std::cout << ME << "Magic: " << magic << std::endl;
    std::cout << ME << "Version: " << version << std::endl;
    std::cout << ME << "Polyploidy: " << polyploidy << std::endl;
}


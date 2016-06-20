/** @file CalqCodec.cc
 *  @brief This file contains the implementations of the CalqEncoder and
 *         CalqDecoder class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: what (who)
 */

#include "Codecs/CalqCodec.h"
#include "cmake_config.h"
#include "definitions.h"
#include "Exceptions.h"
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
    // get reference sequences
    FASTAParser fastaParser;
    for (auto const &fastaFileName : fastaFileNames) {
        std::cout << "Parsing FASTA file: " << fastaFileName << std::endl;
        fastaParser.parseFile(fastaFileName, fastaReferences);
        std::cout << "OK" << std::endl;
    }

    // print info about found reference sequences
    std::cout << "Found the following reference sequence(s):" << std::endl;
    for (auto const &fastaReference : fastaReferences) {
        std::cout << fastaReference.header << " (sequence length: " << fastaReference.sequence.size() << ")" << std::endl;
    }
}

CalqCodec::~CalqCodec(void)
{
    // empty
}

CalqEncoder::CalqEncoder(const std::string &samFileName,
                         const std::string &cqFileName,
                         const std::vector<std::string> &fastaFileNames)
    : CalqCodec(samFileName, cqFileName, fastaFileNames)
    , ofbs()
    , qualEncoder(ofbs, fastaReferences)
    , samParser(samFileName)
{
    // associate the outfile bitstream with the CQ file
    ofbs.open(cqFileName);
}

CalqEncoder::~CalqEncoder(void)
{
    ofbs.close();
}

void CalqEncoder::encode(void)
{
    auto startTime = std::chrono::steady_clock::now();

    size_t uncompressedSize = 0;
    size_t compressedSize = 0;
    size_t numRecords = 0;

    // send all records to the QualEncoder
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
                    // start a new block
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

    // print summary
    auto stopTime = std::chrono::steady_clock::now();
    auto diffTime = stopTime - startTime;
    std::cout << "Took " << std::chrono::duration_cast<std::chrono::milliseconds>(diffTime).count() << " ms"
              << " ~= " << std::chrono::duration_cast<std::chrono::seconds>(diffTime).count() << " s" << std::endl;
    std::cout << "Compressed " << numRecords << " record(s)" << std::endl;
    std::cout << "Uncompressed size: " << uncompressedSize << std::endl;
    std::cout << "Compressed size: " << compressedSize << std::endl;
}

CalqDecoder::CalqDecoder(const std::string &cqFileName,
                         const std::string &qualFileName,
                         const std::vector<std::string> &fastaFileNames)
    : CalqCodec(cqFileName, qualFileName, fastaFileNames)
    , ifbs()
    , ofs()
    , qualDecoder(ifbs, ofs, fastaReferences)
{
    // associate the infile bitstream with the CQ file and the outfile stream
    // with the QUAL file
    ifbs.open(cqFileName);
    ofs.open(qualFileName);
}

CalqDecoder::~CalqDecoder(void)
{
    ifbs.close();
    ofs.close();
}

void CalqDecoder::decode(void)
{
    std::cout << "CalqDecoder::decode() not yet implemented!" << std::endl;
}


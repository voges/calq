/** @file CalqCodec.cc
 *  @brief This file contains the implementations of the CalqEncoder,
 *         CalqDecoder, and CalqInfoTool classes.
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

CalqCodec::CalqCodec(const std::string &inFileName,
                     const std::string &outFileName,
                     const std::vector<std::string> &fastaFileNames)
    : fastaFileNames(fastaFileNames)
    , fastaParser()
    , fastaReferences()
    , inFileName(inFileName)
    , outFileName(outFileName)
{
    readFastaReferences();
}

CalqCodec::~CalqCodec(void)
{
    // empty
}

void CalqCodec::readFastaReferences(void)
{
    // get reference sequences
    for (auto const &fastaFileName : fastaFileNames) {
        std::cout << "Parsing FASTA file: " << fastaFileName << std::endl;
        fastaParser.parseFile(fastaFileName, fastaReferences);
    }

    // print info about found reference sequences
    std::cout << "Found the following reference sequences: " << std::endl;
    for (auto const &fastaReference : fastaReferences) {
        std::cout << fastaReference.header << " (sequence length: " << fastaReference.sequence.size() << ")" << std::endl;
    }
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
    qualEncoder.startBlock();
    SAMRecord *currentRecord = &(samParser.curr);
    do {
        uncompressedSize += strlen(samParser.curr.qual);
        qualEncoder.addRecordToBlock(samParser.curr);
        numRecords++;
    } while (samParser.next());
    compressedSize = qualEncoder.finishBlock();

    // print summary
    auto stopTime = std::chrono::steady_clock::now();
    auto diffTime = stopTime - startTime;
    std::cout << "Took " << std::chrono::duration_cast<std::chrono::seconds>(diffTime).count() << " s" << std::endl;
    std::cout << "Compressed " << numRecords << " record(s)" << std::endl;
    std::cout << "Uncompressed size: " << uncompressedSize << std::endl;
    std::cout << "Compressed size: " << compressedSize << std::endl;
    std::cout << "Compression Ratio (CR): " << (double)uncompressedSize/(double)compressedSize*100 << std::endl;
    std::cout << "Compression Factor (CF): " << (double)compressedSize/(double)uncompressedSize*100 << std::endl;
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


/** @file CalqCodec.cc
 *  @brief This file contains the implementations of the CalqEncoder,
 *         CalqDecoder, and CalqInfoTool classes.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#include "CalqCodec.h"
#include "cmake_config.h"
#include "definitions.h"
#include "Exceptions.h"
#include <chrono>
#include <iostream>
#include <string.h>

CalqEncoder::CalqEncoder(const std::string &samFileName,
                         const std::string &cqFileName,
                         const std::vector<std::string> &fastaFileNames)
    : fastaParser()
    , samFileName()
    , ofbs()
    , cqFileName(cqFileName)
    , qualEncoder(ofbs)
    , fastaFileNames(fastaFileNames)
    , samParser(samFileName)
{
    ofbs.open(cqFileName);

    // get references sequences
    for (auto const &fastaFileName : fastaFileNames) {
        std::cout << "Parsing FASTA file: " << fastaFileName << std::endl;
        fastaParser.parseFile(fastaFileName, fastaReferences);
    }

    // print info about found reference sequences
    std::cout << "Found the following reference sequences: " << std::endl;
    for (auto const &fastaReference : fastaReferences) {
        std::cout << "header: " << fastaReference.header << std::endl;
        std::cout << "sequence: " << fastaReference.sequence << std::endl;
    }
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

    // encode all records
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

CalqDecoder::CalqDecoder(const std::string &inFileName,
                         const std::string &outfileName,
                         const std::vector<std::string> &referenceFileNames)
    //: fastaParser(referenceFileName)
    : ifbs()
    , inFileName(inFileName)
    , ofs()
    , outFileName(outFileName)
    , qualDecoder(ifbs, ofs)
    , referenceFileNames(referenceFileNames)
{
    ifbs.open(inFileName);
    ofs.open(outFileName);
}

CalqDecoder::~CalqDecoder(void)
{
    ifbs.close();
    ofs.close();
}

void CalqDecoder::decode(void)
{

}


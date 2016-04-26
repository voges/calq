#include "CalqCodec.h"
#include "config.h"
#include "Exceptions.h"
#include <chrono>
#include <iostream>
#include <string.h>

CalqEncoder::CalqEncoder(std::string &infileName, std::string &outfileName, size_t &blockSize)
    : infileName(infileName)
    , outfileName(outfileName)
    , samParser(infileName)
    , ofbs()
    , blockSize(blockSize)
    , qualEncoder(ofbs)
{
    ofbs.open(outfileName);
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
    size_t numBlocks = 0;
    size_t numRecords = 0;

    // write file header
    const unsigned char magic[5] = "calq";
    const unsigned char versionMajor = VERSION_MAJOR + '0';
    const unsigned char versionMinor = VERSION_MINOR + '0';
    const unsigned char versionPatch = VERSION_PATCH + '0';
    ofbs << magic << versionMajor << versionMinor << versionPatch;
    //ofbs.seekg(2*sizeof(size_t), std::ios::cur); // space for numBlocks and numRecords

    size_t numRecordsInBlock = 0;
    SAMRecord *currentRecord = &(samParser.curr);

    while (true) {
        bool haveNextRecord = samParser.hasNext();

        if (numRecordsInBlock >= blockSize || haveNextRecord == false) {
            compressedSize += qualEncoder.finishBlock();
            numRecordsInBlock = 0;
            numBlocks++;
            if (haveNextRecord == false) {
                break;
            }
        } else {
            samParser.next();
            uncompressedSize += strlen(currentRecord->qual);
            qualEncoder.encodeRecord(*currentRecord);
            numRecords++;
            numRecordsInBlock++;
        }
    }

    // print summary
    auto stopTime = std::chrono::steady_clock::now();
    auto diffTime = stopTime - startTime;
    std::cout << "Took " << std::chrono::duration_cast<std::chrono::seconds>(diffTime).count() << " s" << std::endl;
    std::cout << "Compressed " << numRecords << " record(s)" << std::endl;
    std::cout << "Wrote " << numBlocks << " block(s)" << std::endl;
    std::cout << "Uncompressed size: " << uncompressedSize << std::endl;
    std::cout << "Compressed size: " << compressedSize << std::endl;
    std::cout << "Compression Ratio (CR): " << (double)uncompressedSize/(double)compressedSize*100 << std::endl;
    std::cout << "Compression Factor (CF): " << (double)compressedSize/(double)uncompressedSize*100 << std::endl;
}


CalqDecoder::CalqDecoder(std::string &infileName, std::string &outfileName)
    : infileName(infileName)
    , outfileName(outfileName)
{
    // Empty
}

CalqDecoder::~CalqDecoder(void)
{
    // Empty
}

void CalqDecoder::decode(void)
{
    throwErrorException("Not yet implemented");
}

CalqInfoTool::CalqInfoTool(std::string& infileName): infileName(infileName)
{
    // Empty
}

CalqInfoTool::~CalqInfoTool(void)
{
    // Empty
}

void CalqInfoTool::extractInfo(void)
{
    throwErrorException("Not yet implemented");
}


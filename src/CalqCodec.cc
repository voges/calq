/** @file CalqCodec.cc
 *  @brief This file contains the implementations of the CalqEncoder,
 *         CalqDecoder, and CalqInfoTool classes.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#include "CalqCodec.h"
#include "config.h"
#include "definitions.h"
#include "Exceptions.h"
#include <chrono>
#include <iostream>
#include <string.h>

static const int BLOCK_SIZE = 10000;

CalqEncoder::CalqEncoder(const std::string &inFileName,
                         const std::string &outFileName,
                         const std::string &referenceFileName)
    : blockPositionList()
    , blockSize(BLOCK_SIZE)
    //, fastaParser(referenceFileName)
    , inFileName(inFileName)
    , ofbs()
    , outFileName(outFileName)
    , qualEncoder(ofbs)
    , referenceFileName(referenceFileName)
    , samParser(inFileName)
{
    ofbs.open(outFileName);
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
    size_t numBlocks = 0;  // total number of blocks generated
    size_t numRecords = 0; // total number of records encoded

    // write file header and reserve 4*8 bytes for blockSize, numBlocks,
    // numRecords, blockListPos
    const BYTE fileHeaderMagicOne = 'C';
    const BYTE fileHeaderMagicTwo = 'Q';
    const BYTE fileHeaderVersionMajor = VERSION_MAJOR + '0';
    const BYTE fileHeaderVersionMinor = VERSION_MINOR + '0';
    const BYTE fileHeaderVersionPatch = VERSION_PATCH + '0';
    ofbs.writeByte(fileHeaderMagicOne);
    ofbs.writeByte(fileHeaderMagicTwo);
    ofbs.writeByte(fileHeaderVersionMajor);
    ofbs.writeByte(fileHeaderVersionMinor);
    ofbs.writeByte(fileHeaderVersionPatch);
    std::streampos fileHeaderPos = ofbs.tellp();
    ofbs.seekp(4*sizeof(uint64_t), std::ios::cur);

    size_t numRecordsInBlock = 0;
    SAMRecord *currentRecord = &(samParser.curr);

    while (true) {
        bool haveNextRecord = samParser.hasNext();

        if (numRecordsInBlock >= blockSize || haveNextRecord == false) {
            compressedSize += qualEncoder.finishBlock();
            //blockPositionList.push_back(ofbs.tellp());
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

    // complete file header
    //std::streampos blockPositionListPos = ofbs.tellp();
    ofbs.seekp(fileHeaderPos, std::ios::beg);
    ofbs.writeUint64(blockSize);
    ofbs.writeUint64(numBlocks);
    ofbs.writeUint64(numRecords);
    //ofbs.writeUint64(blockPositionListPos);
    //ofbs.seekp(blockPositionListPos, std::ios::beg);

    // write block position list
    //for(auto const &blockPosition: blockPositionList) {
    //    ofbs.writeUint64((uint64_t)blockPosition);
    //    std::cout << blockPosition << std::endl;
    //}

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

CalqDecoder::CalqDecoder(const std::string &inFileName,
                         const std::string &outfileName,
                         const std::string &referenceFileName)
    //: fastaParser(referenceFileName)
    : ifbs()
    , inFileName(inFileName)
    , ofs()
    , outFileName(outFileName)
    , qualDecoder(ifbs)
    , referenceFileName(referenceFileName)
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
    // read file header
    BYTE fileHeaderMagicOne = 0x00;
    BYTE fileHeaderMagicTwo = 0x00;
    BYTE fileHeaderVersionMajor = 0x00;
    BYTE fileHeaderVersionMinor = 0x00;
    BYTE fileHeaderVersionPatch = 0x00;
    uint64_t fileHeaderBlockSize = 0ULL;
    uint64_t fileHeaderNumBlocks = 0ULL;
    uint64_t fileHeaderNumRecords = 0ULL;
    //uint64_t fileHeaderBlockPositionListPos = 0ULL;

    ifbs.readByte(fileHeaderMagicOne);
    ifbs.readByte(fileHeaderMagicTwo);
    ifbs.readByte(fileHeaderVersionMajor);
    ifbs.readByte(fileHeaderVersionMinor);
    ifbs.readByte(fileHeaderVersionPatch);
    ifbs.readUint64(fileHeaderBlockSize);
    ifbs.readUint64(fileHeaderNumBlocks);
    ifbs.readUint64(fileHeaderNumRecords);
    //ifbs.readUint64(fileHeaderBlockPositionListPos);

    // print file header
    std::cout << "Magic: " << fileHeaderMagicOne << fileHeaderMagicTwo << std::endl;
    std::cout << "File version: ";
    std::cout << fileHeaderVersionMajor << ".";
    std::cout << fileHeaderVersionMinor << ".";
    std::cout << fileHeaderVersionPatch << std::endl;
    std::cout << "Block size: " << fileHeaderBlockSize << std::endl;
    std::cout << "#blocks: " << fileHeaderNumBlocks << std::endl;
    std::cout << "#records: " << fileHeaderNumRecords << std::endl;
}


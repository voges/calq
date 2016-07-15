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
//#include "cmake_config.h"
#include "Common/Exceptions.h"
#include "Common/debug.h"
#include "Common/definitions.h"
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
        std::cout << ME << "Parsing FASTA file: " << fastaFileName << std::endl;
        fastaParser.parseFile(fastaFileName, fastaReferences);
        std::cout << ME << "OK" << std::endl;
    }

    // print info about found reference sequences
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
                         const int &polyploidy)
    : CalqCodec(samFileName, cqFileName, fastaFileNames)
    , ofbs()
    , qualEncoder(ofbs, fastaReferences, polyploidy)
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
    std::cout << ME << "Took " << std::chrono::duration_cast<std::chrono::milliseconds>(diffTime).count() << " ms"
              << " ~= " << std::chrono::duration_cast<std::chrono::seconds>(diffTime).count() << " s" << std::endl;
    std::cout << ME << "Compressed " << numRecords << " record(s)" << std::endl;
    std::cout << ME << "Uncompressed size: " << uncompressedSize << std::endl;
    std::cout << ME << "Compressed size: " << compressedSize << std::endl;
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
    size_t numBlocks = 0;

    while (true) {
        /* The following ifbs.get() is needed, because the decoder never cheks for eof(). He rather checks, if he
         * read the exact amount of compressed bits + header for every block. Thus the ifbs never gets a failed read
         * attempt, while the decoder is active. After the last decoded block, the ifbs streams position is at the exact
         * end of the file. The eofbit is only set, with a failed get()-attempt. Because this attempt never happens within the 
         * decoder itself, we need to call it upfront and make it a leaving condition for this loop.
         */
        std::streampos cur = ifbs.tellg();
        ifbs.get();
        if(ifbs.eof()){
            break;
        }
        ifbs.seekg(cur);
        //end of fix
        
        
        qualDecoder.decodeBlock();
        numBlocks++;
    }

    std::cout << ME << "Decoded " << numBlocks << " block(s)" << std::endl;
}


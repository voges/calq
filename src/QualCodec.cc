/** @file QualCodec.cc
 *  @brief This file contains the implementations of the QualEncoder and
 *         QualDecoder classes.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "QualCodec.h"
#include "Common/CLIOptions.h"
#include "Common/constants.h"
#include "Common/Exceptions.h"
#include "Common/debug.h"
#include "Compressors/range/range.h"
#include "Compressors/rle/rle.h"
#include <chrono>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdlib.h>
#include <string.h>

static const unsigned int QUANTIZER_STEP_MIN = 2;
static const unsigned int QUANTIZER_STEP_MAX = 8;
static const unsigned int QUANTIZER_NUM = QUANTIZER_STEP_MAX-QUANTIZER_STEP_MIN+1; // 2-8 steps
static const unsigned int QUANTIZER_IDX_MIN = 0;
static const unsigned int QUANTIZER_IDX_MAX = QUANTIZER_NUM-1;

QualEncoder::QualEncoder(File &cqFile, const CLIOptions &cliOptions)
    // Class scope
    : fastaReferences()
    , quantizedPrintout(false)
    , verbose(false)
    , qvMin(cliOptions.qvMin)
    , qvMax(cliOptions.qvMax)
    , cqFile(cqFile)
    , genotyper(cliOptions.polyploidy, QUANTIZER_NUM, QUANTIZER_IDX_MIN, QUANTIZER_IDX_MAX, qvMin, qvMax)
    , uniformQuantizers()
    , uncompressedSize(0)
    , compressedSize(0)
    , numBlocks(0)
    , numRecords(0)
    , numMappedRecords(0)
    , numUnmappedRecords(0)
    // Block scope
    , uncompressedSizeOfBlock(0)
    , compressedSizeOfBlock(0)
    , numRecordsInBlock(0)
    , numMappedRecordsInBlock(0)
    , numUnmappedRecordsInBlock(0)
    , qi("")
    , qvi("")
    , uqv("")
    , referenceName("")
    , reference("")
    , referencePosMin(std::numeric_limits<uint32_t>::max())
    , referencePosMax(std::numeric_limits<uint32_t>::min())
    , mappedRecordQueue()
    , observedNucleotides()
    , observedQualityValues()
    , observedPosMin(std::numeric_limits<uint32_t>::max())
    , observedPosMax(std::numeric_limits<uint32_t>::min())
    , quantizerIndices()
    , quantizerIndicesPosMin(std::numeric_limits<uint32_t>::max())
    , quantizerIndicesPosMax(std::numeric_limits<uint32_t>::min())
{
    if (cliOptions.polyploidy == 0) {
        std::cout << "Polyploidy: " << cliOptions.polyploidy << std::endl;
        throwErrorException("Polyploidy must be greater than zero");
    }

    if (cliOptions.quantizedPrintout == true) {
        std::string qualFileName(cliOptions.outFileName);
        qualFileName += ".qual";
        std::cout << ME << "Printing quantized quality values to file: " << qualFileName << std::endl;
        if ((fileExists(qualFileName) == true) && (cliOptions.force == false)) {
            throwErrorException("Quantized quality values file already exists (use option 'f' to force overwriting)");
        }

        qualFile.open(qualFileName);
        quantizedPrintout = true;
    }

    std::cerr << "rname,pos,depth,entropy,k" << std::endl;

    if (cliOptions.verbose == true) {
        verbose = true;
    }

    // Init uniform quantizers
    std::cout << ME << "Initializing " << QUANTIZER_NUM << " uniform quantizers" << std::endl;
    for (unsigned int i = 0; i < QUANTIZER_NUM; i++) {
        unsigned int numberOfSteps = QUANTIZER_STEP_MIN;
        UniformQuantizer uniformQuantizer(qvMin, qvMax, numberOfSteps);
        //uniformQuantizer.print();
        uniformQuantizers.insert(std::pair<int,UniformQuantizer>(i, uniformQuantizer));
    }
}

QualEncoder::~QualEncoder(void)
{
    if (qualFile.is_open()) {
        qualFile.close();
    }

    if (verbose == true) {
        std::cout << ME << "QualEncoder statistics:" << std::endl;
        std::cout << ME << "  Uncompressed size:  " << std::setw(9) << uncompressedSize << std::endl;
        std::cout << ME << "  Compressed size:    " << std::setw(9) << compressedSize << std::endl;
        std::cout << ME << "  Blocks:             " << std::setw(9) << numBlocks << std::endl;
        std::cout << ME << "  Records:            " << std::setw(9) << numRecords << std::endl;
        std::cout << ME << "    Mapped:           " << std::setw(9) << numMappedRecords << std::endl;
        std::cout << ME << "    Unmapped:         " << std::setw(9) << numUnmappedRecords << std::endl;
        std::cout << ME << "  Compression ratio:  " << std::setw(9) << (double)compressedSize*100/(double)uncompressedSize << "%" << std::endl;
        std::cout << ME << "  Compression factor: " << std::setw(9) << (double)uncompressedSize/(double)compressedSize << std::endl;

    }
}

void QualEncoder::startBlock(void)
{
    std::cout << ME << "Starting block " << numBlocks << std::endl;

    // Reset all variables in block scope
    uncompressedSizeOfBlock = 0;
    compressedSizeOfBlock = 0;
    numRecordsInBlock = 0;
    numMappedRecordsInBlock = 0;
    numUnmappedRecordsInBlock = 0;
    qi = "";
    qvi = "";
    uqv = "";
    reference = "";
    referencePosMin = std::numeric_limits<uint32_t>::max();
    referencePosMax = std::numeric_limits<uint32_t>::min();
    mappedRecordQueue = std::queue<MappedRecord>();
    observedNucleotides.clear();
    observedQualityValues.clear();
    observedPosMin = std::numeric_limits<uint32_t>::max();
    observedPosMax = std::numeric_limits<uint32_t>::min();
    quantizerIndices.clear();
    quantizerIndicesPosMin = std::numeric_limits<uint32_t>::max();
    quantizerIndicesPosMax = std::numeric_limits<uint32_t>::min(); 
}

void QualEncoder::addUnmappedRecordToBlock(const SAMRecord &samRecord)
{
    uncompressedSizeOfBlock += strlen(samRecord.qual);
    uncompressedSize += strlen(samRecord.qual);

    encodeUnmappedQualityValues(std::string(samRecord.qual));

    numUnmappedRecordsInBlock++;
    numRecordsInBlock++;
    numUnmappedRecords++;
    numRecords++;
}

void QualEncoder::addMappedRecordToBlock(const SAMRecord &samRecord)
{
    uncompressedSizeOfBlock += strlen(samRecord.qual);
    uncompressedSize += strlen(samRecord.qual);

    // If the reference is empty, then this is the first mapped record being
    // added to the current block. If so, we load the corresponding reference
    // and initially set the min mapping positions for this block.
    if (reference.empty()) {
        loadFastaReference(std::string(samRecord.rname));
        observedPosMin = samRecord.pos - 1; // SAM format counts from 1
        quantizerIndicesPosMin = observedPosMin;
    }

    // Construct a mapped record from the given SAM record
    MappedRecord mappedRecord(samRecord);

    // Check if we need to update the max mapping position for records in this
    // block. If so, we also need to update the size of our observation vectors.
    if (mappedRecord.posMax > observedPosMax) {
        observedPosMax = mappedRecord.posMax;
        size_t observedSize = observedPosMax - observedPosMin + 1;
        observedNucleotides.resize(observedSize);
        observedQualityValues.resize(observedSize);
    }

    // Parse CIGAR string and assign nucleotides and QVs from the record 
    // (we name them 'observations') to their positions within the reference;
    // then push the current record to the queue.
    mappedRecord.extractObservations(observedPosMin, observedNucleotides, observedQualityValues);
    mappedRecordQueue.push(mappedRecord);

    // Check which genomic positions are ready for processing. Then compute
    // the quantizer indices for these positions and shrink the observation
    // vectors
    while (observedPosMin < mappedRecord.posMin) {
        std::cerr << referenceName << "," << observedPosMin << ",";
        int k = genotyper.computeQuantizerIndex(observedNucleotides[0], observedQualityValues[0]);
        std::cerr << k << std::endl;

        quantizerIndices.push_back(k);
        quantizerIndicesPosMax = observedPosMin;
        observedNucleotides.erase(observedNucleotides.begin());
        observedQualityValues.erase(observedQualityValues.begin());
        observedPosMin++;
    }

    // Check which records from the queue are ready for processing
    while (mappedRecordQueue.front().posMax < observedPosMin) {
        encodeMappedQualityValues(mappedRecordQueue.front());
        uint32_t currFirstPosition = mappedRecordQueue.front().posMin;
        mappedRecordQueue.pop();

        // Check if we can shrink the quantizer indices vector
        if (mappedRecordQueue.empty() == false) {
            uint32_t nextFirstPosition = mappedRecordQueue.front().posMin;
            if (nextFirstPosition > currFirstPosition) {
                quantizerIndices.erase(quantizerIndices.begin(), quantizerIndices.begin()+nextFirstPosition-currFirstPosition);
                quantizerIndicesPosMin += nextFirstPosition-currFirstPosition;
            }
        }
    }

    numMappedRecordsInBlock++;
    numRecordsInBlock++;
    numMappedRecords++;
    numRecords++;
}

size_t QualEncoder::finishBlock(void)
{
    size_t ret = 0;

    // Compute all remaining quantizers
    while (observedPosMin <= observedPosMax) {
        std::cerr << referenceName << "," << observedPosMin << ",";
        int k = genotyper.computeQuantizerIndex(observedNucleotides[0], observedQualityValues[0]);
        std::cerr << k << std::endl;

        quantizerIndices.push_back(k);
        quantizerIndicesPosMax = observedPosMin;
        observedNucleotides.erase(observedNucleotides.begin());
        observedQualityValues.erase(observedQualityValues.begin());
        observedPosMin++;
    }

    // Process all remaining records from queue
    while (mappedRecordQueue.empty() == false) {
        encodeMappedQualityValues(mappedRecordQueue.front());
        uint32_t currFirstPosition = mappedRecordQueue.front().posMin;
        mappedRecordQueue.pop();

        // Check if we can shrink the quantizerIndices vector
        if (mappedRecordQueue.empty() == false) {
            uint32_t nextFirstPosition = mappedRecordQueue.front().posMin;
            if (nextFirstPosition > currFirstPosition) {
                quantizerIndices.erase(quantizerIndices.begin(), quantizerIndices.begin()+nextFirstPosition-currFirstPosition);
                quantizerIndicesPosMin += nextFirstPosition-currFirstPosition;
            }
        }
    }

    // Write block number and numbers of records to output file
    ret += cqFile.writeUint64(numBlocks);
    ret += cqFile.writeUint64(numRecordsInBlock);
    ret += cqFile.writeUint64(numMappedRecordsInBlock);
    ret += cqFile.writeUint64(numUnmappedRecordsInBlock);

    // RLE and range encoding for quantizer indices
    //std::cout << "Quantizer indices: " << qi << std::endl;
    unsigned char *qiBuffer = (unsigned char *)qi.c_str();
    size_t qiBufferSize = qi.length();
    if (qiBufferSize > 0) {
        size_t qiRLESize = 0;
        unsigned char *qiRLE = rle_encode(qiBuffer, qiBufferSize, &qiRLESize, QUANTIZER_NUM, (unsigned char)'0');

        size_t numQiBlocks = (size_t)ceil((double)qiRLESize / (double)MB);
        ret += cqFile.writeUint64(numQiBlocks);

        size_t encodedQiBytes = 0;
        while (encodedQiBytes < qiRLESize) {
            unsigned int bytesToEncode = 0;
            if ((qiRLESize - encodedQiBytes) > MB) {
                bytesToEncode = MB;
            } else {
                bytesToEncode = qiRLESize - encodedQiBytes;
            }

            unsigned int qiRangeSize = 0;
            unsigned char *qiRange = range_compress_o1(qiRLE+encodedQiBytes, (unsigned int)bytesToEncode, &qiRangeSize);

            if (qiRangeSize >= bytesToEncode) {
                ret += cqFile.writeUint8(0);
                ret += cqFile.writeUint32(bytesToEncode);
                ret += cqFile.write(qiRLE+encodedQiBytes, bytesToEncode);
            } else {
                ret += cqFile.writeUint8(1);
                ret += cqFile.writeUint32(qiRangeSize);
                ret += cqFile.write(qiRange, qiRangeSize);
            }

            encodedQiBytes += bytesToEncode;
            free(qiRange);
        }
        free(qiRLE);
    } else {
        ret += cqFile.writeUint64(0);
    }

    // RLE and range encoding for quality value indices
    //std::cout << "Quality value indices: " << qvi << std::endl;
    unsigned char *qviBuffer = (unsigned char *)qvi.c_str();
    size_t qviBufferSize = qvi.length();
    if (qviBufferSize > 0) {
        size_t qviRLESize = 0;
        unsigned char *qviRLE = rle_encode(qviBuffer, qviBufferSize, &qviRLESize, QUANTIZER_STEP_MAX, (unsigned char)'0');

        size_t numQviBlocks = (size_t)ceil((double)qviRLESize / (double)MB);
        ret += cqFile.writeUint64(numQviBlocks);

        size_t encodedQviBytes = 0;
        while (encodedQviBytes < qviRLESize) {
            unsigned int bytesToEncode = 0;
            if ((qviRLESize - encodedQviBytes) > MB) {
                bytesToEncode = MB;
            } else {
                bytesToEncode = qviRLESize - encodedQviBytes;
            }

            unsigned int qviRangeSize = 0;
            unsigned char *qviRange = range_compress_o1(qviRLE+encodedQviBytes, (unsigned int)bytesToEncode, &qviRangeSize);

            if (qviRangeSize >= bytesToEncode) {
                ret += cqFile.writeUint8(0);
                ret += cqFile.writeUint32(bytesToEncode);
                ret += cqFile.write(qviRLE+encodedQviBytes, bytesToEncode);
            } else {
                ret += cqFile.writeUint8(1);
                ret += cqFile.writeUint32(qviRangeSize);
                ret += cqFile.write(qviRange, qviRangeSize);
            }

            encodedQviBytes += bytesToEncode;
            free(qviRange);
        }
        free(qviRLE);
    } else {
        ret += cqFile.writeUint64(0);
    }

    // Range encoding for unmapped quality values
    //std::cout << "Unmapped quality values: " << uqv << std::endl;
    unsigned char *uqvBuffer = (unsigned char *)uqv.c_str();
    size_t uqvBufferSize = uqv.length();
    if (uqvBufferSize > 0) {
        size_t numUqvBlocks = (size_t)ceil((double)uqvBufferSize / (double)MB);
        ret += cqFile.writeUint64(numUqvBlocks);

        size_t encodedUqvBytes = 0;
        while (encodedUqvBytes < uqvBufferSize) {
            unsigned int bytesToEncode = 0;
            if ((uqvBufferSize - encodedUqvBytes) > MB) {
                bytesToEncode = MB;
            } else {
                bytesToEncode = uqvBufferSize - encodedUqvBytes;
            }

            unsigned int uqvRangeSize = 0;
            unsigned char *uqvRange = range_compress_o1(uqvBuffer+encodedUqvBytes, (unsigned int)bytesToEncode, &uqvRangeSize);

            if (uqvRangeSize >= bytesToEncode) {
                ret += cqFile.writeUint8(0);
                ret += cqFile.writeUint32(bytesToEncode);
                ret += cqFile.write(uqvBuffer+encodedUqvBytes, bytesToEncode);
            } else {
                ret += cqFile.writeUint8(1);
                ret += cqFile.writeUint32(uqvRangeSize);
                ret += cqFile.write(uqvRange, uqvRangeSize);
            }
 
            encodedUqvBytes += bytesToEncode;
            free(uqvRange);
        }
    } else {
        ret += cqFile.writeUint64(0);
    }

    // Update sizes & counters
    compressedSizeOfBlock += ret;
    compressedSize += compressedSizeOfBlock;
    numBlocks++;

    // Print statistics for this block
    if (verbose == true) {
        std::cout << ME << "Encoded block " << (numBlocks-1) << std::endl;
        std::cout << ME << "  Uncompressed size:  " << std::setw(9) << uncompressedSizeOfBlock << std::endl;
        std::cout << ME << "  Compressed size:    " << std::setw(9) << compressedSizeOfBlock << std::endl;
        std::cout << ME << "  Records:            " << std::setw(9) << numRecordsInBlock << std::endl;
        std::cout << ME << "    Mapped:           " << std::setw(9) << numMappedRecordsInBlock << std::endl;
        std::cout << ME << "    Unmapped:         " << std::setw(9) << numUnmappedRecordsInBlock << std::endl;
        std::cout << ME << "  Compression ratio:  " << std::setw(9) << (double)compressedSizeOfBlock*100/(double)uncompressedSizeOfBlock << "%" << std::endl;
        std::cout << ME << "  Compression factor: " << std::setw(9) << (double)uncompressedSizeOfBlock/(double)compressedSizeOfBlock << std::endl;
    }

    std::cout << ME << "Finished block " << (numBlocks-1) << " (" << numRecordsInBlock << " records)" << std::endl;

    return ret;
}

void QualEncoder::loadFastaReference(const std::string &rname)
{
    // Find FASTA reference for this RNAME
    if (verbose) {
        std::cout << ME << "Searching FASTA reference for: " << rname << std::endl;
    }
    bool foundFastaReference = false;

    for (auto const &fastaReference : fastaReferences) {
        if (fastaReference.header == rname) {
            if (foundFastaReference == true) {
                throwErrorException("Found multiple FASTA references");
            }
            foundFastaReference = true;
            referenceName = fastaReference.header;
            reference = fastaReference.sequence;
            referencePosMin = 0;
            referencePosMax = (uint32_t)reference.size() - 1;
        }
    }

    if (foundFastaReference == true) {
        if (verbose) {
            std::cout << ME << "Found FASTA reference" << std::endl;
        }
    } else {
        throwErrorException("Could not find FASTA reference");
    }
}

void QualEncoder::encodeMappedQualityValues(const MappedRecord &mappedRecord)
{
    // Iterators
    uint32_t mrIdx = 0; // iterator for quality values in the mapped record
    uint32_t qiIdx = 0; // iterator for quantizer indices

    // Iterate through CIGAR string and code the quality values
    std::string cigar = mappedRecord.cigar;
    size_t cigarIdx = 0;
    size_t cigarLen = mappedRecord.cigar.length();
    size_t opLen = 0; // length of current CIGAR operation

    for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++) {
        if (isdigit(cigar[cigarIdx])) {
            opLen = opLen*10 + (size_t)cigar[cigarIdx] - (size_t)'0';
            continue;
        }
        switch (cigar[cigarIdx]) {
        case 'M':
        case '=':
        case 'X':
            // Encode opLen quality values with computed quantizer indices
            for (size_t i = 0; i < opLen; i++) {
                int qualityValue = (int)mappedRecord.qualityValues[mrIdx++];
                int quantizerIndex = quantizerIndices[qiIdx++];
                int qualityValueIndex = uniformQuantizers.at(quantizerIndex).valueToIndex(qualityValue);
                qi += std::to_string(quantizerIndex);
                qvi += std::to_string(qualityValueIndex);
                if (quantizedPrintout == true) {
                    int qualityValueQuantized = uniformQuantizers.at(quantizerIndex).valueToReconstructionValue(qualityValue);
                    qualFile << (char)qualityValueQuantized;
                }
            }
            break;
        case 'I':
        case 'S':
            // Encode opLen quality values with max quantizer index
            for (size_t i = 0; i < opLen; i++) {
                int qualityValue = (int)mappedRecord.qualityValues[mrIdx++];
                int qualityValueIndex = uniformQuantizers.at(QUANTIZER_IDX_MAX).valueToIndex(qualityValue);
                qi += std::to_string(QUANTIZER_IDX_MAX);
                qvi += std::to_string(qualityValueIndex);
                if (quantizedPrintout == true) {
                    int qualityValueQuantized = uniformQuantizers.at(QUANTIZER_IDX_MAX).valueToReconstructionValue(qualityValue);
                    qualFile << (char)qualityValueQuantized;
                }
            }
            break;
        case 'D':
        case 'N':
            qiIdx += opLen;
            break; // do nothing as these bases are not present
        case 'H':
        case 'P':
            break; // these have been clipped
        default:
            std::cout << ME << "CIGAR string: " << cigar << std::endl;
            throwErrorException("Bad CIGAR string");
        }
        opLen = 0;
    }

    if (quantizedPrintout == true) {
        qualFile << std::endl;
    }
}

void QualEncoder::encodeUnmappedQualityValues(const std::string &qualityValues)
{
    uqv += qualityValues;

    if (quantizedPrintout == true) {
        qualFile << qualityValues << std::endl;
    }
}

QualDecoder::QualDecoder(File &cqFile, File &qualFile, const CLIOptions &cliOptions)
    // Class scope
    : verbose(false)
    , cqFile(cqFile)
    , qualFile(qualFile)
    , uncompressedSize(0)
    , compressedSize(0)
    , numBlocks(0)
    , numRecords(0)
    , numMappedRecords(0)
    , numUnmappedRecords(0)
    // Block scope
    , uncompressedSizeOfBlock(0)
    , compressedSizeOfBlock(0)
    , numRecordsInBlock(0)
    , numMappedRecordsInBlock(0)
    , numUnmappedRecordsInBlock(0)
{
    // Check if this class shall be verbose
    if (cliOptions.verbose == true) {
        verbose = true;
    }
}

QualDecoder::~QualDecoder(void)
{
    if (verbose == true) {
        std::cout << ME << "QualDecoder statistics:" << std::endl;
        //std::cout << ME << "  Uncompressed size: " << std::setw(9) << uncompressedSize << std::endl;
        std::cout << ME << "  Compressed size:   " << std::setw(9) << compressedSize << std::endl;
        std::cout << ME << "  Blocks:            " << std::setw(9) << numBlocks << std::endl;
        std::cout << ME << "  Records:           " << std::setw(9) << numRecords << std::endl;
        std::cout << ME << "    Mapped:          " << std::setw(9) << numMappedRecords << std::endl;
        std::cout << ME << "    Unmapped:        " << std::setw(9) << numUnmappedRecords << std::endl;
    }
}

size_t QualDecoder::decodeBlock(void)
{
    size_t ret = 0;

    std::cout << ME << "Decoding block " << numBlocks << std::endl;

    // Reset all variables in block scope
    uncompressedSizeOfBlock = 0;
    compressedSizeOfBlock = 0;
    numRecordsInBlock = 0;
    numMappedRecordsInBlock = 0;
    numUnmappedRecordsInBlock = 0;

    // Read and check block number
    uint64_t blockNumber = 0;
    ret += cqFile.readUint64(&blockNumber);
    if (numBlocks != blockNumber) {
        std::cout << ME << "Block number: " << blockNumber << std::endl;
        throwErrorException("Corrupted block number in CQ file");
    }

    // Read record numbers for this block
    ret += cqFile.readUint64(&numRecordsInBlock);
    ret += cqFile.readUint64(&numMappedRecordsInBlock);
    ret += cqFile.readUint64(&numUnmappedRecordsInBlock);

    std::string qiRLEString("");
    uint64_t numQiBlocks = 0;
    ret += cqFile.readUint64(&numQiBlocks);
    for (size_t i = 0; i < numQiBlocks; i++) {
        uint8_t qiRangeFlag = 0;
        ret += cqFile.readUint8(&qiRangeFlag);
        if (qiRangeFlag == 1) {
            uint32_t qiRangeSize = 0;
            ret += cqFile.readUint32(&qiRangeSize);
            unsigned char *qiRange = (unsigned char *)malloc(qiRangeSize);
            ret += cqFile.read(qiRange, qiRangeSize);
            unsigned int qiRLESize = 0;
            unsigned char *qiRLE= range_decompress_o1(qiRange, &qiRLESize);
            free(qiRange);
            std::string qiRLEStringTmp((char *)qiRLE, qiRLESize);
            free(qiRLE);
            qiRLEString += qiRLEStringTmp;
        } else { // qiRangeFlag == 0
            unsigned int qiRLESize = 0;
            ret += cqFile.readUint32(&qiRLESize);
            unsigned char *qiRLE = (unsigned char *)malloc(qiRLESize);
            ret += cqFile.read(qiRLE, qiRLESize);
            std::string qiRLEStringTmp((char *)qiRLE, qiRLESize);
            free(qiRLE);
            qiRLEString += qiRLEStringTmp;
        }

        size_t qiSize = 0;
        unsigned char *qi = rle_decode((unsigned char *)qiRLEString.c_str(), qiRLEString.length(), &qiSize, QUANTIZER_NUM, (unsigned int)'0');
        //std::string qiString((char *)qi, qiSize);
        //std::cout << ME << "Decoded quantizer indices: " << qiString << std::endl;
        free(qi);
    }

    std::string qviRLEString("");
    uint64_t numQviBlocks = 0;
    ret += cqFile.readUint64(&numQviBlocks);
    for (size_t i = 0; i < numQviBlocks; i++) {
        uint8_t qviRangeFlag = 0;
        ret += cqFile.readUint8(&qviRangeFlag);
        if (qviRangeFlag == 1) {
            uint32_t qviRangeSize = 0;
            ret += cqFile.readUint32(&qviRangeSize);
            unsigned char *qviRange = (unsigned char *)malloc(qviRangeSize);
            ret += cqFile.read(qviRange, qviRangeSize);
            unsigned int qviRLESize = 0;
            unsigned char *qviRLE= range_decompress_o1(qviRange, &qviRLESize);
            free(qviRange);
            std::string qviRLEStringTmp((char *)qviRLE, qviRLESize);
            free(qviRLE);
            qviRLEString += qviRLEStringTmp;
        } else { // qviRangeFlag == 0
            unsigned int qviRLESize = 0;
            ret += cqFile.readUint32(&qviRLESize);
            unsigned char *qviRLE = (unsigned char *)malloc(qviRLESize);
            ret += cqFile.read(qviRLE, qviRLESize);
            std::string qviRLEStringTmp((char *)qviRLE, qviRLESize);
            free(qviRLE);
            qviRLEString += qviRLEStringTmp;
        }

        size_t qviSize = 0;
        unsigned char *qvi = rle_decode((unsigned char *)qviRLEString.c_str(), qviRLEString.length(), &qviSize, QUANTIZER_STEP_MAX, (unsigned int)'0');
        //std::string qviString((char *)qvi, qviSize);
        //std::cout << ME << "Decoded quality value indices: " << qviString << std::endl;
        free(qvi);
    }

    //std::string uqvString("");
    uint64_t numUqvBlocks = 0;
    ret += cqFile.readUint64(&numUqvBlocks);
    for (size_t i = 0; i < numUqvBlocks; i++) {
        uint8_t uqvRangeFlag = 0;
        ret += cqFile.readUint8(&uqvRangeFlag);
        if (uqvRangeFlag == 1) {
            uint32_t uqvRangeSize = 0;
            ret += cqFile.readUint32(&uqvRangeSize);
            unsigned char *uqvRange = (unsigned char *)malloc(uqvRangeSize);
            ret += cqFile.read(uqvRange, uqvRangeSize);
            unsigned int uqvSize = 0;
            unsigned char *uqv= range_decompress_o1(uqvRange, &uqvSize);
            free(uqvRange);
            //std::string uqvStringTmp((char *)uqv, uqvSize);
            free(uqv);
            //uqvString += uqvStringTmp;
        } else { // uqvRangeFlag == 0
            unsigned int uqvSize = 0;
            ret += cqFile.readUint32(&uqvSize);
            unsigned char *uqv = (unsigned char *)malloc(uqvSize);
            ret += cqFile.read(uqv, uqvSize);
            //std::string uqvStringTmp((char *)uqv, uqvSize);
            free(uqv);
            //uqvString += uqvStringTmp;
        }
        //std::cout << ME << "Decoded unmapped quality values: " << uqvString << std::endl;
    }

    // Update sizes & counters
    compressedSizeOfBlock += ret;
    compressedSize += compressedSizeOfBlock;
    numBlocks++;
    numRecords += numRecordsInBlock;
    numMappedRecords += numMappedRecordsInBlock;
    numUnmappedRecords += numUnmappedRecordsInBlock;

    // Print statistics for this block
    if (verbose == true) {
        std::cout << ME << "Decoded block " << blockNumber << std::endl;
        //std::cout << ME << "  Uncompressed size: " << uncompressedSizeOfBlock << std::endl;
        std::cout << ME << "  Compressed size: " << std::setw(9) << compressedSizeOfBlock << std::endl;
        std::cout << ME << "  Records:         " << std::setw(9) << numRecordsInBlock << std::endl;
        std::cout << ME << "    Mapped:        " << std::setw(9) << numMappedRecordsInBlock << std::endl;
        std::cout << ME << "    Unmapped:      " << std::setw(9) << numUnmappedRecordsInBlock << std::endl;
    }

    return ret;
}


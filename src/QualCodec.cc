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
#include "Common/Exceptions.h"
#include "Common/debug.h"
#include "Compressors/range/range.h"
#include "Compressors/rle/rle.h"
#include <limits>
#include <stdlib.h>

static const int QV_ALPHABET_SIZE = 50;
static const int QV_OFFSET = 33;
static const int QV_MAX = 83;
static const int QV_MIN = QV_OFFSET;

static const unsigned int NUM_QUANTIZERS = 7; // 2-8 steps
static const unsigned int QUANTIZER_IDX_MIN = 0;
static const unsigned int QUANTIZER_IDX_MAX = NUM_QUANTIZERS-1;

QualEncoder::QualEncoder(File &cqFile,
                         const std::vector<FASTAReference> &fastaReferences,
                         const unsigned int &polyploidy)
    : fastaReferences(fastaReferences)
    , numBlocks(0)
    , numMappedRecords(0)
    , numUnmappedRecords(0)
    , cqFile(cqFile)
    , qi("")
    , qvi("")
    , reference("")
    , referencePosMin(std::numeric_limits<uint32_t>::max())
    , referencePosMax(std::numeric_limits<uint32_t>::min())
    , mappedRecordQueue()
    , observedNucleotides()
    , observedQualityValues()
    , observedPosMin(std::numeric_limits<uint32_t>::max())
    , observedPosMax(std::numeric_limits<uint32_t>::min())
    , genotyper(polyploidy, NUM_QUANTIZERS, QUANTIZER_IDX_MIN, QUANTIZER_IDX_MAX, QV_OFFSET)
    , uniformQuantizers()
    , quantizerIndices()
    , quantizerIndicesPosMin(std::numeric_limits<uint32_t>::max())
    , quantizerIndicesPosMax(std::numeric_limits<uint32_t>::min())
{
    // Init uniform quantizers
    std::cout << ME << "Initializing " << NUM_QUANTIZERS << " uniform quantizers" << std::endl;
    for (unsigned int i = 0; i < NUM_QUANTIZERS; i++) {
        unsigned int numberOfSteps = i + 2;
        UniformQuantizer uniformQuantizer(QV_MIN, QV_MAX, numberOfSteps);
        //uniformQuantizer.print();
        uniformQuantizers.insert(std::pair<int,UniformQuantizer>(i, uniformQuantizer));
    }
}

QualEncoder::~QualEncoder(void)
{
    std::cout << ME << "Encoded " << numMappedRecords << " mapped record(s)"
              << " and " << numUnmappedRecords << " unmapped record(s)"
              << " in " << numBlocks << " block(s)" << std::endl;
}

void QualEncoder::startBlock(void)
{
    qi = "";
    qvi = "";
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

    std::cout << ME << "Starting block " << numBlocks << std::endl;
}

void QualEncoder::addUnmappedRecordToBlock(const SAMRecord &samRecord)
{
    encodeUnmappedQualityValues(std::string(samRecord.qual));
    numUnmappedRecords++;
}

void QualEncoder::addMappedRecordToBlock(const SAMRecord &samRecord)
{
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
    // TODO: this loop consumes 99% of the time of this function
    while (observedPosMin < mappedRecord.posMin) {
        //std::cerr << observedPosMin << ",";
        int k = genotyper.computeQuantizerIndex(observedNucleotides[0], observedQualityValues[0]);
        //std::cerr << k << std::endl;
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

    numMappedRecords++;
}

size_t QualEncoder::finishBlock(void)
{
    size_t ret = 0;

    std::cout << ME << "Finishing block " << numBlocks << std::endl;

    // Compute all remaining quantizers
    while (observedPosMin <= observedPosMax) {
        //std::cerr << observedPosMin << ",";
        int k = genotyper.computeQuantizerIndex(observedNucleotides[0], observedQualityValues[0]);
        //std::cerr << k << std::endl;
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

    // RLE encoding for quantizer indices
    std::cout << ME << "Run-length encoding quantizer indices" << std::endl;
    //std::cout << ME << "Quantizer indices: " << qi << std::endl;
    std::string qiRLE("");
    rle_encode(qi, qiRLE, NUM_QUANTIZERS, (unsigned int)'0');

    // Range encoding for quantizer indices
    std::cout << ME << "Range encoding quantizer indices" << std::endl;
    unsigned int qiRangeSize = 0;
    unsigned char *qiRange = range_compress_o1((unsigned char *)qiRLE.c_str(), qiRLE.size(), &qiRangeSize);
    ret += cqFile.writeUint64(qiRangeSize);
    ret += cqFile.write(qiRange, qiRangeSize);
    free(qiRange);

    // RLE encoding for quality value indices
    std::cout << ME << "Run-length encoding quality value indices" << std::endl;
    //std::cout << ME << "Quality value indices: " << qvi << std::endl;
    std::string qviRLE("");
    rle_encode(qvi, qviRLE, NUM_QUANTIZERS+1, (unsigned int)'0');

    // Range encoding for quality value indices
    std::cout << ME << "Range encoding quality value indices" << std::endl;
    //std::cout << ME << "Quality value indices: " << qvi << std::endl;
    unsigned int qviRangeSize = 0;
    unsigned char *qviRange = range_compress_o1((unsigned char *)qviRLE.c_str(), qviRLE.size(), &qviRangeSize);
    ret += cqFile.writeUint64(qviRangeSize);
    ret += cqFile.write(qviRange, qviRangeSize);
    free(qviRange);

    numBlocks++;

    return ret;
}

void QualEncoder::loadFastaReference(const std::string &rname)
{
    // Find FASTA reference for this RNAME
    std::cout << ME << "Searching FASTA reference for: " << rname << std::endl;
    bool foundFastaReference = false;

    for (auto const &fastaReference : fastaReferences) {
        if (fastaReference.header == rname) {
            if (foundFastaReference == true) {
                throwErrorException("Found multiple FASTA references");
            }
            foundFastaReference = true;
            reference = fastaReference.sequence;
            referencePosMin = 0;
            referencePosMax = (uint32_t)reference.size() - 1;
        }
    }

    if (foundFastaReference == true) {
        std::cout << ME << "Found FASTA reference" << std::endl;
    } else {
        throwErrorException("Could not find FASTA reference");
    }
}

void QualEncoder::encodeMappedQualityValues(const MappedRecord &mappedRecord)
{
    //std::cout << ME << "Encoding: " << mappedRecord.posMin << " " << mappedRecord.cigar << " " << mappedRecord.qualityValues << std::endl;

    // iterators
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
            throwErrorException("Bad CIGAR string");
        }
        opLen = 0;
    }
}

void QualEncoder::encodeUnmappedQualityValues(const std::string &qualityValues)
{
    // TODO
}

QualDecoder::QualDecoder(File &cqFile,
                         File &qualFile,
                         const std::vector<FASTAReference> &fastaReferences)
    : cqFile(cqFile)
    , fastaReferences(fastaReferences)
    , qualFile(qualFile)
{

}

QualDecoder::~QualDecoder(void)
{
    // empty
}

void QualDecoder::decodeBlock(void)
{
    std::cout << ME << "Range decoding quantizer indices" << std::endl;
    uint64_t qiRangeSize = 0;
    cqFile.readUint64(&qiRangeSize);
    unsigned char *qiRange = (unsigned char *)malloc(qiRangeSize);
    cqFile.read(qiRange, qiRangeSize);
    unsigned int qiRLESize = 0;
    unsigned char *qiRLE= range_decompress_o1(qiRange, qiRangeSize, &qiRLESize);

    std::cout << ME << "Run-length decoding quantizer indices" << std::endl;
    std::string qi("");
    std::string qiRLEString((char*)qiRLE, qiRLESize);
    rle_decode(qiRLEString, qi, NUM_QUANTIZERS, (unsigned int)'0');
    std::cout << ME << "Decoded quantizer indices: " << qi << std::endl;

    std::cout << ME << "Range decoding quality value indices" << std::endl;
    uint64_t qviRangeSize = 0;
    cqFile.readUint64(&qviRangeSize);
    unsigned char *qviRange = (unsigned char *)malloc(qviRangeSize);
    cqFile.read(qviRange, qviRangeSize);
    unsigned int qviRLESize = 0;
    unsigned char *qviRLE = range_decompress_o1(qviRange, qviRangeSize, &qviRLESize);

    std::cout << ME << "Run-length decoding quality value indices" << std::endl;
    std::string qvi("");
    std::string qviRLEString((char*)qviRLE, qviRLESize);
    rle_decode(qviRLEString, qvi, NUM_QUANTIZERS+1, (unsigned int)'0');
    std::cout << ME << "Decoded quality value indices: " << qvi << std::endl;

    qualFile.write((unsigned char *)qvi.c_str(), qvi.size());
}


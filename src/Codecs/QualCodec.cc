/** @file QualCodec.cc
 *  @brief This file contains the implementations of the QualEncoder and
 *         QualDecoder classes.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: what (who)
 */

#include "Codecs/QualCodec.h"
#include "common.h"
#include "Exceptions.h"
#include <limits>

static const int Q_ALPHABET_SIZE = 50;
static const int Q_OFFSET = 33;
static const int Q_MAX = 83;
static const int Q_MIN = Q_OFFSET;

static const unsigned int NUM_QUANTIZERS = 7; // 2-8 steps

QualEncoder::QualEncoder(ofbitstream &ofbs,
                         const std::vector<FASTAReference> &fastaReferences,
                         const int &polyploidy)
    : fastaReferences(fastaReferences)
    , numBlocks(0)
    , numMappedRecords(0)
    , numUnmappedRecords(0)
    , ofbs(ofbs)
    , reference("")
    , referencePosMin(std::numeric_limits<uint32_t>::max())
    , referencePosMax(std::numeric_limits<uint32_t>::min())
    , mappedRecordQueue()
    , observedNucleotides()
    , observedQualityValues()
    , observedPosMin(std::numeric_limits<uint32_t>::max())
    , observedPosMax(std::numeric_limits<uint32_t>::min())
    , genotyper(polyploidy, NUM_QUANTIZERS)
    , genotyper2()
    , uniformQuantizers()
    , quantizerIndices()
    , quantizerIndicesPosMin(std::numeric_limits<uint32_t>::max())
    , quantizerIndicesPosMax(std::numeric_limits<uint32_t>::min())
{
    // init uniform quantizers
    std::cout << ME << "Initializing " << NUM_QUANTIZERS << " uniform quantizers" << std::endl;
    for (unsigned int i = 0; i < NUM_QUANTIZERS; i++) {
        unsigned int numberOfSteps = i + 2;
        UniformQuantizer uniformQuantizer(Q_MIN, Q_MAX, numberOfSteps);
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

    // construct a mapped record from the given SAM record
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
        std::cout << ME << "Processing locus: " << observedPosMin << std::endl;
        int k = genotyper.computeQuantizerIndex(observedNucleotides[0], observedQualityValues[0]);
        quantizerIndices.push_back(k);
        quantizerIndicesPosMax = observedPosMin;
        observedNucleotides.erase(observedNucleotides.begin());
        observedQualityValues.erase(observedQualityValues.begin());
        observedPosMin++;
    }

    // check which records from the queue are ready for processing
    while (mappedRecordQueue.front().posMax < observedPosMin) {
        encodeMappedQualityValues(mappedRecordQueue.front());
        uint32_t currFirstPosition = mappedRecordQueue.front().posMin;
        mappedRecordQueue.pop();

        // check if we can shrink the quantizer indices vector
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

    // compute all remaining quantizers
    while (observedPosMin <= observedPosMax) {
        std::cout << ME << "Processing locus: " << observedPosMin << std::endl;
        int k = genotyper.computeQuantizerIndex(observedNucleotides[0], observedQualityValues[0]);
        quantizerIndices.push_back(k);
        quantizerIndicesPosMax = observedPosMin;
        observedNucleotides.erase(observedNucleotides.begin());
        observedQualityValues.erase(observedQualityValues.begin());
        observedPosMin++;
    }

    // process all remaining records from queue
    while (mappedRecordQueue.empty() == false) {
        encodeMappedQualityValues(mappedRecordQueue.front());
        uint32_t currFirstPosition = mappedRecordQueue.front().posMin;
        mappedRecordQueue.pop();

        // check if we can shrink the quantizerIndices vector
        if (mappedRecordQueue.empty() == false) {
            uint32_t nextFirstPosition = mappedRecordQueue.front().posMin;
            if (nextFirstPosition > currFirstPosition) {
                quantizerIndices.erase(quantizerIndices.begin(), quantizerIndices.begin()+nextFirstPosition-currFirstPosition);
                quantizerIndicesPosMin += nextFirstPosition-currFirstPosition;
            }
        }
    }

    numBlocks++;

    return ret;
}

void QualEncoder::loadFastaReference(const std::string &rname)
{
    // find FASTA reference for this RNAME
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
    for (auto const &qualityValue : mappedRecord.qualityValues) {
        int qualityValueQuantized = uniformQuantizers.at(4).valueToIndex(qualityValue);
        //int qualityValueQuantizerIndex = uniformQuantizers.at(4).getIndex(qualityValue);
        //caac.encodeSymbol(qualityValueQuantized);
        //caac.updateModel(qualityValueQuantized);
    }

    // DPCM loop
//     for(std::string::size_type i = 0; i < qual.size(); ++i) {
//         int q = (int)qualityValues[i];
//         int qHat = predictor.predict(memory);
//         int e = q - qHat;
//         int eQuantIdx = maxLloyd.quantize(e);
//         int eQuant = maxLloyd.reconstruct(eQuantIdx);
//         int qQuant = qHat + eQuant;
//         predictor.update(memory, qQuant);
//         maxLloyd.update(qQuant);
//         arithmeticEncoder.encodeSymbol(eQuantIdx);
//     }
}

void QualEncoder::encodeUnmappedQualityValues(const std::string &qualityValues)
{
    // TODO
}

QualDecoder::QualDecoder(ifbitstream &ifbs, std::ofstream &ofs, const std::vector<FASTAReference> &fastaReferences)
    : fastaReferences(fastaReferences)
    , ifbs(ifbs)
    , ofs(ofs)
{

}

QualDecoder::~QualDecoder(void)
{
    // empty
}

void QualDecoder::decodeBlock(void)
{
    // DUMMY
    BYTE byte;
    ifbs.readByte(byte);
    // DUMMY

    // caac.decodeBlock();
}


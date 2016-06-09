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

#include "QualCodec.h"
#include "Exceptions.h"
#include <limits>

static const int Q_ALPHABET_SIZE = 50;
static const int Q_OFFSET = 33;
static const int Q_MAX = 83;
static const int Q_MIN = Q_OFFSET;

QualEncoder::QualEncoder(ofbitstream &ofbs, const std::vector<FASTAReference> &fastaReferences)
    : fastaReferences(fastaReferences)
    , numBlocks(0)
    , numMappedRecords(0)
    , numUnmappedRecords(0)
    , ofbs(ofbs)
    , reference("")
    , mappedRecordQueue()
    , maxPosition(std::numeric_limits<uint32_t>::min())
    , minPosition(std::numeric_limits<uint32_t>::max())
    , observedNucleotides()
    , observedQualityValues()
    , observedSize(0)
    , quantizerIndices()
{
    // empty
}

QualEncoder::~QualEncoder(void)
{
    std::cout << "Encoded " << numMappedRecords << " mapped records"
              << " and " << numUnmappedRecords << " unmapped records"
              << " in " << numBlocks << " block(s)" << std::endl;
}

void QualEncoder::startBlock(void)
{
    reference = "";
    mappedRecordQueue = std::queue<MappedRecord>();
    maxPosition = std::numeric_limits<uint32_t>::min();
    minPosition = std::numeric_limits<uint32_t>::max();
    observedNucleotides.clear();
    observedQualityValues.clear();
    observedSize = 0;
    quantizerIndices.clear();

    std::cout << "Starting block " << numBlocks << std::endl;
}

void QualEncoder::addUnmappedRecordToBlock(const SAMRecord &samRecord)
{
    encodeUnmappedRecord(std::string(samRecord.qual));
    numUnmappedRecords++;
}

void QualEncoder::addMappedRecordToBlock(const SAMRecord &samRecord)
{
    // If the reference is empty, then this is the first mapped record being
    // added to the current block. If so, we load the corresponding reference
    // and set the position offset (a.k.a. 'minPosition').
    if (reference.empty()) {
        loadFastaReference(std::string(samRecord.rname));
        minPosition = samRecord.pos - 1; // SAM format counts from 1
    }

    // construct mapped record
    MappedRecord mappedRecord(samRecord, minPosition);

    // Check if we need to update the max mapping position for records in this
    // block. If so, we also need to update the size of our observation vectors.
    if (mappedRecord.lastPosition > maxPosition) {
        maxPosition = mappedRecord.lastPosition;
        observedSize = maxPosition - minPosition + 1;
        observedNucleotides.resize(observedSize);
        observedQualityValues.resize(observedSize);
    }

    // Parse CIGAR string and assign nucleotides and QVs from the record to
    // their positions within the reference; then push the current record to
    // the queue.
    mappedRecord.extractObservations(observedNucleotides, observedQualityValues);
    mappedRecordQueue.push(mappedRecord);
    //uint32_t pos = minPosition;
    //std::cout << "Observations: " << std::endl;
    //for (auto const &observedColumn : observedNucleotides) {
    //    std::cout << pos << ": " << observedColumn << std::endl;
    //    pos++;
    //}

    // Check which genomic positions are ready for processing. Then compute
    // the quantizer indices for these positions and shrink the observation
    // vectors
    while (minPosition < mappedRecord.firstPosition) {
        //int k = computeQuantizerIndex(reference[minPosition], observedNucleotides[minPosition], observedQualityValues[minPosition]);
        //quantizerIndices.push_back(k);
        observedNucleotides.erase(observedNucleotides.begin());
        observedQualityValues.erase(observedQualityValues.begin());
        minPosition++;
    }

    // check which records from the queue are ready for processing
    while (mappedRecordQueue.front().lastPosition < minPosition) {
        encodeMappedRecord(mappedRecordQueue.front());
        mappedRecordQueue.pop();
    }

    numMappedRecords++;
}

size_t QualEncoder::finishBlock(void)
{
    size_t ret = 0;

    std::cout << "Finishing block " << numBlocks << std::endl;

    // compute all remaining quantizers for position(s) minPosition-maxPosition
    while (minPosition <= maxPosition) {
        //int k = computeQuantizerIndex(reference[minPosition], observedNucleotides[minPosition], observedQualityValues[minPosition])
        //quantizerIndices.push_back(k);
        observedNucleotides.erase(observedNucleotides.begin());
        observedQualityValues.erase(observedQualityValues.begin());
        minPosition++;
    }

    // process all remaining records from queue
    while (mappedRecordQueue.empty() == false) {
        encodeMappedRecord(mappedRecordQueue.front());
        mappedRecordQueue.pop();
    }

    numBlocks++;

    return ret;
}

void QualEncoder::loadFastaReference(const std::string &rname)
{
    // find FASTA reference for this RNAME
    std::cout << "Searching FASTA reference for: " << rname << std::endl;
    bool foundFastaReference = false;

    for (auto const &fastaReference : fastaReferences) {
        if (fastaReference.header == rname) {
            if (foundFastaReference == true) {
                throwErrorException("Found multiple FASTA references");
            }
            foundFastaReference = true;
            reference = fastaReference.sequence;
        }
    }

    if (foundFastaReference == true) {
        std::cout << "Found FASTA reference" << std::endl;
    } else {
        throwErrorException("Could not find FASTA reference");
    }
}

void QualEncoder::encodeMappedRecord(const MappedRecord &mappedRecord)
{
    // Markov-chain prediction for current current qual vector
//     std::vector<int> err;
//     std::vector<int> memory(PREDICTOR_MEMORY_SIZE, -1);
// 
//     for(std::string::size_type i = 0; i < qual.size(); ++i) {
//         int q = (int)qual[i];
// 
//         if (i < PREDICTOR_MEMORY_SIZE) {
//             std::fill(memory.begin(), memory.end(), -1); // this is not necessary
//             err.push_back(q);
//         } else {
//             // fill memory
//             for (size_t m = 0; m < PREDICTOR_MEMORY_SIZE; m++) {
//                 memory[m] = ((int)qual[i-1-m]);
//             }
// 
//             // predict current q value, update predictor, and compute prediction
//             // error
//             int qHat = predictor.predict(memory);
//             predictor.update(memory, q);
//             int e = q - qHat;
//             err.push_back(e);
//         }
//     }

    // TODO: DPCM loop
//     for(std::string::size_type i = 0; i < qual.size(); ++i) {
//         int q = (int)qual[i];
//         int qHat = predictor.predict();
//         int e = q - qHat;
//         eQuant = maxLloyd.quantize(e);
//         int qQuant = qHat + e;
//         maxLloyd.updateModel(e);
//         predictor.update(qQuant);
//     }
}

void QualEncoder::encodeUnmappedRecord(const std::string &qual)
{
    // empty
}

QualDecoder::QualDecoder(ifbitstream &ifbs, std::ofstream &ofs, const std::vector<FASTAReference> &fastaReferences)
    : fastaReferences(fastaReferences)
    , ifbs(ifbs)
    , ofs(ofs)
{

}

QualDecoder::~QualDecoder(void)
{

}


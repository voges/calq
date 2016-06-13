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

static const char alleleAlphabet[] = {'A','C','G','T','N'};
static const unsigned int alleleAlphabetLen = 5;
// TODO: what about polyploidy?

static void allele2genotype(std::map<std::string, double> &genotypeLikelihoods,
                            std::map<char, double> baseLikelihoods,
                            const unsigned int &depth,
                            const unsigned int &offset)
{
    static std::vector<char> genotype;
    double p = 0.0;
    std::string s("");

    if (depth == 0 ) {
        // we are using the log likelihood to avoid numerical problems
        // TODO: replace w/ std::log1p (more accurate)
        p = log((baseLikelihoods[genotype[0]] + baseLikelihoods[genotype[1]]) / 2.0);
        s = std::string() + genotype[0]+genotype[1];
        genotypeLikelihoods[s] += p;
        return;
    }

    for (unsigned int i = offset; i <= (alleleAlphabetLen - depth); i++) {
        genotype.push_back(alleleAlphabet[i]);
        allele2genotype(genotypeLikelihoods, baseLikelihoods, depth-1, i);
        genotype.pop_back();
    }
}

static int computeQuantizerIndex(const char &reference,
                                 const std::string &observedNucleotides,
                                 const std::string &observedQualityValues)
{
    // a map containing the likelihood of each of the possible alleles
    // NOTE: can we work with integers (counts) instead??
    std::map<char, double> baseLikelihoods;

    // a map containing the likelihood of each of the possible genotypes
    // NOTE: can we work with integers (counts) instead??
    std::map<std::string, double> genotypeLikelihoods;

    for (size_t i = 0; i < observedNucleotides.length(); i++) {
        char base = (char)observedNucleotides[i];
        double qualityValue = (double)(observedQualityValues[i]) - Q_OFFSET;

        double p = pow(10.0, -qualityValue/10.0);

        for(unsigned int a = 0; a < alleleAlphabetLen; a++) {
            if (alleleAlphabet[a] == base) {
                baseLikelihoods.insert(std::pair<char, double>(alleleAlphabet[a], 1-p));
            } else {
                baseLikelihoods.insert(std::pair<char, double>(alleleAlphabet[a], p));
            }
        }

        allele2genotype(genotypeLikelihoods, baseLikelihoods, 2, 0);
    }

    // normalize the genotype likelihoods applying the softmax
    double cum = 0.0;
    for (auto &genotypeLikelihood : genotypeLikelihoods) {
        // TODO: replace with std::expm1 (more accurate)
        genotypeLikelihood.second = exp(genotypeLikelihood.second);
        cum += genotypeLikelihood.second;
    }
    for (auto &genotypeLikelihood : genotypeLikelihoods) {
        genotypeLikelihood.second /= cum;
    }

    // compute the entropy
    double entropy = 0.0;
    for (auto &genotypeLikelihood: genotypeLikelihoods) {
        if (genotypeLikelihood.second != 0)
            entropy -= genotypeLikelihood.second * log(genotypeLikelihood.second);
    }

    return 0;
}

QualEncoder::QualEncoder(ofbitstream &ofbs, const std::vector<FASTAReference> &fastaReferences)
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
    , quantizerIndices()
    , quantizerIndicesPosMin(std::numeric_limits<uint32_t>::max())
    , quantizerIndicesPosMax(std::numeric_limits<uint32_t>::min())
{
    // empty
}

QualEncoder::~QualEncoder(void)
{
    std::cout << "Encoded " << numMappedRecords << " mapped record(s)"
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

    std::cout << "Starting block " << numBlocks << std::endl;
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
    // (we name then 'observations') to their positions within the reference;
    // then push the current record to the queue.
    mappedRecord.extractObservations(observedPosMin, observedNucleotides, observedQualityValues);
    mappedRecordQueue.push(mappedRecord);

    // Check which genomic positions are ready for processing. Then compute
    // the quantizer indices for these positions and shrink the observation
    // vectors
    // TODO: this loop consumes 99% of the time of this function
    while (observedPosMin < mappedRecord.posMin) {
        int k = computeQuantizerIndex(reference[observedPosMin], observedNucleotides[0], observedQualityValues[0]);
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

    std::cout << "Finishing block " << numBlocks << std::endl;

    // compute all remaining quantizers
    while (observedPosMin <= observedPosMax) {
        int k = computeQuantizerIndex(reference[observedPosMin], observedNucleotides[0], observedQualityValues[0]);
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
    std::cout << "Searching FASTA reference for: " << rname << std::endl;
    bool foundFastaReference = false;

    for (auto const &fastaReference : fastaReferences) {
        if (fastaReference.header == rname) {
            if (foundFastaReference == true) {
                throwErrorException("Found multiple FASTA references");
            }
            foundFastaReference = true;
            reference = fastaReference.sequence;
            referencePosMin = 0;
            referencePosMax = reference.size() - 1;
        }
    }

    if (foundFastaReference == true) {
        std::cout << "Found FASTA reference" << std::endl;
    } else {
        throwErrorException("Could not find FASTA reference");
    }
}

void QualEncoder::encodeMappedQualityValues(const MappedRecord &mappedRecord)
{
    //std::cout << "Encoding mapped QV vector ranging from " << mappedRecord.posMin << "-" << mappedRecord.posMax << ":" << std::endl;
    //std::cout << mappedRecord.qualityValues << std::endl;
    //std::cout << "Quantizer indices ranging from " << quantizerIndicesPosMin << "-" << quantizerIndicesPosMax << ":" << std::endl;
    //for (auto const &quantizerIndex : quantizerIndices) {
    //     std::cout << quantizerIndex;
    //}
    //std::cout << std::endl;
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

}


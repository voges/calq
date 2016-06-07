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

static const int Q_ALPHABET_SIZE = 50;
static const int Q_OFFSET = 33;
static const int Q_MAX = 83;
static const int Q_MIN = Q_OFFSET;

QualEncoder::QualEncoder(ofbitstream &ofbs, const std::vector<FASTAReference> &fastaReferences)
    : g_fastaReferences(fastaReferences)
    , g_numBlocks(0)
    , g_numMappedRecords(0)
    , g_numUnmappedRecords(0)
    , g_ofbs(ofbs)
    , numMappedRecords(0)
    , numUnmappedRecords(0)
    , reference("")
    , mappedRecordQueue()
    , nextReferencePos(0)
    , observedNucleotides()
    , observedQualityValues()
    , quantizerIndices()
{
    // empty
}

QualEncoder::~QualEncoder(void)
{
    std::cout << "Encoded " << g_numMappedRecords << " mapped records";
    std::cout << " and " << g_numUnmappedRecords << " unmapped records";
    std::cout << " in " << g_numBlocks << " block(s)" << std::endl;
}

void QualEncoder::startBlock(const std::string &rname)
{
    // reset member variables that are used per block
    numMappedRecords = 0;
    numUnmappedRecords = 0;
    reference = "";
    mappedRecordQueue = std::queue<MappedRecord>();
    nextReferencePos = 0;
    observedNucleotides.clear();
    observedQualityValues.clear();
    quantizerIndices.clear();

    // find FASTA reference for this RNAME
    std::cout << "Searching FASTA reference for: " << rname << std::endl;
    bool foundFastaReference = false;

    for (auto const &fastaReference : g_fastaReferences) {
        if (fastaReference.header.find(rname) != std::string::npos) {
            if (foundFastaReference == true) {
                throwErrorException("Found multiple FASTA references");
            }
            foundFastaReference = true;
            reference = fastaReference.sequence;
        }
    }

    if (foundFastaReference == true) {
        std::cout << "Started new block with FASTA reference: " << reference << std::endl;
    } else {
        throwErrorException("Could not find FASTA reference");
    }
}

void QualEncoder::addUnmappedRecordToBlock(const SAMRecord &samRecord)
{
    // deconstruct SAM record
    const std::string rname(samRecord.rname);
    const uint32_t pos = samRecord.pos;
    const std::string cigar(samRecord.cigar);
    const std::string seq(samRecord.seq);
    const std::string qual(samRecord.qual);

    // debug output
    std::cout << "Adding unmapped record to block: " << g_numBlocks << std::endl;
    std::cout << "  rname: " << rname << std::endl;
    std::cout << "  pos:   " << pos << std::endl;
    std::cout << "  cigar: " << cigar << std::endl;
    std::cout << "  seq:   " << seq << std::endl;
    std::cout << "  qual:  " << qual << std::endl;

    numMappedRecords++;
    numUnmappedRecords++;
}

void QualEncoder::addMappedRecordToBlock(const SAMRecord &samRecord)
{
    // deconstruct SAM record
    const std::string rname(samRecord.rname);
    const long pos = samRecord.pos - 1; // SAM format counts from 1
    const std::string cigar(samRecord.cigar);
    const std::string seq(samRecord.seq);
    const std::string qual(samRecord.qual);

    // debug output
    std::cout << "Adding mapped record to block: " << g_numBlocks << std::endl;
    std::cout << "  rname: " << rname << std::endl;
    std::cout << "  pos:   " << pos << std::endl;
    std::cout << "  cigar: " << cigar << std::endl;
    std::cout << "  seq:   " << seq << std::endl;
    std::cout << "  qual:  " << qual << std::endl;

    // check if this alignment is complete
//     if (   (pos == 1)
//         || (cigar.length() == 0 || cigar.compare("*") == 0)
//         || (seq.length() == 0 || seq.compare("*") == 0)
//         || (qual.length() == 0 || qual.compare("*") == 0)) {
//         throwErrorException("Incomplete alignment");
//     }

    // parse CIGAR string and assign nucleotides and QVs to their positions
    // within the reference
    MappedRecord mappedRecord(pos, cigar, seq, qual);
    //assignObservationsToReference(mappedRecord, observedNucleotides, observedQualityValues);
    mappedRecordQueue.push(mappedRecord);

    // check which genomic positions are ready for processing
    while (nextReferencePos < pos) {
        std::cout << "Computing quantizer for position " << nextReferencePos << std::endl;
        //quantizerIndices[nextReferencePos] = computeQuantizerIndex(nextReferencePos);
        nextReferencePos++;
    }

    // check which records in 'recordBuffer' are ready for processing
    while (mappedRecordQueue.front().lastPos < nextReferencePos) {
        encodeMappedRecord(mappedRecordQueue.front());
        mappedRecordQueue.pop();
    }

    numUnmappedRecords++;
    numMappedRecords++;
}

void QualEncoder::encodeMappedRecord(const MappedRecord &mappedRecord)
{
    // TODO: accumulate depths
    //quantizer1.updateDepths(pos, cigar);

    // TODO: quantize the added quality score vector
    // TODO: design/update nested quantizers for genomic columns
    /*for(std::string::size_type i = 0; i < qual.size(); ++i) {
        int q = (int)qual[i];
        int qTilde1 = quantizer1.quantize(q, pos+i);
        qualTilde1.push_back(qTilde1);
    }*/

    // TODO: apply and update quality score mask for current vector
    /*for(size_t i = 0; i < qualTilde1.size(); ++i) {
        int qTilde1 = (int)qual[i];
        int qTilde1Bar = mask.apply(qTilde1);
        qualTilde1Bar.push_back(qTilde1Bar);
    }*/

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
    /*for(std::string::size_type i = 0; i < qual.size(); ++i) {
        int q = (int)qual[i];
        int qHat = predictor.predict();
        int e = q - qHat;
        eQuant = maxLloyd.quantize(e);
        int qQuant = qHat + e;
        maxLloyd.updateModel(e);
        predictor.update(qQuant);
    }*/

    // TODO: adaptive clustering of QS vectors

    // TODO: run-length encoding?

    // TODO: pass QS vector to arithmetic or Huffman coder
    /*for(size_t i = 0; i < tmp.size(); ++i) {
        int q = (int)tmp[i];
        compress();
    }*/

    // ALTERNATIVE: try to predict quality score at a certain position and
    // then quantize the prediction error
}

void QualEncoder::encodeUnmappedRecord(const std::string &qual)
{
    // empty
}

size_t QualEncoder::finishBlock(void)
{
    size_t ret = 0;

    std::cout << "Finished block " << g_numBlocks << std::endl;

    g_numBlocks++;
    g_numMappedRecords += numMappedRecords;
    g_numUnmappedRecords += numUnmappedRecords;

    return ret;
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


/** @file QualEncoder.cc
 *  @brief This file contains the implementation of the QualEncoder class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "QualCodec/QualEncoder.h"
#include "Common/constants.h"
#include "Common/Exceptions.h"
#include "Compressors/range/range.h"
#include "Compressors/rle/rle.h"
#include <math.h>

// static static size_t rangeEncodeAndWrite(cq::CQFile &cqFile, unsigned char *buffer, const size_t &bufferSize)
// {
//     
// }

static size_t rangeEncodeAndWrite(cq::CQFile &cqFile, unsigned char *buffer, const size_t &bufferSize)
{
    size_t ret = 0;

    size_t nrBlocks = (size_t)ceil((double)bufferSize / (double)1*MB);
    ret = cqFile.writeUint64(nrBlocks);

    size_t encodedBytes = 0;
    while (encodedBytes < bufferSize) {
        unsigned int bytesToEncode = 0;
        if ((bufferSize - encodedBytes) > 1*MB) {
            bytesToEncode = 1*MB;
        } else {
            bytesToEncode = bufferSize - encodedBytes;
        }

        unsigned int compressedSize = 0;
        unsigned char *compressed = range_compress_o1(buffer+encodedBytes, (unsigned int)bytesToEncode, &compressedSize);

        if (compressedSize >= bytesToEncode) {
            ret += cqFile.writeUint8(0);
            ret += cqFile.writeUint32(bytesToEncode);
            ret += cqFile.write(buffer+encodedBytes, bytesToEncode);
        } else {
            ret += cqFile.writeUint8(1);
            ret += cqFile.writeUint32(compressedSize);
            ret += cqFile.write(compressed, compressedSize);
        }

        encodedBytes += bytesToEncode;
        free(compressed);
    }

    return ret;
}

cq::QualEncoder::QualEncoder(const unsigned int &polyploidy,
                             const int &qMin,
                             const int &qMax,
                             const std::map<int,Quantizer> &quantizers)
    : startTime_(std::chrono::steady_clock::now())
    , stopTime_(std::chrono::steady_clock::now())

    , compressedMappedQualSize_(0)
    , compressedUnmappedQualSize_(0)
    , nrMappedRecords_(0)
    , nrUnmappedRecords_(0)
    , uncompressedMappedQualSize_(0)
    , uncompressedUnmappedQualSize_(0)

    , posOff_(0)

    , unmappedQual_("")
    , mappedQuantizerIndices_("")
    , mappedQualIndices_("")

    , samPileupDeque_()

    , genotyper_(polyploidy, qMin, qMax, quantizers.size())

    , quantizers_(quantizers)

    , samRecordDeque_()
{
    if (polyploidy == 0) {
       throwErrorException("Polyploidy must be greater than zero");
    }
    if (qMin > qMax) {
        throwErrorException("qMin is greater than qMax");
    }
    if (quantizers.empty() == true) {
        throwErrorException("quantizers is empty");
    }
}

cq::QualEncoder::~QualEncoder(void)
{
    // empty
}

void cq::QualEncoder::startBlock(void)
{
    startTime_ = std::chrono::steady_clock::now();
    stopTime_ = std::chrono::steady_clock::now();

    compressedMappedQualSize_ = 0;
    compressedUnmappedQualSize_ = 0;
    nrMappedRecords_ = 0;
    nrUnmappedRecords_ = 0;
    uncompressedMappedQualSize_ = 0;
    uncompressedUnmappedQualSize_ = 0;

    posOff_ = 0;

    unmappedQual_ = "";
    mappedQuantizerIndices_ = "";
    mappedQualIndices_ = "";

    samPileupDeque_.clear();

    // nothing to be resetted for genotyper_

    // nothing to be resetted for quantizers_

    samRecordDeque_.clear();
}

void cq::QualEncoder::addUnmappedRecordToBlock(const SAMRecord &samRecord)
{
    uncompressedUnmappedQualSize_ += samRecord.qual.length();
    encodeUnmappedQual(samRecord.qual);
    nrUnmappedRecords_++;
}

void cq::QualEncoder::addMappedRecordToBlock(const SAMRecord &samRecord)
{
    uncompressedMappedQualSize_ += samRecord.qual.length();

    if (nrMappedRecords() == 0) {
        posOff_ = samRecord.posMin;
        samPileupDeque_.setPosMin(samRecord.posMin);
        samPileupDeque_.setPosMax(samRecord.posMax);
    }

    if (samRecord.posMax > samPileupDeque_.posMax()) {
        samPileupDeque_.setPosMax(samRecord.posMax);
    }

    samRecord.addToPileupQueue(samPileupDeque_);
    samRecordDeque_.push_back(samRecord);

    while (samPileupDeque_.posMin() < samRecord.posMin) {
        int k = genotyper_.computeQuantizerIndex(samPileupDeque_.front().seq, samPileupDeque_.front().qual);
        mappedQuantizerIndices_ += std::to_string(k);
//         std::cout << samPileupDeque_.front().pos << " ";
//         std::cout << samPileupDeque_.front().seq << " ";
//         std::cout << samPileupDeque_.front().qual << " > " << k << std::endl;
        samPileupDeque_.pop_front();
    }

    while (samRecordDeque_.front().posMax < samPileupDeque_.posMin()) {
        encodeMappedQual(samRecordDeque_.front());
        samRecordDeque_.pop_front();
    }

    nrMappedRecords_++;
}

size_t cq::QualEncoder::finishAndWriteBlock(CQFile &cqFile)
{
    // Compute all remaining quantizers
    while (samPileupDeque_.empty() == false) {
        int k = genotyper_.computeQuantizerIndex(samPileupDeque_.front().seq, samPileupDeque_.front().qual);
        mappedQuantizerIndices_ += std::to_string(k);
//         std::cout << samPileupDeque_.front().pos << " ";
//         std::cout << samPileupDeque_.front().seq << " ";
//         std::cout << samPileupDeque_.front().qual << " > " << k << std::endl;
        samPileupDeque_.pop_front();
    }

    // Process all remaining records from queue
    while (samRecordDeque_.empty() == false) {
        encodeMappedQual(samRecordDeque_.front());
        samRecordDeque_.pop_front();
    }

    CQ_LOG("Finishing and writing unmapped quality values");
    unsigned char *uq = (unsigned char *)unmappedQual_.c_str();
    size_t uqSize = unmappedQual_.length();
    compressedUnmappedQualSize_ += rangeEncodeAndWrite(cqFile, uq, uqSize);

    CQ_LOG("Finishing and writing mapped quantizer indices");
    unsigned char *mqi = (unsigned char *)mappedQuantizerIndices_.c_str();
    size_t mqiSize = mappedQuantizerIndices_.length();
    size_t mqiRLESize = 0;
    unsigned char *mqiRLE = rle_encode(mqi, mqiSize, &mqiRLESize, NR_QUANTIZERS, (unsigned char)'0');
    compressedMappedQualSize_ += rangeEncodeAndWrite(cqFile, mqiRLE, mqiRLESize);
    free(mqiRLE);

    CQ_LOG("Finishing and writing mapped quality value indices");
    unsigned char *mq = (unsigned char *)mappedQualIndices_.c_str();
    size_t mqSize = mappedQualIndices_.length();
    size_t mqRLESize = 0;
    unsigned char *mqRLE = rle_encode(mq, mqSize, &mqRLESize, QUANTIZER_STEPS_MAX, (unsigned char)'0');
    compressedMappedQualSize_ += rangeEncodeAndWrite(cqFile, mqiRLE, mqRLESize);
    free(mqRLE);

    stopTime_ = std::chrono::steady_clock::now();
    auto diffTime = stopTime_ - startTime_;
    auto diffTimeMs = std::chrono::duration_cast<std::chrono::milliseconds>(diffTime).count();
    auto diffTimeS = std::chrono::duration_cast<std::chrono::seconds>(diffTime).count();

    CQ_LOG("Finished block (%zu mapped + %zu unmapped = %zu records)", nrMappedRecords(), nrUnmappedRecords(), nrRecords());
    CQ_LOG("Took %ld ms ~= %ld s", diffTimeMs, diffTimeS);
    CQ_LOG("Speed (uncompressed size/time): %.2f MB/s", ((double)(uncompressedQualSize()/MB))/(double)((double)diffTimeMs/1000));
    CQ_LOG("Bits per quality value: %.4f", ((double)compressedQualSize() * 8)/(double)uncompressedQualSize());
    CQ_LOG("  Mapped: %.4f", ((double)compressedMappedQualSize() * 8)/(double)uncompressedMappedQualSize());
    CQ_LOG("  Unmapped: %.4f", ((double)compressedUnmappedQualSize() * 8)/(double)uncompressedUnmappedQualSize());

    return compressedQualSize();
}

size_t cq::QualEncoder::compressedMappedQualSize(void) const { return compressedMappedQualSize_; }
size_t cq::QualEncoder::compressedUnmappedQualSize(void) const { return compressedUnmappedQualSize_; }
size_t cq::QualEncoder::compressedQualSize(void) const { return compressedMappedQualSize_ + compressedUnmappedQualSize_; }
size_t cq::QualEncoder::nrMappedRecords(void) const { return nrMappedRecords_; }
size_t cq::QualEncoder::nrUnmappedRecords(void) const { return nrUnmappedRecords_; }
size_t cq::QualEncoder::nrRecords(void) const { return nrMappedRecords_ + nrUnmappedRecords_; }
size_t cq::QualEncoder::uncompressedMappedQualSize(void) const { return uncompressedMappedQualSize_; }
size_t cq::QualEncoder::uncompressedUnmappedQualSize(void) const { return uncompressedUnmappedQualSize_; }
size_t cq::QualEncoder::uncompressedQualSize(void) const { return uncompressedMappedQualSize_ + uncompressedUnmappedQualSize_; }

void cq::QualEncoder::encodeMappedQual(const SAMRecord &samRecord)
{
    size_t cigarIdx = 0;
    size_t cigarLen = samRecord.cigar.length();
    size_t opLen = 0; // length of current CIGAR operation
    size_t qualIdx = 0;
    size_t mqiIdx = samRecord.posMin - posOff_;

    for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++) {
       if (isdigit(samRecord.cigar[cigarIdx])) {
           opLen = opLen*10 + (size_t)samRecord.cigar[cigarIdx] - (size_t)'0';
           continue;
       }

       switch (samRecord.cigar[cigarIdx]) {
       case 'M':
       case '=':
       case 'X':
           // Encode opLen quality values with computed quantizer indices
           for (size_t i = 0; i < opLen; i++) {
               int q = (int)samRecord.qual[qualIdx++];
               int mappedQuantizerIndex = mappedQuantizerIndices_[mqiIdx++] - '0';
               int mappedQualIndex = quantizers_.at(mappedQuantizerIndex).valueToIndex(q);
               mappedQualIndices_ += std::to_string(mappedQualIndex);
           }
           break;
       case 'I':
       case 'S':
           // Encode opLen quality values with max quantizer index
           for (size_t i = 0; i < opLen; i++) {
               int q = (int)samRecord.qual[qualIdx++];
               int mappedQualIndex = quantizers_.at(quantizers_.size()-1).valueToIndex(q);
               mappedQualIndices_ += std::to_string(mappedQualIndex);
           }
           break;
       case 'D':
       case 'N':
           mqiIdx += opLen;
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

void cq::QualEncoder::encodeUnmappedQual(const std::string &qual)
{
    unmappedQual_ += qual;
}


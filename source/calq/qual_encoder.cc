#include "calq/qual_encoder.h"

#include <fstream>
#include <iostream>
#include <algorithm>

#include "calq/fasta_file.h"
#include "calq/error_exception_reporter.h"

namespace calq {

QualEncoder::QualEncoder(const Options &options, const std::map<int, Quantizer> &quant)
        : compressedMappedQualSize_(0),
          compressedUnmappedQualSize_(0),
          nrMappedRecords_(0),
          nrUnmappedRecords_(0),
          uncompressedMappedQualSize_(0),
          uncompressedUnmappedQualSize_(0),
          NR_QUANTIZERS(static_cast<int>(options.quantizationMax - options.quantizationMin + 1)),

          qualityValueOffset_(static_cast<int>(options.qualityValueOffset)),
          posOffset_(0),

          unmappedQualityValues_(""),
          mappedQuantizerIndices_(),
          mappedQualityValueIndices_(),

          samPileupDeque_(),

          haplotyper_(options.filterSize, options.polyploidy, options.qualityValueOffset,
                      static_cast<size_t>(NR_QUANTIZERS), 50, 7, 50, options.debug, options.squash, options.filterType),
          genotyper_(static_cast<const int &>(options.polyploidy), static_cast<const int &>(options.qualityValueOffset), NR_QUANTIZERS),
          posCounter(0),

          quantizers_(quant),

          samRecordDeque_(),
          debugOut(options.debug),

          version_(options.version) {
    // Initialize a buffer for mapped quality value indices per quantizer
    for (int i = 0; i < NR_QUANTIZERS; ++i) {
        mappedQualityValueIndices_.emplace_back();
    }
}

QualEncoder::~QualEncoder() = default;

void QualEncoder::addUnmappedRecordToBlock(const SAMRecord &samRecord) {
    encodeUnmappedQual(samRecord.qual);
    uncompressedUnmappedQualSize_ += samRecord.qual.length();
    nrUnmappedRecords_++;
}

void QualEncoder::addMappedRecordToBlock(const SAMRecord &samRecord, const FASTAFile &fasta) {
    if (nrMappedRecords() == 0) {
        posOffset_ = samRecord.posMin;
        samPileupDeque_.setPosMin(samRecord.posMin);
        samPileupDeque_.setPosMax(samRecord.posMax);
    }

    if (samRecord.posMax > samPileupDeque_.posMax()) {
        samPileupDeque_.setPosMax(samRecord.posMax);
    }

    samRecord.addToPileupQueue(&samPileupDeque_, static_cast<size_t>(qualityValueOffset_));
    samRecordDeque_.push_back(samRecord);

    if (version_ == Options::Version::V2) {
        while (samPileupDeque_.posMin() < samRecord.posMin) {
            auto k = static_cast<int>(haplotyper_.push(samPileupDeque_.front().seq, samPileupDeque_.front().qual, samPileupDeque_.front().hq_softcounter,
                                                       fasta.references.at(samRecord.rname)[samPileupDeque_.posMin()]));
            ++posCounter;
            // Start not until pipeline is full
            if (posCounter > haplotyper_.getOffset()) {
                mappedQuantizerIndices_.push_back(k);
            }
            samPileupDeque_.pop_front();
        }
        // Start not until pipeline is full
        while (samRecordDeque_.front().posMax + haplotyper_.getOffset() < samPileupDeque_.posMin()) {
            encodeMappedQual(samRecordDeque_.front());
            samRecordDeque_.pop_front();
        }
    } else {
        while (samPileupDeque_.posMin() < samRecord.posMin) {
            int l = genotyper_.computeQuantizerIndex(samPileupDeque_.front().seq, samPileupDeque_.front().qual);
            mappedQuantizerIndices_.push_back(l);
            samPileupDeque_.pop_front();
        }

        while (samRecordDeque_.front().posMax < samPileupDeque_.posMin()) {
            encodeMappedQual(samRecordDeque_.front());
            samRecordDeque_.pop_front();
        }
    }

    uncompressedMappedQualSize_ += samRecord.qual.length();
    nrMappedRecords_++;
}

void QualEncoder::finishBlock(const FASTAFile &fasta, const std::string &section) {
    // Compute all remaining quantizers
    if (version_ == Options::Version::V2) {
        while (!samPileupDeque_.empty()) {
            auto k = static_cast<int>(haplotyper_.push(samPileupDeque_.front().seq, samPileupDeque_.front().qual, samPileupDeque_.front().hq_softcounter,
                                                       fasta.references.at(section)[samPileupDeque_.posMin()]));
            ++posCounter;
            if (posCounter > haplotyper_.getOffset()) {
                mappedQuantizerIndices_.push_back(k);
            }
            samPileupDeque_.pop_front();
        }

        // Empty pipeline
        size_t offset = std::min(posCounter, haplotyper_.getOffset());
        for (size_t i = 0; i < offset; ++i) {
            int k = static_cast<int>(haplotyper_.push("", "", 0, 'N'));
            mappedQuantizerIndices_.push_back(k);
        }
    } else {
        // Compute all remaining quantizers
        while (!samPileupDeque_.empty()) {
            int k = genotyper_.computeQuantizerIndex(samPileupDeque_.front().seq, samPileupDeque_.front().qual);
            mappedQuantizerIndices_.push_back(k);
            samPileupDeque_.pop_front();
        }
    }

    // TODO(muenteferi): borders of blocks probably too low activity

    // Process all remaining records from queue
    while (!samRecordDeque_.empty()) {
        encodeMappedQual(samRecordDeque_.front());
        samRecordDeque_.pop_front();
    }

    posCounter = 0;
}

size_t QualEncoder::writeBlock(CQFile* cqFile) {
    compressedMappedQualSize_ = 0;
    compressedUnmappedQualSize_ = 0;

    // Write block parameters
    compressedMappedQualSize_ += cqFile->writeUint32(posOffset_);
    compressedMappedQualSize_ += cqFile->writeUint32((uint32_t) qualityValueOffset_);

    // Write inverse quantization LUTs
    // compressedMappedQualSize_ += cqFile->writeQuantizers(quantizers_);

    if (debugOut) {
        std::cerr << "New block. Quantizers:" << std::endl;
        for (auto &q : quantizers_) {
            std::cerr << "Quantizer " << q.first << ", " << q.second.inverseLut().size() << " steps:" << std::endl;
            for (auto &lut : q.second.inverseLut()) {
                std::cerr << lut.first << ": " << lut.second << std::endl;
            }
        }
    }

    // Write unmapped quality values
    auto* uqv = (unsigned char*) unmappedQualityValues_.c_str();
    size_t uqvSize = unmappedQualityValues_.length();
    if (uqvSize > 0) {
        compressedUnmappedQualSize_ += cqFile->writeUint8(0x01);


        std::cerr << "unmapped qvalues:" << std::endl;

        for(auto c : unmappedQualityValues_) {
            std::cerr << uint32_t (c) << std::endl;
        }

        std::cerr <<  std::endl;

        compressedUnmappedQualSize_ += cqFile->writeQualBlock(uqv, uqvSize);
    } else {
        compressedUnmappedQualSize_ += cqFile->writeUint8(0x00);
    }

    // Write mapped quantizer indices
    std::string mqiString;
    for (auto const &mappedQuantizerIndex : mappedQuantizerIndices_) {
        mqiString += std::to_string(mappedQuantizerIndex);
    }
    auto* mqi = (unsigned char*) mqiString.c_str();
    size_t mqiSize = mqiString.length();
    if (mqiSize > 0) {
        compressedMappedQualSize_ += cqFile->writeUint8(0x01);
        std::cerr << "quantizer indices:" << std::endl;

        for(auto c : mappedQuantizerIndices_) {
            std::cerr << uint32_t (c) << std::endl;
        }

        std::cerr <<  std::endl;
        compressedMappedQualSize_ += cqFile->writeQualBlock(mqi, mqiSize);
    } else {
        compressedMappedQualSize_ += cqFile->writeUint8(0x00);
    }

    // Write mapped quality value indices
    for (int i = 0; i < NR_QUANTIZERS; ++i) {
        std::deque<int> mqviStream = mappedQualityValueIndices_[i];
        std::string mqviString;
        for (auto const &mqviInt : mqviStream) {
            mqviString += std::to_string(mqviInt);
        }
        auto* mqvi = (unsigned char*) mqviString.c_str();
        size_t mqviSize = mqviString.length();

        std::cerr << "Step indices" << i << ":" << std::endl;


        for(auto c : mqviStream) {
            std::cerr << uint32_t (c) << std::endl;
        }

        std::cerr <<  std::endl;

        if (mqviSize > 0) {
            compressedMappedQualSize_ += cqFile->writeUint8(0x01);
            compressedMappedQualSize_ += cqFile->writeQualBlock(mqvi, mqviSize);
            // compressedMappedQualSize_ += cqFile->write(mqvi, mqviSize);
        } else {
            compressedMappedQualSize_ += cqFile->writeUint8(0x00);
        }
    }

    return compressedQualSize();
}

size_t QualEncoder::compressedMappedQualSize() const { return compressedMappedQualSize_; }

size_t QualEncoder::compressedUnmappedQualSize() const { return compressedUnmappedQualSize_; }

size_t QualEncoder::compressedQualSize() const { return (compressedMappedQualSize_ + compressedUnmappedQualSize_); }

size_t QualEncoder::nrMappedRecords() const { return nrMappedRecords_; }

size_t QualEncoder::nrUnmappedRecords() const { return nrUnmappedRecords_; }

size_t QualEncoder::nrRecords() const { return (nrMappedRecords_ + nrUnmappedRecords_); }

size_t QualEncoder::uncompressedMappedQualSize() const { return uncompressedMappedQualSize_; }

size_t QualEncoder::uncompressedUnmappedQualSize() const { return uncompressedUnmappedQualSize_; }

size_t QualEncoder::uncompressedQualSize() const { return (uncompressedMappedQualSize_ + uncompressedUnmappedQualSize_); }

void QualEncoder::encodeMappedQual(const SAMRecord &samRecord) {
    size_t cigarIdx = 0;
    size_t cigarLen = samRecord.cigar.length();
    size_t opLen = 0;  // length of current CIGAR operation
    size_t qualIdx = 0;
    size_t quantizerIndicesIdx = samRecord.posMin - posOffset_;

    for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++) {
        if (isdigit(samRecord.cigar[cigarIdx])) {
            opLen = opLen * 10 + (size_t) samRecord.cigar[cigarIdx] - (size_t) '0';
            continue;
        }

        switch (samRecord.cigar[cigarIdx]) {
            case 'M':
            case '=':
            case 'X':
                // Encode opLen quality values with computed quantizer indices
                for (size_t i = 0; i < opLen; i++) {
                    int q = static_cast<int>(samRecord.qual[qualIdx++]) - qualityValueOffset_;
                    int quantizerIndex = mappedQuantizerIndices_[quantizerIndicesIdx++];
                    int qualityValueIndex = quantizers_.at(quantizerIndex).valueToIndex(q);
                    mappedQualityValueIndices_.at(static_cast<size_t>(quantizerIndex)).push_back(qualityValueIndex);
                }
                break;
            case 'I':
            case 'S':
                // Encode opLen quality values with max quantizer index
                for (size_t i = 0; i < opLen; i++) {
                    int q = static_cast<int>(samRecord.qual[qualIdx++]) - qualityValueOffset_;
                    int qualityValueIndex = quantizers_.at(NR_QUANTIZERS - 1).valueToIndex(q);
                    mappedQualityValueIndices_.at(static_cast<size_t>(NR_QUANTIZERS - 1)).push_back(qualityValueIndex);
                }
                break;
            case 'D':
            case 'N':
                quantizerIndicesIdx += opLen;
                break;  // do nothing as these bases are not present
            case 'H':
            case 'P':
                break;  // these have been clipped
            default:
                throwErrorException("Bad CIGAR string");
        }
        opLen = 0;
    }
}

void QualEncoder::encodeUnmappedQual(const std::string &qual) {
    unmappedQualityValues_ += qual;
}

}  // namespace calq

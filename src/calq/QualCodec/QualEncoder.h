/** @file QualEncoder.h
 *  @brief This file contains the definition of the QualEncoder class.
 */

// Copyright 2015-2017 Leibniz Universitaet Hannover

#ifndef CALQ_QUALCODEC_QUALENCODER_H_
#define CALQ_QUALCODEC_QUALENCODER_H_

#include <chrono>
#include <deque>
#include <map>
#include <string>
#include <vector>

#include "IO/CQ/CQFile.h"
#include "IO/SAM/SAMPileupDeque.h"
#include "IO/SAM/SAMRecord.h"
#include "QualCodec/Genotyper.h"
#include "QualCodec/Quantizers/Quantizer.h"
#include "Haplotyper.h"

namespace calq {

class QualEncoder {
 public:
    explicit QualEncoder(const Options &options, const std::map<int, Quantizer> &quant);
    ~QualEncoder();

    void addUnmappedRecordToBlock(const SAMRecord &samRecord);
    void addMappedRecordToBlock(const SAMRecord &samRecord, const FASTAFile &fasta);
    void finishBlock(const FASTAFile &fasta, const std::string &section);
    size_t writeBlock(CQFile* cqFile);

    size_t compressedMappedQualSize() const;
    size_t compressedUnmappedQualSize() const;
    size_t compressedQualSize() const;
    size_t nrMappedRecords() const;
    size_t nrUnmappedRecords() const;
    size_t nrRecords() const;
    size_t uncompressedMappedQualSize() const;
    size_t uncompressedUnmappedQualSize() const;
    size_t uncompressedQualSize() const;

 private:
    void encodeMappedQual(const SAMRecord &samRecord);
    void encodeUnmappedQual(const std::string &qual);

 private:
    // Sizes & counters
    size_t compressedMappedQualSize_;
    size_t compressedUnmappedQualSize_;
    size_t nrMappedRecords_;
    size_t nrUnmappedRecords_;
    size_t uncompressedMappedQualSize_;
    size_t uncompressedUnmappedQualSize_;

    int NR_QUANTIZERS;

    // Quality value offset for this block
    int qualityValueOffset_;

    // 0-based position offset of this block
    uint32_t posOffset_;

    // Buffers
    std::string unmappedQualityValues_;
    std::deque<int> mappedQuantizerIndices_;
    std::vector<std::deque<int> > mappedQualityValueIndices_;

    // Pileup
    SAMPileupDeque samPileupDeque_;

    Haplotyper haplotyper_;

    Genotyper genotyper_;

    size_t posCounter;

    // Quantizers
    std::map<int, Quantizer> quantizers_;

    // Double-ended queue holding the SAM records; records get popped when they
    // are finally encoded
    std::deque<SAMRecord> samRecordDeque_;

    bool debugOut;

    Options::Version version_;
};

}  // namespace calq

#endif  // CALQ_QUALCODEC_QUALENCODER_H_


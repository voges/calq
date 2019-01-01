#ifndef CALQ_QUAL_ENCODER_H_
#define CALQ_QUAL_ENCODER_H_

#include <chrono>
#include <deque>
#include <map>
#include <string>
#include <vector>

#include "cq_file.h"
#include "sam_pileup_deque.h"
#include "sam_record.h"
#include "genotyper.h"
#include "quantizer.h"
#include "haplotyper.h"

namespace calq {

class QualEncoder {
 public:
    explicit QualEncoder(const EncodingOptions &options, const std::map<int, Quantizer> &quant);
    ~QualEncoder();

    void addUnmappedRecordToBlock(const EncodingRead &samRecord);
    void addMappedRecordToBlock(const EncodingRead &samRecord, const std::string& ref);
   void finishBlock(uint32_t pos, const std::string& ref);
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
    void encodeMappedQual(const EncodingRead &samRecord);
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
    std::deque<EncodingRead> samRecordDeque_;

    bool debugOut;

    EncodingOptions::Version version_;
};

}  // namespace calq

#endif  // CALQ_QUAL_ENCODER_H_

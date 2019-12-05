#ifndef CALQ_QUAL_ENCODER_H_
#define CALQ_QUAL_ENCODER_H_

#include <chrono>
#include <deque>
#include <map>
#include <string>
#include <vector>
#include "genotyper.h"
#include "haplotyper.h"
#include "quantizer.h"
#include "sam-pileup.h"

namespace calq {

struct EncodingOptions;
struct DecodingBlock;
enum struct Version;

class QualEncoder {
   public:
    QualEncoder() = delete;
    explicit QualEncoder(const EncodingOptions &options, std::map<int, Quantizer> quantizers, DecodingBlock *out);
    QualEncoder(const QualEncoder &) = delete;
    QualEncoder &operator=(const QualEncoder &) = delete;
    QualEncoder(QualEncoder &&) = delete;
    QualEncoder &operator=(QualEncoder &&) = delete;
    ~QualEncoder() = default;

    void addMappedRecordToBlock(const MinSamRecord &samRecord);
    void finishBlock();
    size_t nrMappedRecords() const;

   private:
    void encodeMappedQual(const MinSamRecord &samRecord);

    size_t nrMappedRecords_;
    int NR_QUANTIZERS;
    uint8_t qualityValueOffset_;  // Quality value offset for this block
    uint32_t posOffset_;          // 0-based position offset of this block
    SamPileupDeque samPileupDeque_;
    Haplotyper haplotyper_;
    Genotyper genotyper_;
    DecodingBlock *out;
    size_t posCounter;
    uint8_t hqSoftClipThreshold;
    std::map<int, Quantizer> quantizers_;

    // Double-ended queue holding the SAM records; records get popped when they
    // are finally encoded
    std::deque<MinSamRecord> samRecordDeque_;

    Version version_;
};

}  // namespace calq

#endif  // CALQ_QUAL_ENCODER_H_

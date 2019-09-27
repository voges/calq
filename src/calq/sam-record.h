/**
 * @file sam-record.h
 */

#ifndef CALQ_SAM_RECORD_H_
#define CALQ_SAM_RECORD_H_

#include <cinttypes>
#include <string>
#include <vector>

namespace calq {

class SamRecord {
   public:
    SamRecord() = delete;
    explicit SamRecord(const std::vector<std::string> &fields);
    SamRecord(const SamRecord &) = default;
    SamRecord &operator=(const SamRecord &) = delete;
    SamRecord(SamRecord &&) = delete;
    SamRecord &operator=(SamRecord &&) = delete;
    ~SamRecord() = default;

    std::string str() const;

    std::string qname;  // Query template NAME
    uint16_t flag;      // bitwise FLAG (uint16_t)
    std::string rname;  // Reference sequence NAME
    uint32_t pos;       // 1-based leftmost mapping POSition (uint32_t)
    uint8_t mapq;       // MAPping Quality (uint8_t)
    std::string cigar;  // CIGAR string
    std::string rnext;  // Ref. name of the mate/NEXT read
    uint32_t pnext;     // Position of the mate/NEXT read (uint32_t)
    int64_t tlen;       // observed Template LENgth (int64_t)
    std::string seq;    // segment SEQuence
    std::string qual;   // QUALity scores
    std::string opt;    // OPTional information
};

}  // namespace calq

#endif  // CALQ_SAM_RECORD_H_

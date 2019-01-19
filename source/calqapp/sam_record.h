#ifndef CALQ_SAM_RECORD_H_
#define CALQ_SAM_RECORD_H_

// -----------------------------------------------------------------------------

#include <cinttypes>

// -----------------------------------------------------------------------------

#include <deque>
#include <string>

// -----------------------------------------------------------------------------

namespace calq {

// -----------------------------------------------------------------------------

class FASTAFile;

// -----------------------------------------------------------------------------

class SAMRecord
{
 public:
    static const int NUM_FIELDS = 12;

    explicit SAMRecord(char *fields[NUM_FIELDS]);
    ~SAMRecord();

    bool isMapped() const;
    void printLong() const;
    void printShort() const;
    void printSeqWithPositionOffset() const;

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

    uint32_t posMin;    // 0-based leftmost mapping position
    uint32_t posMax;    // 0-based rightmost mapping position

 private:
    void check();
    bool mapped_;
};

// -----------------------------------------------------------------------------

}  // namespace calq

// -----------------------------------------------------------------------------

#endif  // CALQ_SAM_RECORD_H_

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

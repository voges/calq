#include "calq/sam_record.h"

#include <queue>

#include "calq/error_exception_reporter.h"
#include "calq/log.h"
#include "calq/fasta_file.h"

namespace calq {

SAMRecord::SAMRecord(char* fields[NUM_FIELDS])
        : qname(fields[0]),
          flag((uint16_t) strtol(fields[1], nullptr, 10)),
          rname(fields[2]),
          pos((uint32_t) strtol(fields[3], nullptr, 10)),
          mapq((uint8_t) strtol(fields[4], nullptr, 10)),
          cigar(fields[5]),
          rnext(fields[6]),
          pnext((uint32_t) strtol(fields[7], nullptr, 10)),
          tlen((int64_t) strtol(fields[8], nullptr, 10)),
          seq(fields[9]),
          qual(fields[10]),
          opt(fields[11]),
          posMin(0),
          posMax(0),
          mapped_(false) {
    check();

    if (mapped_) {
        // Compute 0-based first position and 0-based last position this record
        // is mapped to on the reference used for alignment
        posMin = pos - 1;
        posMax = pos - 1;

        size_t cigarIdx = 0;
        size_t cigarLen = cigar.length();
        uint32_t opLen = 0;  // length of current CIGAR operation

        for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++) {
            if (isdigit(cigar[cigarIdx])) {
                opLen = opLen * 10 + (uint32_t) cigar[cigarIdx] - (uint32_t) '0';
                continue;
            }
            switch (cigar[cigarIdx]) {
                case 'M':
                case '=':
                case 'X':
                    posMax += opLen;
                    break;
                case 'I':
                case 'S':
                    break;
                case 'D':
                case 'N':
                    posMax += opLen;
                    break;
                case 'H':
                case 'P':
                    break;  // these have been clipped
                default:
                    throwErrorException("Bad CIGAR string");
            }
            opLen = 0;
        }
        posMax -= 1;
    }
}

SAMRecord::~SAMRecord() = default;

void SAMRecord::addToPileupQueue(SAMPileupDeque* samPileupDeque_, size_t qualityOffset) const {
    if (samPileupDeque_->empty()) {
        throwErrorException("samPileupQueue is empty");
    }
    if ((samPileupDeque_->posMin() > posMin) || (samPileupDeque_->posMax() < posMax)) {
        throwErrorException("samPileupQueue does not overlap record");
    }

    size_t cigarIdx = 0;
    size_t cigarLen = cigar.length();
    size_t opLen = 0;  // length of current CIGAR operation
    size_t idx = 0;
    size_t pileupIdx = posMin - samPileupDeque_->posMin();

    size_t softclips = 0;

    for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++) {
        if (isdigit(cigar[cigarIdx])) {
            opLen = opLen * 10 + (size_t) cigar[cigarIdx] - (size_t) '0';
            continue;
        }

        const auto HQ_SOFTCLIP_THRESHOLD = static_cast<const char>(29 + qualityOffset);


        switch (cigar[cigarIdx]) {
            case 'M':
            case '=':
            case 'X':
                for (size_t i = 0; i < opLen; i++) {
                    samPileupDeque_->pileups_[pileupIdx].pos
                        = static_cast<uint32_t>(samPileupDeque_->posMin() + pileupIdx);
                    samPileupDeque_->pileups_[pileupIdx].seq += seq[idx];
                    samPileupDeque_->pileups_[pileupIdx].qual += qual[idx];

                    idx++;
                    pileupIdx++;
                }
                break;

            case 'S':
                for (int l = 0; l < static_cast<int>(opLen); ++l) {
                    if (this->qual[idx + l] >= HQ_SOFTCLIP_THRESHOLD) {
                        ++softclips;
                    }
                }
                /* fall through */
            case 'I':
                idx += opLen;
                break;
            case 'D':
            case 'N':
                pileupIdx += opLen;
                break;
            case 'H':
            case 'P':
                break;  // these have been clipped
            default:
                throwErrorException("Bad CIGAR string");
        }

        opLen = 0;
    }

    // Write clips
    for (size_t i = posMin - samPileupDeque_->posMin(); i <= pileupIdx; ++i) {
        samPileupDeque_->pileups_[i].hq_softcounter += softclips;
    }
}

bool SAMRecord::isMapped() const {
    return mapped_;
}

void SAMRecord::printLong() const {
    printShort();
    printf("isMapped: %d, ", mapped_);
    printf("posMin: %d, ", posMin);
    printf("posMax: %d\n", posMax);
}

void SAMRecord::printShort() const {
    printf("%s\t", qname.c_str());
    printf("%d\t", flag);
    printf("%s\t", rname.c_str());
    printf("%d\t", pos);
    printf("%d\t", mapq);
    printf("%s\t", cigar.c_str());
    printf("%s\t", rnext.c_str());
    printf("%d\t", pnext);
    printf("%" PRId64 "\t", tlen);
    printf("%s\t", seq.c_str());
    printf("%s\t", qual.c_str());
    printf("%s\t", opt.c_str());
    printf("\n");
}

void SAMRecord::printSeqWithPositionOffset() const {
    printf("%s %6d-%6d|", rname.c_str(), posMin, posMax);
    for (unsigned int i = 0; i < posMin; i++) {
        printf(" ");
    }
    printf("%s\n", seq.c_str());
}

void SAMRecord::check() {
    // Check all fields
    if (qname.empty()) {
        throwErrorException("qname is empty");
    }
    // flag
    if (rname.empty()) {
        throwErrorException("rname is empty");
    }
    // pos
    // mapq
    if (cigar.empty()) {
        throwErrorException("cigar is empty");
    }
    if (rnext.empty()) {
        throwErrorException("rnext is empty");
    }
    // pnext
    // tlen
    if (seq.empty()) {
        throwErrorException("seq is empty");
    }
    if (qual.empty()) {
        throwErrorException("qual is empty");
    }
    if (opt.empty()) {
        CALQ_LOG("opt is empty");
    }

    // Check if this record is mapped
    if ((flag & 0x4) != 0) { //NOLINT
        mapped_ = false;
    } else {
        mapped_ = true;
        if (rname == "*" || pos == 0 || cigar == "*" || seq == "*" || qual == "*") {
            throwErrorException("Corrupted record");
        }
    }
}

}  // namespace calq

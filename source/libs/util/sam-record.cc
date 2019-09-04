#include "sam-record.h"
#include <queue>
#include <sstream>
#include "calq/errors.h"  // TODO

namespace util {

SamRecord::SamRecord(const std::vector<std::string> &fields)
    : qname(fields[0]),
      flag(static_cast<uint16_t>(std::stoi(fields[1]))),
      rname(fields[2]),
      pos(static_cast<uint32_t>(std::stoi(fields[3]))),
      mapq(static_cast<uint8_t>(std::stoi(fields[4]))),
      cigar(fields[5]),
      rnext(fields[6]),
      pnext(static_cast<uint32_t>(std::stoi(fields[7]))),
      tlen(static_cast<int64_t>(std::stoi(fields[8]))),
      seq(fields[9]),
      qual(fields[10]),
      opt(fields[11]) {}

SamRecord::~SamRecord() = default;

std::string SamRecord::str() const {
    std::stringstream ss;

    ss << qname << " ";
    ss << flag << " ";
    ss << rname << " ";
    ss << pos << " ";
    ss << static_cast<int>(mapq) << " ";  // cast needed because std::cout is overloaded for uint8_t
    ss << cigar << " ";
    ss << rnext << " ";
    ss << pnext << " ";
    ss << tlen << " ";
    ss << seq << " ";
    ss << qual << " ";
    ss << opt;

    return ss.str();
}

}  // namespace util

//    SamRecord::SamRecord(char *fields[NUM_FIELDS])
//    : qname(fields[0]),
//      flag((uint16_t)strtol(fields[1], nullptr, 10)),
//      rname(fields[2]),
//      pos((uint32_t)strtol(fields[3], nullptr, 10)),
//      mapq((uint8_t)strtol(fields[4], nullptr, 10)),
//      cigar(fields[5]),
//      rnext(fields[6]),
//      pnext((uint32_t)strtol(fields[7], nullptr, 10)),
//      tlen((int64_t)strtol(fields[8], nullptr, 10)),
//      seq(fields[9]),
//      qual(fields[10]),
//      opt(fields[11]),
//      posMin(0),
//      posMax(0) {
//    if (isMapped()) {
//        // Compute 0-based first position and 0-based last position this record
//        // is mapped to on the reference used for alignment
//        posMin = pos - 1;
//        posMax = pos - 1;
//
//        size_t cigarIdx = 0;
//        size_t cigarLen = cigar.length();
//        uint32_t opLen = 0;  // Length of current CIGAR operation
//
//        for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++) {
//            if (isdigit(cigar[cigarIdx])) {
//                opLen = opLen * 10 + (uint32_t)cigar[cigarIdx] - (uint32_t)'0';
//                continue;
//            }
//            switch (cigar[cigarIdx]) {
//                case 'M':
//                case '=':
//                case 'X':
//                    posMax += opLen;
//                    break;
//                case 'I':
//                case 'S':
//                    break;
//                case 'D':
//                case 'N':
//                    posMax += opLen;
//                    break;
//                case 'H':
//                case 'P':
//                    break;  // These have been clipped
//                default:
//                    throwErrorException("Bad CIGAR string");
//            }
//            opLen = 0;
//        }
//        posMax -= 1;
//    }
//}
//
// bool SamRecord::isMapped() const {
//    return (flag & 0x4) == 0;  // NOLINT(hicpp-signed-bitwise)
//}

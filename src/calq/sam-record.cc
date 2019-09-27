/**
 * @file sam-record.cc
 */

#include "sam-record.h"
#include <cassert>
#include <queue>
#include <sstream>

namespace calq {

SamRecord::SamRecord(const std::vector<std::string> &fields)
    : qname(""),
      flag(0),
      rname(""),
      pos(0),
      mapq(0),
      cigar(""),
      rnext(""),
      pnext(0),
      tlen(0),
      seq(""),
      qual(""),
      opt("") {
    assert(fields.size() == 12);

    qname = fields[0];
    flag = static_cast<uint16_t>(std::stoi(fields[1]));
    rname = fields[2];
    pos = static_cast<uint32_t>(std::stoi(fields[3]));
    mapq = static_cast<uint8_t>(std::stoi(fields[4]));
    cigar = fields[5];
    rnext = fields[6];
    pnext = static_cast<uint32_t>(std::stoi(fields[7]));
    tlen = static_cast<int64_t>(std::stoi(fields[8]));
    seq = fields[9];
    qual = fields[10];
    opt = fields[11];
}

std::string SamRecord::str() const {
    std::stringstream ss;

    ss << qname << "\t";
    ss << flag << "\t";
    ss << rname << "\t";
    ss << pos << "\t";
    ss << static_cast<int>(mapq) << "\t";  // Cast needed because std::cout is overloaded for uint8_t
    ss << cigar << "\t";
    ss << rnext << "\t";
    ss << pnext << "\t";
    ss << tlen << "\t";
    ss << seq << "\t";
    ss << qual << "\t";
    ss << opt;

    return ss.str();
}

}  // namespace calq

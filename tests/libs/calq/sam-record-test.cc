#include "calq/sam-record.h"
#include <gtest/gtest.h>

TEST(SamRecord, String) {  // NOLINT(cert-err58-cpp)
    std::vector<std::string> fields;

    std::string qname = "QNAME";
    std::string flag = "0";
    std::string rname = "RNAME";
    std::string pos = "0";
    std::string mapq = "0";
    std::string cigar = "CIGAR";
    std::string rnext = "RNEXT";
    std::string pnext = "0";
    std::string tlen = "0";
    std::string seq = "SEQ";
    std::string opt = "OPT";

    fields.push_back(qname);
    fields.push_back(flag);
    fields.push_back(rname);
    fields.push_back(pos);
    fields.push_back(mapq);
    fields.push_back(cigar);
    fields.push_back(rnext);
    fields.push_back(pnext);
    fields.push_back(tlen);
    fields.push_back(seq);
    fields.push_back(opt);

    calq::SamRecord record(fields);
}

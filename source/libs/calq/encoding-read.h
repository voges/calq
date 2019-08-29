/**
 * @file encoding-read.h
 */

#ifndef CALQ_ENCODING_READ_H_
#define CALQ_ENCODING_READ_H_

struct EncodingRead {
    uint32_t posMin;
    uint32_t posMax;
    std::string qvalues;
    std::string cigar;
    std::string sequence;
    std::string reference;
};

#endif  // CALQ_ENCODING_READ_H_

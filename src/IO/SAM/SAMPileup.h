/** @file SAMPileup.h
 *  @brief This file contains the definition of the SAMPileup class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef CQ_SAMPILEUP_H
#define CQ_SAMPILEUP_H

#include <string>

namespace cq {

class SAMPileup {
public:
    explicit SAMPileup(void);
    ~SAMPileup(void);

    bool empty(void) const;
    void clear(void);

    void print(void) const;
    void printQual(void) const;
    void printSeq(void) const;

public:
    uint32_t pos; // 0-based position of this pileup
    std::string qual;
    std::string seq;
};

}

#endif // CQ_SAMPILEUP_H


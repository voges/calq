/** @file Pileup.h
 *  @brief This file contains the definitions of the Pileup and PileupQueue
 *         classes.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef CQ_PILEUP_H
#define CQ_PILEUP_H

#include <deque>
#include <string>

namespace cq {

class Pileup {
public:
    explicit Pileup(void);
    ~Pileup(void);

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

class PileupQueue {
    friend class SAMRecord;

public:
    explicit PileupQueue(void);
    ~PileupQueue(void);

    bool empty(void) const;
    void clear(void);
    const Pileup & front(void) const;
    size_t length(void) const;
    uint32_t posMax(void) const;
    uint32_t posMin(void) const;
    void pop_front(void);
    void setPosMax(const uint32_t &posMax);
    void setPosMin(const uint32_t &posMin);

private:
    std::deque<Pileup> pileups_;
    uint32_t posMax_;
    uint32_t posMin_;
};

}

#endif // CQ_PILEUP_H


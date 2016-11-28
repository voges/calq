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

    void print(void) const;

public:
    std::string seq;
    std::string qual;
};

class PileupQueue {
public:
    explicit PileupQueue(void);
    ~PileupQueue(void);

    void clear(void);
    size_t length(void) const;
    uint32_t posMax(void) const;
    uint32_t posMin(void) const;
    void pop_front(void);
    void setPosMax(const uint32_t &posMax);
    void setPosMin(const uint32_t &posMin);

public:
    std::deque<Pileup> pileups;

private:
    uint32_t posMax_;
    uint32_t posMin_;
};

}

#endif // CQ_PILEUP_H


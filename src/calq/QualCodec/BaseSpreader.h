#ifndef BASESPREADER_H
#define BASESPREADER_H

#include <stddef.h>
#include <vector>

#include "CircularBuffer.h"

class BaseSpreader {
private:
    std::vector<std::pair<size_t, double>> forwardSpread;
    CircularBuffer<double> buffer;
    const size_t MAX_PROPAGATION;
    const size_t MIN_HQ_SOFTCLIPS;
public:
    double push(double score, size_t softclips);
    size_t getOffset() const;

    BaseSpreader(size_t max_prop, size_t min_hq_clips);
};

#endif

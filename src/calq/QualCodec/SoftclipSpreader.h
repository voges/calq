#ifndef BASESPREADER_H
#define BASESPREADER_H

#include <stddef.h>
#include <vector>

#include "CircularBuffer.h"

//Detects areas with high-quality-softclips and spreads their activity score
class SoftclipSpreader {
private:
    std::vector<std::pair<size_t, double>> forwardSpread; //Memory in order to spread to future positions
    CircularBuffer<double> buffer; //Buffer in order to spread to past positions
    const size_t MAX_PROPAGATION; //How far activity scores are supposed to be spread at maximum
    const size_t MIN_HQ_SOFTCLIPS; //How many softclips have to be there to trigger the spreading
public:
    double push(double score, size_t softclips); //Push an activity score and average number of softclips at that position into buffer.
                                                 //Returns the oldest activity score which can't be influenced by new clips anymore.
    size_t getOffset() const;                    //Delay between push of score and processed output

    SoftclipSpreader(size_t max_prop, size_t min_hq_clips);
};

#endif

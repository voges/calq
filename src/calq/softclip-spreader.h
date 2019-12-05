#ifndef CALQ_SOFTCLIP_SPREADER_H_
#define CALQ_SOFTCLIP_SPREADER_H_

#include <cstddef>
#include <utility>
#include <vector>
#include "circular-buffer.h"

namespace calq {

/**
 * Detects areas with high-quality soft clips and spreads their activity score
 */
class SoftclipSpreader {
   public:
    SoftclipSpreader() = delete;
    SoftclipSpreader(size_t maxPropagation, size_t minHqSoftclips, bool isSquashed);
    SoftclipSpreader(const SoftclipSpreader &) = delete;
    SoftclipSpreader &operator=(const SoftclipSpreader &) = delete;
    SoftclipSpreader(SoftclipSpreader &&) = delete;
    SoftclipSpreader &operator=(SoftclipSpreader &&) = delete;
    ~SoftclipSpreader() = default;

    /**
     * Push an activity score and average number of soft clips at that position into buffer. Returns the oldest activity
     * score which can't be influenced by new clips anymore.
     */
    double push(double score, size_t softclips);

    /**
     * Delay between push of score and processed output
     */
    size_t getOffset() const;

   private:
    /**
     * Memory to spread to future positions
     */
    std::vector<std::pair<size_t, double>> forwardSpread;

    /**
     * Buffer to spread to past positions
     */
    CircularBuffer<double> buffer;

    /**
     * Buffer to save raw scores
     */
    CircularBuffer<double> original;

    /**
     * How far activity scores are supposed to be spread at maximum
     */
    const size_t MAX_PROPAGATION;

    /**
     * How many soft clips have to be there to trigger the spreading
     */
    const size_t MIN_HQ_SOFTCLIPS;

    bool squashed;

    /**
     * Squashes activity score between 0 and 1 using the un-squashed activity score and an "anti score" measuring
     * probability of no variant occurring
     */
    static double squash(double activity, double antiActivity) ;
};

}  // namespace calq

#endif  // CALQ_SOFTCLIP_SPREADER_H_

#include "softclip-spreader.h"
#include <algorithm>

namespace calq {

SoftclipSpreader::SoftclipSpreader(const size_t maxPropagation, const size_t minHqSoftclips, const bool isSquashed)
    : buffer(maxPropagation, 0.0),
      original(maxPropagation, 0.0),
      MAX_PROPAGATION(maxPropagation),
      MIN_HQ_SOFTCLIPS(minHqSoftclips),
      squashed(isSquashed) {}

double SoftclipSpreader::push(const double score, const size_t softclips) {
    // Trigger spreading
    if (softclips >= MIN_HQ_SOFTCLIPS) {
        // Radius
        auto clipped = static_cast<int>(std::min(softclips, MAX_PROPAGATION));

        // Remember for future positions
        forwardSpread.emplace_back(clipped + 1, score);

        // Change past positions
        for (auto i = static_cast<int>(buffer.size() - clipped); i < static_cast<int>(buffer.size()); ++i) {
            buffer[i] += score;
        }
    }

    double ownScore = score;

    // Apply all changes to current position in memory from other softclips
    for (auto it = forwardSpread.begin(); it != forwardSpread.end();) {
        ownScore += (*it).second;
        --(*it).first;
        if ((*it).first == 0) {
            it = forwardSpread.erase(it);
        } else {
            ++it;
        }
    }

    double orig = std::min(original.push(score), 1.0);

    return squashed ? squash(buffer.push(ownScore), 1.0 - orig) : buffer.push(ownScore);
}

size_t SoftclipSpreader::getOffset() const { return MAX_PROPAGATION; }

double SoftclipSpreader::squash(const double activity, const double antiActivity) const {
    return activity / (activity + antiActivity);
}

}  // namespace calq

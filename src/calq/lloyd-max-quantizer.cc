#include "lloyd-max-quantizer.h"
#include <algorithm>
#include <utility>

namespace calq {

LloydMaxQuantizer::LloydMaxQuantizer(const size_t steps) {
    steps_ = steps;
    borders_.resize(steps, 0);
    values_.resize(steps, 0);
}

void LloydMaxQuantizer::build(const ProbabilityDistribution& pdf) {
    calcBorders(pdf);
    fillLUT(pdf);
}

void LloydMaxQuantizer::fillLUT(const ProbabilityDistribution& pdf) {
    size_t index = 0;
    size_t pos = pdf.rangeMin();

    while (index < borders_.size()) {
        while (pos <= pdf.rangeMax()) {
            if (pos >= borders_[index]) {
                break;
            }
            this->lut_[pos] = std::pair<int, int>(static_cast<const int&>(index), std::round(values_[index]));
            pos += 1;
        }
        index += 1;
    }

    for (index = 0; index < borders_.size(); ++index) {
        inverseLut_[index] = static_cast<int>(std::round(values_[index]));
    }
}

double LloydMaxQuantizer::centroid(const size_t left, size_t right, const ProbabilityDistribution& pdf) {
    double sum = 0.0;
    double weightSum = 0.0;

    const double THRESHOLD = 0.01;  // Avoid division by zero

    if (right == pdf.rangeMax() + 1) {
        sum += pdf[pdf.size() - 1];
        weightSum += static_cast<double>(right) * pdf[pdf.size() - 1];
        right -= 1;
    }

    for (size_t i = left; i <= right; ++i) {
        sum += pdf[i - pdf.rangeMin()];
        weightSum += static_cast<double>(i) * pdf[i - pdf.rangeMin()];
    }

    return (sum > THRESHOLD) ? (weightSum / sum) : (std::floor(static_cast<double>(left + right) / 2.0));
}

void LloydMaxQuantizer::calcBorders(const ProbabilityDistribution& pdf) {
    // Step 1: init
    double stepSize = pdf.size() / static_cast<double>(steps_);
    for (size_t i = 0; i < steps_; ++i) {
        borders_[i] = pdf.rangeMin() + static_cast<double>(i + 1) * stepSize;
        values_[i] = pdf.rangeMin() + i * stepSize + stepSize / 2.0;
    }

    double change = 0.0;

    // Step 2: Lloyd's II. algorithm
    for (int k = 0; k < static_cast<int>(borders_.size()); ++k) {
        double left = (k == 0) ? pdf.rangeMin() : borders_[k - 1];
        double right = borders_[k];

        // Calc centroid
        double c = centroid(static_cast<size_t>(ceil(left)), static_cast<size_t>(floor(right)), pdf);

        change = std::max(change, std::abs(values_[k] - c));
        values_[k] = c;

        if (k == static_cast<int>(borders_.size() - 1)) {
            constexpr double EPSILON = 0.05;
            if (change < EPSILON) {
                break;
            }
            k = -1;
            change = 0.0;
            continue;
        }

        borders_[k] = (values_[k] + values_[k + 1]) / 2;
    }
}

}  // namespace calq

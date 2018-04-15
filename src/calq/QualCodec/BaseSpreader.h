#include <stddef.h>
#include <vector>

class BaseSpreader {
private:
    constexpr size_t MAX_PROPAGATION=50;
    constexpr size_t MIN_HQ_SOFTCLIPS=7;
    constexpr size_t HQ_SOFTCLIP_THRESHOLD=29;
    std::vector<std::pair<size_t, double>> forwardSpread;
public:
    double push(double score, const std::string& cigar) {

    }


};

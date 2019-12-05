#ifndef CALQ_QUANTIZER_H_
#define CALQ_QUANTIZER_H_

#include <map>
#include <utility>

namespace calq {

class Quantizer {
   public:
    Quantizer();
    explicit Quantizer(std::map<int, int> inverseLut);
    Quantizer(const Quantizer &) = delete;
    Quantizer &operator=(const Quantizer &) = delete;
    Quantizer(Quantizer &&) = default;
    Quantizer &operator=(Quantizer &&) = delete;
    virtual ~Quantizer() = default;

    int valueToIndex(int value) const;
    int indexToReconstructedValue(int index) const;
    const std::map<int, int> &inverseLut() const;

   protected:
    std::map<int, std::pair<int, int>> lut_;  // value->(index,reconstructedValue)
    std::map<int, int> inverseLut_;           // index->reconstructedValue
};

}  // namespace calq

#endif  // CALQ_QUANTIZER_H_

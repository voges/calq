/**
 * @file genotyper.h
 */

#ifndef CALQ_GENOTYPER_H_
#define CALQ_GENOTYPER_H_

#include <map>
#include <string>
#include <vector>

namespace calq {

class Genotyper {
   public:
    Genotyper() = delete;
    Genotyper(int polyploidy, int qualOffset, int numQuantizers);
    Genotyper(const Genotyper&) = delete;
    Genotyper& operator=(const Genotyper&) = delete;
    Genotyper(Genotyper&&) = delete;
    Genotyper& operator=(Genotyper&&) = delete;
    ~Genotyper() = default;

    int computeQuantizerIndex(const std::string& seqPileup, const std::string& qualPileup);
    const std::map<std::string, double>& getGenotypeLikelihoods(const std::string& seqPileup,
                                                                const std::string& qualPileup);

   private:
    void initLikelihoods();
    void resetLikelihoods();
    void computeGenotypeLikelihoods(const std::string& seqPileup, const std::string& qualPileup, size_t depth);

    const std::vector<char> ALLELE_ALPHABET = {'A', 'C', 'G', 'T'};
    const size_t ALLELE_ALPHABET_SIZE = 4;
    const std::vector<char> alleleAlphabet_;
    std::map<char, double> alleleLikelihoods_;
    std::vector<std::string> genotypeAlphabet_;
    std::map<std::string, double> genotypeLikelihoods_;
    const int numQuantizers_;
    const int polyploidy_;
    const int qualOffset_;
    // const bool DEBUG;
};

}  // namespace calq

#endif  // CALQ_GENOTYPER_H_

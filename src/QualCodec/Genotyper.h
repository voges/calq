/** @file Genotyper.h
 *  @brief This file contains the definition of the Genotyper class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#ifndef CALQ_QUALCODEC_GENOTYPER_H_
#define CALQ_QUALCODEC_GENOTYPER_H_

#include <map>
#include <string>
#include <vector>

namespace calq {

class Genotyper {
public:
    Genotyper(const int &polyploidy,
              const int &qualMin,
              const int &qualMax,
              const int &nrQuantizers);
    ~Genotyper(void);

//     void computeAdjustedQual(std::string &adjustedQual,
//                              const std::string &seqPileup,
//                              const std::string &qualPileup);
    double computeEntropy(const std::string &seqPileup,
                          const std::string &qualPileup);
    int computeQuantizerIndex(const std::string &seqPileup,
                              const std::string &qualPileup);
private:
    void initLikelihoods(void);
    void resetLikelihoods(void);
    void computeGenotypeLikelihoods(const std::string &seqPileup,
                                    const std::string &qualPileup,
                                    const size_t &depth);

    const std::vector<char> ALLELE_ALPHABET = {'A', 'C', 'G', 'T'};
    const size_t ALLELE_ALPHABET_SIZE = 4;

    const std::vector<char> alleleAlphabet_;
    std::map<char, double> alleleLikelihoods_;
    std::vector<std::string> genotypeAlphabet_;
    std::map<std::string, double> genotypeLikelihoods_;
    const int nrQuantizers_;
    const int polyploidy_;
    const int qualMin_;
    const int qualMax_;
};

} // namespace calq

#endif // CALQ_QUALCODEC_GENOTYPER_H_


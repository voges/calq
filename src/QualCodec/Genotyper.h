/** @file Genotyper.h
 *  @brief This file contains the definition of the Genotyper class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef CQ_GENOTYPER_H
#define CQ_GENOTYPER_H

#include <map>
#include <string>
#include <vector>

namespace cq {

class Genotyper {
public:
    Genotyper(const unsigned int &polyploidy,
              const unsigned int &quantizerIdxMin,
              const unsigned int &quantizerIdxMax,
              const unsigned int &qualMin,
              const unsigned int &qualMax);
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
private:
    const std::vector<char> ALLELE_ALPHABET = {'A','C','G','T'};
    const size_t ALLELE_ALPHABET_SIZE = 4;

private:
    const std::vector<char> m_alleleAlphabet;
    std::map<char,double> m_alleleLikelihoods;
    std::vector<std::string> m_genotypeAlphabet;
    std::map<std::string,double> m_genotypeLikelihoods;
    const unsigned int m_numQuantizers;
    const unsigned int m_polyploidy;
    const unsigned int m_qualMin;
    const unsigned int m_qualMax;
    const unsigned int m_quantizerIdxMin;
    const unsigned int m_quantizerIdxMax;
};

}

#endif // CQ_GENOTYPER_H


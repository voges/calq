/** @file Genotyper.h
 *  @brief This file contains the definitions of the Genotyper class.
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
              const unsigned int &numQuantizers,
              const unsigned int &quantizerIdxMin,
              const unsigned int &quantizerIdxMax,
              const unsigned int &qualityValueMin,
              const unsigned int &qualityValueMax);
    ~Genotyper(void);

private:
    void initLikelihoods(void);
    void resetLikelihoods(void);

private:
    void computeGenotypeLikelihoods(const std::string &observedNucleotides,
                                    const std::string &observedQualityValues);

public:
    double computeEntropy(const std::string &observedNucleotides,
                          const std::string &observedQualityValues);
    int computeQuantizerIndex(const std::string &observedNucleotides,
                              const std::string &observedQualityValues);
//     void computeAdjustedQualityValues(std::string &adjustedQualityValues,
//                                       const std::string &observedNucleotides,
//                                       const std::string &observedQualityValues);

private:
    const std::vector<char> alleleAlphabet;
    std::map<char,double> alleleLikelihoods;
    std::vector<std::string> genotypeAlphabet;
    std::map<std::string,double> genotypeLikelihoods;
    const unsigned int numQuantizers;
    const unsigned int polyploidy;
    const unsigned int qualityValueMin;
    const unsigned int qualityValueMax;
    const unsigned int quantizerIdxMin;
    const unsigned int quantizerIdxMax;
    bool stats;
};

}

#endif // CQ_GENOTYPER_H


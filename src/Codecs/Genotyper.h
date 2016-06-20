/** @file Genotyper.h
 *  @brief This file contains the definitions of the Genotyper class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: what (who)
 */

#ifndef GENOTYPER_H
#define GENOTYPER_H

#include <map>
#include <string>
#include <vector>

/** @brief Class: Genotyper
 *
 *  The Genotyper class provides methods for computing quantizer indices for
 *  genomic positions given the observed data.
 */
class Genotyper {
public:
    Genotyper(const unsigned int &polyploidy);
    ~Genotyper(void);

    int computeQuantizerIndex(const char &reference,
                              const std::string &observedNucleotides,
                              const std::string &observedQualityValues);

private:
    const std::vector<char> alleleAlphabet;
    std::map<char,double> alleleLikelihoods;
    std::vector<std::string>genotypeAlphabet;
    std::map<std::string,double>genotypeLikelihoods;
    const unsigned int polyploidy;

    void reset(void);
};

#endif // GENOTYPER_H


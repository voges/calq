/** @file Genotyper2.cc
 *  @brief This file contains the implementation of the Genotyper2 class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: what (who)
 */

#include "Codecs/Genotyper2.h"
#include "Exceptions.h"
#include <math.h>

static const int Q_OFFSET = 33;
static const char alleleAlphabet[] = {'A','C','G','T','N'};
static const unsigned int alleleAlphabetLen = 5;

static void allele2genotype(std::map<std::string, double> &genotypeLikelihoods,
                            std::map<char, double> &baseLikelihoods,
                            const unsigned int &depth,
                            const unsigned int &offset)
{
    static std::vector<char> genotype;
    double p = 0.0;
    std::string s("");

    if (depth == 0 ) {
        // we are using the log likelihood to avoid numerical problems
        // TODO: replace w/ std::log1p (more accurate)
        p = log((baseLikelihoods[genotype[0]] + baseLikelihoods[genotype[1]]) / 2.0);
        s = std::string() + genotype[0]+genotype[1];
        genotypeLikelihoods[s] += p;
        return;
    }

    for (unsigned int i = offset; i <= (alleleAlphabetLen - depth); i++) {
        genotype.push_back(alleleAlphabet[i]);
        allele2genotype(genotypeLikelihoods, baseLikelihoods, depth-1, i);
        genotype.pop_back();
    }
}

int Genotyper2::computeQuantizerIndex(const char &reference,
                                      const std::string &observedNucleotides,
                                      const std::string &observedQualityValues)
{
    // sequencing depth at this position
    size_t depth = observedNucleotides.length();

    if (depth != observedQualityValues.length()) {
        throwErrorException("Observation lengths do not match");
    }
    if (depth == 0) {
        return 0; // should return QUANTIZER_INDEX_MAX;
    }

    // a map containing the likelihood of each of the possible alleles
    // NOTE: can we work with integers (counts) instead??
    std::map<char,double> baseLikelihoods;

    // a map containing the likelihood of each of the possible genotypes
    // NOTE: can we work with integers (counts) instead??
    std::map<std::string,double> genotypeLikelihoods;

    for (size_t i = 0; i < depth; i++) {
        char base = (char)observedNucleotides[i];
        double qualityValue = (double)(observedQualityValues[i]) - Q_OFFSET;

        double p = pow(10.0, -qualityValue/10.0);

        for (unsigned int a = 0; a < alleleAlphabetLen; a++) {
            if (alleleAlphabet[a] == base) {
                baseLikelihoods.insert(std::pair<char, double>(alleleAlphabet[a], 1-p));
            } else {
                baseLikelihoods.insert(std::pair<char, double>(alleleAlphabet[a], p));
            }
        }

        allele2genotype(genotypeLikelihoods, baseLikelihoods, 2, 0);
    }

    // normalize the genotype likelihoods applying the softmax
    double cum = 0.0;
    for (auto &genotypeLikelihood : genotypeLikelihoods) {
        // TODO: replace with std::expm1 (more accurate)
        genotypeLikelihood.second = exp(genotypeLikelihood.second);
        cum += genotypeLikelihood.second;
    }
    for (auto &genotypeLikelihood : genotypeLikelihoods) {
        genotypeLikelihood.second /= cum;
    }

    // compute the entropy
    double entropy = 0.0;
    for (auto &genotypeLikelihood: genotypeLikelihoods) {
        if (genotypeLikelihood.second != 0) {
            entropy -= genotypeLikelihood.second * log(genotypeLikelihood.second);
        }
    }

    // debug output
    std::cerr << entropy;

    return 0;
}


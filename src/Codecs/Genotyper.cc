/** @file Genotyper.cc
 *  @brief This file contains the implementation of the Genotyper class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: what (who)
 */

#include "Codecs/Genotyper.h"
#include "Exceptions.h"
#include <math.h>

static const int Q_OFFSET = 33;

static unsigned int combinationsWithoutRepetitions(std::vector<std::string> &genotypeAlphabet,
                                                   const std::vector<char> &alleleAlphabet,
                                                   int *got,
                                                   int n_chosen,
                                                   int len,
                                                   int at,
                                                   int max_types)
{
    int i;
    long count = 0;
    if (n_chosen == len) {
        if (!got) { return 1; }
        std::string tmp("");
        for (i = 0; i < len; i++) {
            tmp += alleleAlphabet[got[i]];
        }
        genotypeAlphabet.push_back(tmp);
        return 1;
    }

    for (i = at; i < max_types; i++) {
        if (got) { got[n_chosen] = i; }
        count += combinationsWithoutRepetitions(genotypeAlphabet, alleleAlphabet, got, n_chosen + 1, len, i, max_types);
    }
    return count;
}

Genotyper::Genotyper(const unsigned int &polyploidy)
    : alleleAlphabet({'A','C','G','T','N'})
    , alleleLikelihoods()
    , genotypeAlphabet()
    , genotypeLikelihoods()
    , polyploidy(polyploidy)
{
    std::cout << "Using polyploidy: " << polyploidy << std::endl;

    for (auto const &allele : alleleAlphabet) {
        alleleLikelihoods.insert(std::pair<char,double>(allele, 0.0));
    }

    size_t alleleAlphabetSize = alleleAlphabet.size();
    int chosen[alleleAlphabetSize];
    combinationsWithoutRepetitions(genotypeAlphabet, alleleAlphabet, chosen, 0, polyploidy, 0, alleleAlphabetSize);

    std::cout << "Using genotype alphabet: ";
    for (auto &genotype : genotypeAlphabet) {
        std::cout << genotype << " ";
        genotypeLikelihoods.insert(std::pair<std::string,double>(genotype, 0.0));
    }
    std::cout << std::endl;
}

Genotyper::~Genotyper(void)
{
    // empty
}

int Genotyper::computeQuantizerIndex(const char &reference,
                                     const std::string &observedNucleotides,
                                     const std::string &observedQualityValues)
{
    const size_t depth = observedNucleotides.length();

    if (depth != observedQualityValues.length()) {
        throwErrorException("Observation lengths do not match");
    }
    if (depth == 0) { return -1; }
    if (depth == 1) { std::cerr << "0"; return 0; } // should return min index

    // a map containing the likelihood of each of the possible alleles
    std::map<char,double> alleleLikelihoods;

    for (size_t d = 0; d < depth; d++) {
        char y = (char)observedNucleotides[d];
        double q = (double)(observedQualityValues[d] - Q_OFFSET);
        double p = pow(10.0, -q/10.0);

        for (auto const &allele : alleleAlphabet) {
            if (allele == y) {
                alleleLikelihoods.insert(std::pair<char,double>(allele, 1-p));
            } else {
                alleleLikelihoods.insert(std::pair<char,double>(allele, p));
            }
        }

        for (auto const &genotype : genotypeAlphabet) {
            // we are using the log likelihood to avoid numerical problems
            // TODO: replace with std::log1p (more accurate)
            double p = 0.0;
            for (size_t i = 0; i < polyploidy; i++) {
                p += alleleLikelihoods[genotype[i]];
            }
            p = log(p / polyploidy);
            genotypeLikelihoods[genotype] += p;
        }
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

    reset();

    return 0;
}

void Genotyper::reset(void)
{
    for (auto &genotypeLikelihood : genotypeLikelihoods) {
        genotypeLikelihood.second = 0.0;
    }

    for (auto &alleleLikelihood : alleleLikelihoods) {
        alleleLikelihood.second = 0.0;
    }
}


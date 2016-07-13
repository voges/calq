/** @file Genotyper.cc
 *  @brief This file contains the implementation of the Genotyper class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: what (who)
 */

#include "Helpers/Genotyper.h"
#include "common.h"
#include "Exceptions.h"
#include <math.h>

// quality value offset
static const int Q_OFFSET = 33; // TODO: is this alway true?

// allele alphabet
static const std::vector<char> ALLELE_ALPHABET = {'A','C','G','T','N'};

static unsigned int combinationsWithoutRepetitions(std::vector<std::string> &genotypeAlphabet,
                                                   const std::vector<char> &alleleAlphabet,
                                                   int *got,
                                                   int nChosen,
                                                   int len,
                                                   int at,
                                                   int maxTypes)
{
    if (nChosen == len) {
        if (!got) { return 1; }
        std::string tmp("");
        for (int i = 0; i < len; i++) {
            tmp += alleleAlphabet[got[i]];
        }
        genotypeAlphabet.push_back(tmp);
        return 1;
    }

    long count = 0;
    for (int i = at; i < maxTypes; i++) {
        if (got) { got[nChosen] = i; }
        count += combinationsWithoutRepetitions(genotypeAlphabet, alleleAlphabet, got, nChosen+1, len, i, maxTypes);
    }

    return count;
}

Genotyper::Genotyper(const unsigned int &polyploidy, const unsigned int &numQuantizers)
    : alleleAlphabet(ALLELE_ALPHABET)
    , alleleLikelihoods()
    , genotypeAlphabet()
    , genotypeLikelihoods()
    , numQuantizers(numQuantizers)
    , polyploidy(polyploidy)
{
    initLikelihoods();
}

Genotyper::~Genotyper(void)
{
    // empty
}

void Genotyper::initLikelihoods(void)
{
    // init map containing the allele likelihoods
    std::cout << ME << "Allele alphabet: ";
    for (auto const &allele : alleleAlphabet) {
        std::cout << allele << " ";
        alleleLikelihoods.insert(std::pair<char,double>(allele, 0.0));
    }
    std::cout << std::endl;

    // init map containing the genotype likelihoods
    int chosen[ALLELE_ALPHABET.size()];
    combinationsWithoutRepetitions(genotypeAlphabet, alleleAlphabet, chosen, 0, polyploidy, 0, ALLELE_ALPHABET.size());

    std::cout << ME << "Genotype alphabet ";
    std::cout << "(" << genotypeAlphabet.size() << " possible genotypes): ";
    for (auto &genotype : genotypeAlphabet) {
        std::cout << genotype << " ";
        genotypeLikelihoods.insert(std::pair<std::string,double>(genotype, 0.0));
    }
    std::cout << std::endl;
}

void Genotyper::resetLikelihoods(void)
{
    for (auto &genotypeLikelihood : genotypeLikelihoods) {
        genotypeLikelihood.second = 0.0;
    }

    for (auto &alleleLikelihood : alleleLikelihoods) {
        alleleLikelihood.second = 0.0;
    }
}

bool Genotyper::computeGenotypeLikelihoods(const std::string &observedNucleotides,
                                           const std::string &observedQualityValues)
{
    resetLikelihoods();

    const size_t depth = observedNucleotides.length();
    if (depth != observedQualityValues.length()) {
        throwErrorException("Observation lengths do not match");
    }
    if (depth == 0) { return false; } // no computation possible for depth=0
    if (depth == 1) { return false; } // finest quantization required

    std::cout << ME << "Depth: " << depth << std::endl;
    std::cout << ME << "Observed nucleotides: " << observedNucleotides << std::endl;
    std::cout << ME << "Observed quality values: " << observedQualityValues << std::endl;

    for (size_t d = 0; d < depth; d++) {
        char y = (char)observedNucleotides[d];
        double q = (double)(observedQualityValues[d] - Q_OFFSET);
        double pStrike = 1 - pow(10.0, -q/10.0);
        double pError = (1-pStrike) / (ALLELE_ALPHABET_SIZE-1);

        for (auto const &allele : alleleAlphabet) {
            if (allele == y) {
                alleleLikelihoods[allele] = pStrike;
            } else {
                alleleLikelihoods[allele] = pError;
            }
        }

        for (auto const &genotype : genotypeAlphabet) {
            double p = 0.0;
            for (size_t i = 0; i < polyploidy; i++) {
                p += alleleLikelihoods[genotype[i]];
            }
            p /= polyploidy;

            // we are using the log likelihood to avoid numerical problems
            genotypeLikelihoods[genotype] += log(p);
        }
    }

    // normalize the genotype likelihoods
    double cum = 0.0;
    for (auto &genotypeLikelihood : genotypeLikelihoods) {
        genotypeLikelihood.second = exp(genotypeLikelihood.second);
        cum += genotypeLikelihood.second;
    }
    for (auto &genotypeLikelihood : genotypeLikelihoods) {
        genotypeLikelihood.second /= cum;
    }

    return true;
}

double Genotyper::computeGenotypeEntropy(const std::string &observedNucleotides,
                                         const std::string &observedQualityValues)
{
    if (!computeGenotypeLikelihoods(observedNucleotides, observedQualityValues)) {
        return -1;
    }

    double entropy = 0.0;
    for (auto &genotypeLikelihood: genotypeLikelihoods) {
        if (genotypeLikelihood.second != 0) {
            entropy -= genotypeLikelihood.second * log(genotypeLikelihood.second);
        }
    }

    std::cout << ME << "Entropy: " << entropy << std::endl;
    return entropy;
}

double Genotyper::computeGenotypeConfidence(const std::string &observedNucleotides,
                                           const std::string &observedQualityValues)
{
    if (!computeGenotypeLikelihoods(observedNucleotides, observedQualityValues)) {
        return -1;
    }

    double largestGenotypeLikelihood = 0.0;
    double secondLargestGenotypeLikelihood = 0.0;
    for (auto &genotypeLikelihood : genotypeLikelihoods) {
        if (genotypeLikelihood.second > secondLargestGenotypeLikelihood) {
            secondLargestGenotypeLikelihood = genotypeLikelihood.second;
        }
        if (secondLargestGenotypeLikelihood > largestGenotypeLikelihood) {
            secondLargestGenotypeLikelihood = largestGenotypeLikelihood;
            largestGenotypeLikelihood = genotypeLikelihood.second;
        }
    }

    double confidence = largestGenotypeLikelihood - secondLargestGenotypeLikelihood;

    std::cout << ME << "Confidence: " << confidence << std::endl;
    return confidence;
}

int Genotyper::computeQuantizerIndex(const std::string &observedNucleotides,
                                     const std::string &observedQualityValues)
{
    double confidence = computeGenotypeConfidence(observedNucleotides, observedQualityValues);
    
    if (confidence == -1 ) {
        std::cout << ME << "Select finest quantizer" << std::endl;
        return -1;
    }

    std::cout << ME << "Confidence: " << confidence << std::endl;
    return (int)confidence;
}


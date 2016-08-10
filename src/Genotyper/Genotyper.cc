/** @file Genotyper.cc
 *  @brief This file contains the implementation of the Genotyper class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "Genotyper.h"
#include "Common/debug.h"
#include "Common/Exceptions.h"
#include "Common/fileSystemHelpers.h"
#include <math.h>

// Allele alphabet
static const std::vector<char> GENOTYPER_ALLELE_ALPHABET = {'A','C','G','T'};

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

Genotyper::Genotyper(const unsigned int &polyploidy, 
                     const unsigned int &numQuantizers,
                     const unsigned int &quantizerIdxMin,
                     const unsigned int &quantizerIdxMax,
                     const unsigned int &qualityValueMin,
                     const unsigned int &qualityValueMax)
    : alleleAlphabet(GENOTYPER_ALLELE_ALPHABET)
    , alleleLikelihoods()
    , genotypeAlphabet()
    , genotypeLikelihoods()
    , numQuantizers(numQuantizers)
    , polyploidy(polyploidy)
    , qualityValueMin(qualityValueMin)
    , qualityValueMax(qualityValueMax)
    , quantizerIdxMin(quantizerIdxMin)
    , quantizerIdxMax(quantizerIdxMax)
{
    initLikelihoods();
}

Genotyper::~Genotyper(void)
{
    // empty
}

void Genotyper::initLikelihoods(void)
{
    // Init map containing the allele likelihoods
    std::cout << ME << "Allele alphabet: ";
    for (auto const &allele : alleleAlphabet) {
        std::cout << allele << " ";
        alleleLikelihoods.insert(std::pair<char,double>(allele, 0.0));
    }
    std::cout << std::endl;

    // Init map containing the genotype likelihoods
    int chosen[GENOTYPER_ALLELE_ALPHABET.size()];
    combinationsWithoutRepetitions(genotypeAlphabet, alleleAlphabet, chosen, 0, polyploidy, 0, GENOTYPER_ALLELE_ALPHABET.size());

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

void Genotyper::computeGenotypeLikelihoods(const std::string &observedNucleotides,
                                           const std::string &observedQualityValues)
{
    int i = 0;
    resetLikelihoods();

    const size_t depth = observedNucleotides.length();
    if (depth != observedQualityValues.length()) {
        throwErrorException("Observation lengths do not match");
    }
    if (depth == 0  || depth == 1) {
        throwErrorException("Depth must be greater than one");
    }

    tempGenotypeLL = (double*)calloc(sizeof(double), 1024);
    for (size_t d = 0; d < depth; d++) {
        char y = (char)observedNucleotides[d];
        double q = (double)(observedQualityValues[d] - qualityValueMin);
        if ((q > qualityValueMax-qualityValueMin) || q < 0) {
            throwErrorException("Quality value out of range");
        }
        double pStrike = 1 - pow(10.0, -q/10.0);

        double pError = (1-pStrike) / (GENOTYPER_ALLELE_ALPHABET.size()-1);
        int i = 0;
        for (auto const &genotype : genotypeAlphabet) {
            double p = 0.0;
            for (size_t i = 0; i < polyploidy; i++) {
                p += (y==genotype[i])? pStrike:pError;
            }
            p /= polyploidy;

            // We are using the log likelihood to avoid numerical problems
            tempGenotypeLL[i++] += log(p);
        }
    }
    for (auto const &genotype : genotypeAlphabet)
        genotypeLikelihoods[genotype] = tempGenotypeLL[i++];
    free(tempGenotypeLL);

    // Nnormalize the genotype likelihoods
    double cum = 0.0;
    for (auto &genotypeLikelihood : genotypeLikelihoods) {
        genotypeLikelihood.second = exp(genotypeLikelihood.second);
        cum += genotypeLikelihood.second;
    }
    for (auto &genotypeLikelihood : genotypeLikelihoods) {
        genotypeLikelihood.second /= cum;
    }
}

double Genotyper::computeGenotypeEntropy(const std::string &observedNucleotides,
                                         const std::string &observedQualityValues)
{
    const size_t depth = observedNucleotides.length();
    if (depth == 0 || depth == 1) {
        return -1.0;
    }

    computeGenotypeLikelihoods(observedNucleotides, observedQualityValues);

    double entropy = 0.0;
    for (auto &genotypeLikelihood: genotypeLikelihoods) {
        if (genotypeLikelihood.second != 0) {
            entropy -= genotypeLikelihood.second * log(genotypeLikelihood.second);
        }
    }

    std::cerr << entropy << ",";
    return entropy;
}

int Genotyper::computeQuantizerIndex(const std::string &observedNucleotides,
                                     const std::string &observedQualityValues)
{
    const size_t depth = observedNucleotides.length();
    //std::cerr << depth << ",";
    //std::cerr << observedNucleotides << ",";
    //std::cerr << observedQualityValues << ",";
    if (depth == 0) {
        return -1; // computation of quantizer index not possible
    }
    if (depth == 1) {
        return quantizerIdxMax;
    }

    computeGenotypeLikelihoods(observedNucleotides, observedQualityValues);

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
    //std::cerr << confidence << ",";

    //if (confidence == 1) return quantizerIdxMin;
    return (int)((1-confidence)*(numQuantizers-1));
}


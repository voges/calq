/** @file Genotyper.cc
 *  @brief This file contains the implementation of the Genotyper class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "QualCodec/Genotyper.h"
#include "Common/CLIOptions.h"
#include "Common/helpers.h"
#include "Common/Exceptions.h"
#include "Common/helpers.h"
#include <math.h>

// Allele alphabet
static const std::vector<char> GENOTYPER_ALLELE_ALPHABET = {'A','C','G','T'};
static const size_t GENOTYPER_ALLELE_ALPHABET_SIZE = 4;

static unsigned int combinationsWithRepetitions(std::vector<std::string> &genotypeAlphabet,
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
        count += combinationsWithRepetitions(genotypeAlphabet, alleleAlphabet, got, nChosen+1, len, i, maxTypes);
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
    , stats(false)
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
    for (auto const &allele : alleleAlphabet) {
        alleleLikelihoods.insert(std::pair<char,double>(allele, 0.0));
    }

    // Init map containing the genotype likelihoods
    int chosen[GENOTYPER_ALLELE_ALPHABET_SIZE];
    combinationsWithRepetitions(genotypeAlphabet, alleleAlphabet, chosen, 0, polyploidy, 0, GENOTYPER_ALLELE_ALPHABET_SIZE);

    CQ_LOG("Initializing genotype alphabet with %zu possible genotypes", genotypeAlphabet.size());
    for (auto &genotype : genotypeAlphabet) {
        genotypeLikelihoods.insert(std::pair<std::string,double>(genotype, 0.0));
    }
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
    resetLikelihoods();

    const size_t depth = observedNucleotides.length();
    if (depth != observedQualityValues.length()) {
        throwErrorException("Observation lengths do not match");
    }
    if (depth == 0  || depth == 1) {
        throwErrorException("Depth must be greater than one");
    }

    double *tempGenotypeLikelihoods = (double *)calloc(genotypeAlphabet.size(), sizeof(double));
    unsigned int itr = 0;
    for (size_t d = 0; d < depth; d++) {
        char y = (char)observedNucleotides[d];
        double q = (double)(observedQualityValues[d] - qualityValueMin);
        if ((q > qualityValueMax-qualityValueMin) || q < 0) {
            //LOG("min,curr,max=%d,%d,%d", 0, (int)q, qualityValueMax-qualityValueMin);
            throwErrorException("Quality value out of range");
        }
        double pStrike = 1 - pow(10.0, -q/10.0);

        double pError = (1-pStrike) / (GENOTYPER_ALLELE_ALPHABET.size()-1);
        itr = 0;
        for (auto const &genotype : genotypeAlphabet) {
            double p = 0.0;
            for (unsigned int i = 0; i < polyploidy; i++) {
                p += (y == genotype[i]) ? pStrike : pError;
            }
            p /= polyploidy;

            // We are using the log likelihood to avoid numerical problems
            tempGenotypeLikelihoods[itr++] += log(p);
        }
    }
    itr = 0;
    for (auto const &genotype : genotypeAlphabet) {
        genotypeLikelihoods[genotype] = tempGenotypeLikelihoods[itr++];
    }
    free(tempGenotypeLikelihoods);

    // Normalize the genotype likelihoods
    double cum = 0.0;
    for (auto &genotypeLikelihood : genotypeLikelihoods) {
        genotypeLikelihood.second = exp(genotypeLikelihood.second);
        cum += genotypeLikelihood.second;
    }
    for (auto &genotypeLikelihood : genotypeLikelihoods) {
        genotypeLikelihood.second /= cum;
    }
}

double Genotyper::computeEntropy(const std::string &observedNucleotides,
                                 const std::string &observedQualityValues)
{
    const size_t depth = observedNucleotides.length();

    if (depth == 0) {
        return -1.0; // computation of entropy not possible
    }
    if (depth == 1) {
        return 0.0; // no information content for one symbol
    }

    computeGenotypeLikelihoods(observedNucleotides, observedQualityValues);

    double entropy = 0.0;
    for (auto &genotypeLikelihood: genotypeLikelihoods) {
        if (genotypeLikelihood.second != 0) {
            entropy -= genotypeLikelihood.second * log(genotypeLikelihood.second);
        }
    }

    return entropy;
}

int Genotyper::computeQuantizerIndex(const std::string &observedNucleotides,
                                     const std::string &observedQualityValues)
{
    const size_t depth = observedNucleotides.length();

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

    return (int)((1-confidence)*(numQuantizers-1));
}

// void Genotyper::computeAdjustedQualityValues(std::string &adjustedQualityValues,
//                                              const std::string &observedNucleotides,
//                                              const std::string &observedQualityValues)
// {
//     //haploid: base N has Q and P -- select max as new QV
//     //diploid: base N has Q and P1,2,3,4 -- select max as new QV
// 
//     for (auto &genotypeLikelihood : genotypeLikelihoods) {
//         if (genotypeLikelihood.second > largestGenotypeLikelihood) {
//             largestGenotypeLikelihood = genotypeLikelihood.second;
//         }
//     }
//     // Convert the probability to phred score
//     double gtLLPhred = -log10(1-largestGenotypeLikelihood);
// 
//     adjustedQualityValues = "";
//     const size_t depth = observedNucleotides.length();
// 
//     if ((depth == 0) || (depth == 1)) {
//         adjustedQualityValues = observedQualityValues;
//     }
// 
//     computeGenotypeLikelihoods(observedNucleotides, observedQualityValues);
// 
//     // Get likelihoods for the alleles
//     std::map<char,std::vector<double>> nucleotideQualityValues;
//     for (auto const &allele : alleleAlphabet) {
//         for (auto const &genotype : genotypeAlphabet) {
//             if (genotype.find(allele) != std::string::npos) {
//                 nucleotideQualityValues[allele].append(genotypeLikelihoods[genotype];
//             }
//         }
//     }
// 
//     // Get the largest likelihood
//     for (auto const &nucleotideQualityValue : nucleotideQualityValues) {
//         //
//     }
// }


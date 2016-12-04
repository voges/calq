/** @file Genotyper.cc
 *  @brief This file contains the implementation of the Genotyper class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#include "QualCodec/Genotyper.h"

#include <math.h>

#include <utility>

#include "Common/Exceptions.h"
#include "Common/log.h"

namespace calq {

static int combinationsWithRepetitions(std::vector<std::string> *genotypeAlphabet,
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
        genotypeAlphabet->push_back(tmp);
        return 1;
    }

    int count = 0;
    for (int i = at; i < maxTypes; i++) {
        if (got) { got[nChosen] = i; }
        count += combinationsWithRepetitions(genotypeAlphabet, alleleAlphabet, got, nChosen+1, len, i, maxTypes);
    }

    return count;
}

Genotyper::Genotyper(const int &polyploidy,
                     const int &qualMin,
                     const int &qualMax,
                     const int &nrQuantizers)
    : alleleAlphabet_(ALLELE_ALPHABET),
      alleleLikelihoods_(),
      genotypeAlphabet_(),
      genotypeLikelihoods_(),
      nrQuantizers_(nrQuantizers),
      polyploidy_(polyploidy),
      qualMin_(qualMin),
      qualMax_(qualMax)
{
    initLikelihoods();
}

Genotyper::~Genotyper(void) {}

double Genotyper::computeEntropy(const std::string &seqPileup,
                                 const std::string &qualPileup)
{
    const size_t depth = seqPileup.length();

    if (depth != qualPileup.length()) {
        throwErrorException("Lengths of seqPileup and qualPileup differ");
    }

    if (depth == 0) {
        return -1.0; // computation of entropy not possible
    }
    if (depth == 1) {
        return 0.0; // no information content for one symbol
    }

    computeGenotypeLikelihoods(seqPileup, qualPileup, depth);

    double entropy = 0.0;
    for (auto &genotypeLikelihood : genotypeLikelihoods_) {
        if (genotypeLikelihood.second != 0) {
            entropy -= genotypeLikelihood.second * log(genotypeLikelihood.second);
        }
    }

    return entropy;
}

int Genotyper::computeQuantizerIndex(const std::string &seqPileup,
                                     const std::string &qualPileup)
{
    const size_t depth = seqPileup.length();

    if (depth != qualPileup.length()) {
        throwErrorException("Lengths of seqPileup and qualPileup differ");
    }

    if (depth == 0) {
        return -1; // computation of quantizer index not possible
    }
    if (depth == 1) {
        return (nrQuantizers_ - 1); // no inference can be made, stay safe
    }

    computeGenotypeLikelihoods(seqPileup, qualPileup, depth);

    double largestGenotypeLikelihood = 0.0;
    double secondLargestGenotypeLikelihood = 0.0;
    for (auto &genotypeLikelihood : genotypeLikelihoods_) {
        if (genotypeLikelihood.second > secondLargestGenotypeLikelihood) {
            secondLargestGenotypeLikelihood = genotypeLikelihood.second;
        }
        if (secondLargestGenotypeLikelihood > largestGenotypeLikelihood) {
            secondLargestGenotypeLikelihood = largestGenotypeLikelihood;
            largestGenotypeLikelihood = genotypeLikelihood.second;
        }
    }

    double confidence = largestGenotypeLikelihood - secondLargestGenotypeLikelihood;
//     std::cout << "confidence: " << confidence << std::endl;
//     std::cout << "k: " << (int)((1-confidence)*(nrQuantizers_-1)) << std::endl;
    return (int)((1-confidence)*(nrQuantizers_-1));
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

void Genotyper::initLikelihoods(void)
{
    // Init map containing the allele likelihoods
    for (auto const &allele : alleleAlphabet_) {
        alleleLikelihoods_.insert(std::pair<char, double>(allele, 0.0));
    }

    // Init map containing the genotype likelihoods
    int chosen[ALLELE_ALPHABET_SIZE];
    combinationsWithRepetitions(&genotypeAlphabet_, alleleAlphabet_, chosen, 0, polyploidy_, 0, ALLELE_ALPHABET_SIZE);

    CALQ_LOG("Initializing genotype alphabet with %zu possible genotypes", genotypeAlphabet_.size());
    for (auto &genotype : genotypeAlphabet_) {
        genotypeLikelihoods_.insert(std::pair<std::string, double>(genotype, 0.0));
    }
}

void Genotyper::resetLikelihoods(void)
{
    for (auto &genotypeLikelihood : genotypeLikelihoods_) {
        genotypeLikelihood.second = 0.0;
    }

    for (auto &alleleLikelihood : alleleLikelihoods_) {
        alleleLikelihood.second = 0.0;
    }
}

void Genotyper::computeGenotypeLikelihoods(const std::string &seqPileup,
                                           const std::string &qualPileup,
                                           const size_t &depth)
{
    resetLikelihoods();

    double *tempGenotypeLikelihoods = (double *)calloc(genotypeAlphabet_.size(), sizeof(double));
    int itr = 0;
    for (size_t d = 0; d < depth; d++) {
        char y = (char)seqPileup[d];
//         double q = (double)(qualPileup[d] - qualMin_);
        double q = (double)(qualPileup[d] - 33);
//         std::cout << "qmin, q: " << qualMin_ << "," << q << std::endl;
        if ((q > qualMax_-qualMin_) || q < 0) {
//             CALQ_LOG("min,curr,max=%d,%d,%d", 0, (int)q, qualMax_-qualMin_);
            throwErrorException("Quality value out of range");
        }
        double pStrike = 1 - pow(10.0, -q/10.0);

        double pError = (1-pStrike) / (ALLELE_ALPHABET_SIZE-1);
        itr = 0;
        for (auto const &genotype : genotypeAlphabet_) {
            double p = 0.0;
            for (int i = 0; i < polyploidy_; i++) {
                p += (y == genotype[i]) ? pStrike : pError;
            }
            p /= polyploidy_;

            // We are using the log likelihood to avoid numerical problems
            tempGenotypeLikelihoods[itr++] += log(p);
        }
    }
    itr = 0;
    for (auto const &genotype : genotypeAlphabet_) {
        genotypeLikelihoods_[genotype] = tempGenotypeLikelihoods[itr++];
    }
    free(tempGenotypeLikelihoods);

    // Normalize the genotype likelihoods
    double cum = 0.0;
    for (auto &genotypeLikelihood : genotypeLikelihoods_) {
        genotypeLikelihood.second = exp(genotypeLikelihood.second);
        cum += genotypeLikelihood.second;
    }
    for (auto &genotypeLikelihood : genotypeLikelihoods_) {
        genotypeLikelihood.second /= cum;
    }
}

} // namespace calq

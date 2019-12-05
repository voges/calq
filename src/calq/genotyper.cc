#include "genotyper.h"
#include <cassert>
#include <cmath>
#include <utility>

namespace calq {

static int combinationsWithRepetitions(std::vector<std::string>* const genotypeAlphabet,
                                       const std::vector<char>& alleleAlphabet, int* const got, const int nChosen,
                                       const int len, const int at, const int maxTypes) {
    if (nChosen == len) {
        if (!got) {
            return 1;
        }
        std::string tmp;
        for (int i = 0; i < len; i++) {
            tmp += alleleAlphabet[got[i]];
        }
        genotypeAlphabet->push_back(tmp);
        return 1;
    }

    int count = 0;
    for (int i = at; i < maxTypes; i++) {
        if (got) {
            got[nChosen] = i;
        }
        count += combinationsWithRepetitions(genotypeAlphabet, alleleAlphabet, got, nChosen + 1, len, i, maxTypes);
    }

    return count;
}

Genotyper::Genotyper(const int polyploidy, const int qualOffset, const int numQuantizers)
    : alleleAlphabet_(ALLELE_ALPHABET),
      alleleLikelihoods_(),
      genotypeAlphabet_(),
      genotypeLikelihoods_(),
      numQuantizers_(numQuantizers),
      polyploidy_(polyploidy),
      qualOffset_(qualOffset) {
    initLikelihoods();
}

int Genotyper::computeQuantizerIndex(const std::string& seqPileup, const std::string& qualPileup) {
    const size_t depth = seqPileup.length();

    assert(depth == qualPileup.length());

    if (depth == 0) {
        return numQuantizers_;  // Computation of quantizer index not possible
    }
    if (depth == 1) {
        return (numQuantizers_ - 1);  // No inference can be made, stay safe
    }

    computeGenotypeLikelihoods(seqPileup, qualPileup, depth);

    double largestGenotypeLikelihood = 0.0;
    double secondLargestGenotypeLikelihood = 0.0;
    for (auto& genotypeLikelihood : genotypeLikelihoods_) {
        if (genotypeLikelihood.second > secondLargestGenotypeLikelihood) {
            secondLargestGenotypeLikelihood = genotypeLikelihood.second;
        }
        if (secondLargestGenotypeLikelihood > largestGenotypeLikelihood) {
            secondLargestGenotypeLikelihood = largestGenotypeLikelihood;
            largestGenotypeLikelihood = genotypeLikelihood.second;
        }
    }

    double confidence = largestGenotypeLikelihood - secondLargestGenotypeLikelihood;

    auto quant = static_cast<int>((1 - confidence) * (numQuantizers_ - 1));

    return quant;
}

void Genotyper::initLikelihoods() {
    // Initialize map containing the allele likelihoods
    for (auto const& allele : alleleAlphabet_) {
        alleleLikelihoods_.insert(std::pair<char, double>(allele, 0.0));
    }

    // Initialize map containing the genotype likelihoods
    int chosen[ALLELE_ALPHABET_SIZE];
    combinationsWithRepetitions(&genotypeAlphabet_, alleleAlphabet_, chosen, 0, polyploidy_, 0,
                                static_cast<int>(ALLELE_ALPHABET_SIZE));

    // Initialize genotype alphabet
    for (auto& genotype : genotypeAlphabet_) {
        genotypeLikelihoods_.insert(std::pair<std::string, double>(genotype, 0.0));
    }
}

void Genotyper::resetLikelihoods() {
    for (auto& genotypeLikelihood : genotypeLikelihoods_) {
        genotypeLikelihood.second = 0.0;
    }

    for (auto& alleleLikelihood : alleleLikelihoods_) {
        alleleLikelihood.second = 0.0;
    }
}

void Genotyper::computeGenotypeLikelihoods(const std::string& seqPileup, const std::string& qualPileup,
                                           const size_t depth) {
    resetLikelihoods();

    auto* tempGenotypeLikelihoods = static_cast<double*>(calloc(genotypeAlphabet_.size(), sizeof(double)));
    int itr = 0;
    for (size_t d = 0; d < depth; d++) {
        auto y = static_cast<char>(seqPileup[d]);
        auto q = static_cast<double>(qualPileup[d] - qualOffset_);

        double pStrike = 1 - pow(10.0, -q / 10.0);
        double pError = (1 - pStrike) / static_cast<double>(ALLELE_ALPHABET_SIZE - 1);

        itr = 0;
        for (auto const& genotype : genotypeAlphabet_) {
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
    for (auto const& genotype : genotypeAlphabet_) {
        genotypeLikelihoods_[genotype] = tempGenotypeLikelihoods[itr++];
    }
    free(tempGenotypeLikelihoods);

    // Normalize the genotype likelihoods
    double cum = 0.0;
    for (auto& genotypeLikelihood : genotypeLikelihoods_) {
        genotypeLikelihood.second = exp(genotypeLikelihood.second);
        cum += genotypeLikelihood.second;
    }
    for (auto& genotypeLikelihood : genotypeLikelihoods_) {
        genotypeLikelihood.second /= cum;
    }
}

const std::map<std::string, double>& Genotyper::getGenotypeLikelihoods(const std::string& seqPileup,
                                                                       const std::string& qualPileup) {
    computeGenotypeLikelihoods(seqPileup, qualPileup, qualPileup.size());
    return genotypeLikelihoods_;
}

}  // namespace calq

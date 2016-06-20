/** @file Genotyper2.h
 *  @brief This file contains the definitions of the Genotyper2 class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: what (who)
 */

#ifndef GENOTYPER2_H
#define GENOTYPER2_H

#include <map>
#include <string>
#include <vector>

/** @brief Class: Genotyper2
 *
 *  The Genotyper class provides methods for computing quantizer indices for
 *  genomic positions given the observed data.
 */
class Genotyper2 {
public:
    Genotyper2(void) {};
    ~Genotyper2(void) {};

    int computeQuantizerIndex(const char &reference,
                              const std::string &observedNucleotides,
                              const std::string &observedQualityValues);
};

#endif // GENOTYPER_H


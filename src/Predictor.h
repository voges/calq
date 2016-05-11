/** @file Predictor.h
 *  @brief This file contains the definition of the Predictor class.
 *  @author Philipp Schelske (schelske)
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  2016-05-10: edited coding style, added const specifiers (voges)
 */

#ifndef PREDICTOR_H
#define PREDICTOR_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

/** @brief Class: Predictor
 *
 *  The Predictor class provides two methods as main interface:
 *  - predict: Memory is used to predict the next value in a sequence.
 *  - update: The current symbol x is used to update the predictor.
 */
class Predictor {
public:
    Predictor(const size_t &alphabetSize, const size_t &memorySize, const size_t &offset);
    ~Predictor(void);

    int predict(const std::vector<int> &memory) const;
    void update(const std::vector<int> &memory, const int &x);

    void reset(void);

    void writeFrequencyTable(std::ostream &os);
    void createCSVFile(void);

private:
    const size_t alphabetSize;
    const int alphabetMax;
    const int alphabetMin;
    std::vector<std::vector<int>> frequencyTable;
    const size_t memorySize;
    size_t numCols;
    size_t numRows;
    const unsigned int offset;

    int computeState(const std::vector<int> &memory) const;
    int findMax(const std::vector<int> &array) const;
    void initFrequencyTable(void);
};

#endif // PREDICTOR_H


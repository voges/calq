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

class Predictor {
public:
    Predictor(const size_t &alphabetSize, const size_t &memorySize);
    ~Predictor(void);

    int predict(const std::vector<int> &memory);
    void update(const std::vector<int> &memory, const int &q);
    unsigned int getOffset(void) const;
    void setOffset(const unsigned int &offset);
    void createCSVFile(void);

private:
    size_t alphabetSize;
    std::vector<std::vector<int>> frequencies;
    size_t memorySize;
    size_t numCols;
    size_t numRows;
    unsigned int offset;

    int findMax(const std::vector<int> &array);
};

#endif // PREDICTOR_H


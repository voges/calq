/** @file Predictor.c
 *  @brief This file contains the implementation of the Predictor class.
 *  @author Philipp Schelske (schelske)
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  2016-05-10: edited coding style; added const specifiers; CSV filename is
 *              now generated from SAM input filename (voges)
 *  2016-05-11: fixed update() and predict() functions to work with arbitrary
 *              memory sizes; added public method writeFrequencyTable() (voges)
 */

#include "Predictor.h"
#include "Exceptions.h"
#include "common.h"
#include <climits>
#include <math.h>

Predictor::Predictor(const size_t &alphabetSize, const size_t &memorySize, const size_t &offset)
    : alphabetSize(alphabetSize)
    , alphabetMax(offset+alphabetSize)
    , alphabetMin(offset)
    , frequencyTable()
    , memorySize(memorySize)
    , numCols(0)
    , numRows(0)
    , offset(offset)
{
    numCols = alphabetSize;
    numRows = (size_t)pow(alphabetSize, memorySize);
    initFrequencyTable();
}

Predictor::~Predictor(void) 
{
    //writeFrequencyTable(std::cout);
}

int Predictor::predict(const std::vector<int> &memory) const
{
    if (memory.size() != memorySize) {
        throwErrorException("Memory sizes do not match");
    }

    if (memory[0] == -1) {
        throwErrorException("Memory is not initialized");
    }

    int row = computeState(memory);

    return findMax(frequencyTable[row]) + offset;
}

void Predictor::update(const std::vector<int> &memory, const int &x)
{
    if (x < alphabetMin || x > alphabetMax) {
        throwErrorException("Symbol out of range");
    }

    if (memory.size() != memorySize) {
        throwErrorException("Memory sizes do not match");
    }

    if (memory[0] == -1) {
        throwErrorException("Memory is not initialized");
    }

    int row = computeState(memory);
    frequencyTable[row][x-offset]++;
}

void Predictor::reset(void)
{
    initFrequencyTable();
}

void Predictor::writeFrequencyTable(std::ostream &os)
{
    for (size_t i = 0; i < numRows; ++i) {
        for (size_t j = 0; j < numCols; ++j) {
            if (j == (alphabetSize - 1)) {
                os << frequencyTable[i][j] << "\n";
            } else {
                os << frequencyTable[i][j] << ", ";
            }
        }
    }
}

void Predictor::createCSVFile(void)
{
    std::string filename = cliOptions.infileName + ".csv";
    std::ofstream ofs;

    if (fileExists(filename) && cliOptions.force == false) {
        std::cout << "CSV file already exists: " << filename << std::endl;
        std::cout << "Do you want to overwrite it? ";
        if (!yesno()) {
            throwUserException("Exited because we do not overwrite output files (add '-f' to force overwriting)");
        }
    }

    ofs.open(filename);
    writeFrequencyTable(ofs);
    ofs.close();

    std::cout << "Wrote prediction table to: " << filename << std::endl;
}

int Predictor::computeState(const std::vector<int> &memory) const
{
    if (memory.size() != memorySize) {
        throwErrorException("Memory sizes do not match");
    }

    // compute state (i.e. row) from memory values
    int state = memory[0] - offset;
    for (size_t i = 1; i < memorySize; ++i) {
        state += (memory[i] - offset) * alphabetSize * i;
    }

    return state;
}

int Predictor::findMax(const std::vector<int> &array) const
{
    int max = INT_MIN;
    int pos = 0;

    for (size_t i = 0; i < array.size(); ++i) {
        if (array[i] > max) {
            max = array[i];
            pos = (int)i;
        }
    }

    return pos;
}

void Predictor::initFrequencyTable(void) 
{
    // resize frequency table and init it with zeros
    frequencyTable.resize(numRows, std::vector<int>(numCols, 0));

    // construct unity matrices
    for (size_t i = 0; i < numRows; ++i) {
        frequencyTable[i][i%alphabetSize] = 1;
    }
}


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
 */

#include "Predictor.h"
#include "Exceptions.h"
#include "common.h"
#include <climits>
#include <math.h>

static const int OFFSET = 33;

Predictor::Predictor(const size_t &alphabetSize, const size_t &memorySize)
    : alphabetSize(alphabetSize)
    , frequencies(pow(alphabetSize, memorySize), std::vector<int>(alphabetSize, 0))
    , memorySize(memorySize)
    , numCols(0)
    , numRows(0)
    , offset(OFFSET)
{
    numCols = alphabetSize;
    numRows = (size_t)pow(alphabetSize, memorySize);

    for (size_t i = 0; i < numRows; ++i) {
        frequencies[i][i%alphabetSize] = 1;
    }
}

Predictor::~Predictor(void) 
{
    // empty
}

int Predictor::predict(const std::vector<int> &memory)
{
    if (memory.size() != memorySize) {
        throwErrorException("Memory sizes do not match");
    }

    //std::cout << "memory: " << std::endl;
    //for (int i = 0; i < memorySize; i++)
    //    std::cout << (char)memory[i] << " ";
    //std::cout << std::endl;

    // check wether memory is initialized (otherwise it is filled with '-1')
    if (memory[0] != -1) {
        // TODO!!
        int val1 = memory[0];
        val1 -= offset;
        val1 *= alphabetSize - 1;
        int val2 = memory[1];
        val2 -= offset;

        int row = val1 + val2;
        return (int)findMax(frequencies[row]) + offset;
    } else {
        return memory[0];
    }
}

void Predictor::update(const std::vector<int> &memory, const int &q)
{
    if (memory.size() != memorySize) {
        throwErrorException("Memory sizes do not match");
    }

    // TODO only update frequency table if memory is initialized
    if (memory[0] != -1) {
        int row = ((int)memory[0]-offset)*40 + (int)memory[1];
        frequencies[row][(int)q-offset] = frequencies[row][(int)q-offset] + 1;
    }
}

unsigned int Predictor::getOffset(void) const
{
    return offset;
}

void Predictor::setOffset(const unsigned int &offset)
{
    this->offset = offset;
}

int Predictor::findMax(const std::vector<int> &array)
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

void Predictor::createCSVFile(void)
{
    std::string filename = cliOptions.infileName + ".csv";
    std::ofstream ofs;

    if (fileExists(filename) && cliOptions.force == false) {
        std::cout << "CSV file already exists: " << filename << std::endl;
        std::cout << "Do you want to overwrite it? ";
        if (!yesno()) {
            throwUserException("Exited because we do not overwrite the output file");
        }
    }

    ofs.open(filename);

    for (size_t i = 0; i < numRows; ++i) {
        for (size_t j = 0; j < numCols; ++j) {
            if (j == (alphabetSize-1)) {
                ofs << frequencies[i][j]<<"\n";
            } else {
                ofs << frequencies[i][j]<< ", ";
            }
        }
    }

    ofs.close();
}


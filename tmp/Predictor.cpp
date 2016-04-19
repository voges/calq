#include "qualcodec/Predictor.hpp"
#include <iostream>

Predictor::Predictor(size_t memorySize, size_t alphabetSize)
{
    this->memorySize = memorySize;
    this->alphabetSize = alphabetSize;
    this->frequencies.resize(this->memorySize, this->alphabetSize);
}

Predictor::~Predictor()
{

}

void Predictor::predict(std::string qual)
{

}


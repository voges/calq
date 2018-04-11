#include <cmath>
#include <algorithm>

#include "GaussKernel.h"

GaussKernel::GaussKernel(double sigma) : SIGMA(sigma), INV_SQRT_SIGMA_2PI(1.0/(std::sqrt(2.0*PI)*sigma)){
}
    
double GaussKernel::calcValue(size_t pos, size_t size){
    const double MEAN = (size-1)/2;
    double exponent = (pos - MEAN)/SIGMA;
    exponent = exponent * exponent * (-0.5);
    return INV_SQRT_SIGMA_2PI * std::pow(EULER, exponent);
}

size_t GaussKernel::calcMinSize(double threshold, size_t maximum){
    threshold /= INV_SQRT_SIGMA_2PI;
    threshold = std::log(threshold);
    threshold *= -2.0;
    threshold = std::sqrt(threshold) * SIGMA; // Euler now reversed
        
    size_t size = std::ceil(threshold)*2+1; // + 1 to make sure it is odd.
    return std::min(size, maximum);        
}


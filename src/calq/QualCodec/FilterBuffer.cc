#include "FilterBuffer.h"
#include "Common/Exceptions.h"

#include <cmath>
#include <algorithm>

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



//New activity score in pipeline
void FilterBuffer::push (double activityScore) {
    buffer[bufferPos] = activityScore;
    bufferPos = (bufferPos+1) % buffer.size();
}

//Calculate filter score at offset position
double FilterBuffer::filter() const{
    double result = 0.0;
    for(size_t i=0;i<kernel.size();++i) {
        size_t index = (i+ bufferPos) % buffer.size();
        result += kernel[i] * buffer[index];

    }
    return result;
}

//Initialize buffer and 
FilterBuffer::FilterBuffer(const std::function<double(size_t, size_t)>& kernelBuilder, size_t kernelSize){
    if(!(kernelSize % 2)) {
        throwErrorException("Kernel size must be an odd number");
    }
    bufferPos = 0;
    buffer.resize(kernelSize, 0.0);
    kernel.resize(kernelSize, 0.0);

    for(size_t i=0;i<kernel.size();++i){
        kernel[i] = kernelBuilder(i, kernelSize);
    }
}

FilterBuffer::FilterBuffer(){
    bufferPos = 0;
}

size_t FilterBuffer::getSize() const{
    return buffer.size();
} 

size_t FilterBuffer::getOffset() const{
    return (buffer.size()-1)/2;
}

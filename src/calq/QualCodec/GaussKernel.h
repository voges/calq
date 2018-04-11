#include "stddef.h"

class GaussKernel {
private: 
    const double SIGMA;

    constexpr static double PI = 3.14159265359;
    constexpr static double EULER = 2.71828182846;
    const double INV_SQRT_SIGMA_2PI;

public: 
    //init
    GaussKernel(double sigma=1.0);
    
    //Get gauss value at position pos and buffersize size with mean=size/2
    double calcValue(size_t pos, size_t size);

    //Calcutes how big a buffer must be to contain all values above threshold. No size greater than maximum is returned.
    size_t calcMinSize(double threshold, size_t maximum = 128);
};


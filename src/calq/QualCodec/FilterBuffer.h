#include <vector>
#include <stddef.h>
#include <functional>


//----------------------------------------------------------------------------------------------------------------------

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

//----------------------------------------------------------------------------------------------------------------------

class FilterBuffer {
private:
    std::vector<double> kernel;
    std::vector<double> buffer;
    size_t bufferPos; //Index pointer to oldest value

public:
    //New activity score in pipeline
    void push (double activityScore);

    //Calculate filter score at offset position
    double filter() const;

    //Initialize buffer and 
    FilterBuffer(const std::function<double(size_t, size_t)>& kernelBuilder, size_t kernelSize);

    //Create dummy buffer
    FilterBuffer();

    //Buffer size
    size_t getSize() const;

    //DIstance between buffer center and borders
    size_t getOffset() const;
};

//----------------------------------------------------------------------------------------------------------------------

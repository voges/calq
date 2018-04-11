#include <functional>
#include <cmath>

#include "FilterBuffer.h"
#include "GaussKernel.h"

class Haplotyper {
private:
    //Saving kernel and old raw quality scores
    FilterBuffer buffer;

    GaussKernel kernel;
    
    //Function to calulate raw quality score from seqPile, qualPile, position
    std::function<double(std::string, std::string, size_t pos)> scoreCalc;

public:

    //Init
    Haplotyper(const std::function<double(std::string, std::string, size_t pos)>& calcScore);

    //Returns offset between activity scores' position and front
    size_t getOffset() const;

    //Pushes new activity score calculated using parameters and returns filtered acticityscore for (pos-offset)
    double insertData(const std::string& seqPile, const std::string& qualPile, size_t pos);
};

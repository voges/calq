#include "Haplotyper.h"

static constexpr double SIGMA = 17.0;

//Returns score, takes seq and qual pileup and position
Haplotyper::Haplotyper(const std::function<double(std::string, std::string, size_t pos)>& calcScore) : scoreCalc(calcScore), kernel(SIGMA){
    double THRESHOLD = 0.001;
    size_t size = this->kernel.calcMinSize(THRESHOLD);
 
    buffer = FilterBuffer([this](size_t pos, size_t size) -> double{return this->kernel.calcValue(pos,size);}, size);
}

size_t Haplotyper::getOffset() const {
    return buffer.getOffset();
}

double Haplotyper::insertData(const std::string& seqPile, const std::string& qualPile, size_t pos){
    buffer.push(scoreCalc(seqPile, qualPile, pos));
    return buffer.filter();
}

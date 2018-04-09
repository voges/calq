//Returns score, takes seq and qual pileup and position
Haplotyper::Haplotyper(const std::function<double, std::string, std::string, size_t pos>& calcScore) : scoreCalc(calcScore){
    double SIGMA = 17.0;
    double THRESHOLD = 0.001;
    GaussKernel kernel = GaussKernel(SIGMA);
    size_t size = kernel.calqMinSize(THRESHOLD)
 
    buffer = FilterBuffer(kernel.calcValue, size);       
}

size_t Haplotyper::getOffset() const {
    return buffer.getOffset();
}

double Haplotyper::insertData(const std::string& seqPile, const std::string& qualPile, size_t pos){
    buffer.push(coreCalc(seqPile, qualPile, pos));
    return buffer.filter();
}

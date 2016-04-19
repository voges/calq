#ifndef CQ_PREDICTOR_HPP
#define CQ_PREDICTOR_HPP

#include "Eigen/Dense"
using namespace Eigen;

class Predictor
{
  private:
    size_t memorySize;
    size_t alphabetSize;
    Matrix<size_t, Dynamic, Dynamic> frequencies;
  public:
    Predictor();
    Predictor(size_t memorySize, size_t alphabetSize);
    ~Predictor();

    void predict(std::string qual);
    void reset(void);
};

#endif // CQ_PREDICTOR_HPP


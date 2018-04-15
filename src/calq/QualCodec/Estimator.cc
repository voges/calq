#include "Estimator.h"


void Estimator::push(double obs) {
    this->num += 1;
    this->mean += (obs - this->mean) / this->num;
}

double Estimator::getMean() {
    return this->mean;
}

#include <stddef.h>

class Estimator {
private:
    double mean;
    size_t num;
public:
    Estimator();

    void push(double obs);

    double getMean();
};

#include "FilterBuffer.h"

#include <iostream>


void equals(double a, double b, double EPSILON = 0.0001){
        if ((a-EPSILON) < b && (a+EPSILON) > b) {
            std::cout << "Test passed, " << a << " equals " << b << " !" << std::endl;
        } else {
            std::cout << "Test failed, " << a << " not equal to " << b << " !" << std::endl;
    }
}


void gaussKernelTest () {
        //Test variance 1
        GaussKernel k(1.0);

        double buffer[63];
        for(size_t i=0;i<63;++i){
            buffer[i] = k.calcValue(i,63);
        }

        equals(buffer[31], 0.3989);
        equals(buffer[32], 0.2419);
        equals(buffer[30], 0.2419);
        equals(buffer[33], 0.0539);
        equals(buffer[29], 0.0539);

        //Test variance 17
        GaussKernel k2(17.0);
        double buffer2[127];
        for(size_t i=0;i<127;++i){
            buffer2[i] = k2.calcValue(i,127);
        }

        equals(buffer2[63], 0.0234);
        equals(buffer2[62], 0.0234);
        equals(buffer2[64], 0.0234);
        equals(buffer2[61], 0.0233);
        equals(buffer2[65], 0.0233);
        equals(buffer2[41], 0.0101);
        equals(buffer2[85], 0.0101);

        //Test length limits
        equals(k2.calcMinSize(0.01),47);
        equals(k.calcMinSize(0.01),7);
}

void filterBufferTest(){
    FilterBuffer buffer([](size_t pos, size_t size)->double{return pos/(double)(size-1);}, 3);

    //Test "Convolution"
    equals(buffer.filter(), 0);
    buffer.push(1);
    equals(buffer.filter(), 1);
    buffer.push(2);
    equals(buffer.filter(), 2.5);
    buffer.push(3);
    equals(buffer.filter(), 4);

    //Test offset
    equals(buffer.getOffset(), 1);
}


void haplotyperTest(){
    std::cout << "****Starting haplotyper test suite****" << std::endl;
    std::cout << "-> Gauss kernel tests" << std::endl;
    gaussKernelTest();
    std::cout << "-> Filter buffer tests" << std::endl;
    filterBufferTest();
    std::cout << "****Haplotyper test suite finished****" << std::endl;
}



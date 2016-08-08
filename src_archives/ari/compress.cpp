/*
The MIT License (MIT)

Copyright (c) 2014 Mark Thomas Nelson

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.

 This code was written to illustrate the article:
 Data Compression With Arithmetic Coding
 by Mark Nelson
 published at: http://marknelson.us/2014/10/19/data-compression-with-arithmetic-coding

*/
#include <iostream>
#include <fstream>

#include "modelA.h"
#include "compressor.h"

int main(int argc, char* argv[])
{
  if ( argc < 3 ) {
    std::cerr << "missing command line arguments\n";
    return 255;
  }
  try {
    std::ifstream input(argv[1], std::ifstream::binary);
    std::ofstream output(argv[2], std::ofstream::binary);
    modelA<int, 16, 14> cmodel;
    cmodel.dump("cmodel", std::clog);

    std::cout << "compressing...\n";
    compress(input, output, cmodel);
    std::cout << cmodel.m_bytesProcessed << "\n";
    return 0;
  }
  catch (std::exception &ex)
  {
    std::cerr << "Failed with exception: " << ex.what() << "\n";
  }
  return 255;
}

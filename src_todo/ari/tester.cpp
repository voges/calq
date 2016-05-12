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
#include <cstdio>
#include <string>

#include "modelA.h"
#include "compressor.h"
#include "decompressor.h"

struct my_model : public modelA<int>
{
  void pacify(){}
  void frozen(){}
};

int validate( const std::string &input_file,
              const std::string &compressed_file,
              const std::string &output_file,
              double &bpb );

int main(int argc, char* argv[])
{
  try {
    std::cout << "compressing " << argv[1] << "... " << std::flush;
    std::ofstream output1("temp.ari", std::ofstream::binary);
    std::ifstream input1(argv[1], std::ifstream::binary);
    my_model cmodel;

    compress(input1, output1, cmodel);
    output1.close();

    std::ifstream input2("temp.ari", std::ifstream::binary);
    std::ofstream output2("temp.out",std::ofstream::binary);
    cmodel.reset();
    decompress(input2, output2, cmodel );
    output2.close();

    double bpb;
    int result = validate(argv[1], "temp.ari", "temp.out", bpb);
    std::cout << bpb << "\n";
    return 0;
  }
  catch (std::exception &ex)
  {
    std::cerr << "Failed with exception: " << ex.what() << "\n";
  }
  return 255;
}

int validate( const std::string &input_file,
              const std::string &compressed_file,
              const std::string &output_file,
              double &bpb )
{
  bool verbose = false; //might turn this on in some contexts

  std::ifstream in(input_file.c_str(), std::ifstream::binary);
  if ( !in ) {
    std::cout << "validate error opening inptut file: " << input_file << "\n";
    return 255;
  }
  std::ifstream compressed(compressed_file.c_str(), std::ifstream::binary);
  if ( !compressed ) {
    std::cout << "validate error opening compressed file: " << compressed_file << "\n";
    return 255;
  }
  std::ifstream out(output_file.c_str(), std::ifstream::binary);
  if ( !out ) {
    std::cout << "validate error opening output file: " << output_file << "\n";
    return 255;
  }
  in.seekg(0,std::ios::end);
  out.seekg(0,std::ios::end);
  compressed.seekg(0,std::ios::end);
  auto in_length = in.tellg();
  auto out_length = out.tellg();
  auto compressed_length = compressed.tellg();
  in.seekg(0,std::ios::beg);
  out.seekg(0,std::ios::beg);
  if ( verbose )
    std::cout << "input length: " << in_length << "\n"
              << "output length: " << out_length << "\n"
              << "compressed length: " << compressed_length << "\n";
  if ( in_length != out_length ) {
    std::cout << "Error, input file and output file have different lengths\n";
    return 255;
  }
  if ( static_cast<long long>(in_length) == 0 )
    bpb = 8.0;
  else
    bpb = compressed_length * 8.0 / in_length;
  if ( verbose )
    std::cout << "Compressed to " << bpb << " bits per byte\n";
  int c1;
  int c2;
  while ( c1 = in.get(), c2 = out.get(), c1 != -1 || c2 != -1 ) {
    if ( c1 != c2 ) {
      int i = 1;
      std::cerr << "Error comparing at position: " << (std::streamoff(in.tellg()) - 1) << "\n";
      return 255;
    }
  }
  if ( verbose )
    std::cout << "Comparison passed!\n";
  return 0;
}

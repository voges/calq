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
#include <iomanip>
#include <string>

struct {
  static std::pair<double,double> getProbability( char c )
  {
    if (c >= 'A' && c <= 'Z')
      return std::make_pair( (c - 'A') * .01, (c - 'A') * .01 + .01);
    else if (c >= 'a' && c <= 'z')
      return std::make_pair( (c - 'a') * .02 + 0.30, (c - 'a') * .02 + .02 + 0.30);
    else
      throw "character out of range";
  }
  static char getSymbol( double d)
  {
    if ( d >= 0.0 && d < 0.26)
      return 'A' + static_cast<int>(d*100);
    else if ( d >= 0.3 && d < 0.82)
      return 'a' + static_cast<int>((d-0.3)*50);
    else
      throw "message out of range";
  }
} model;

double compress( std::string s)
{
  double high = 1.0;
  double low = 0.0;
  for ( char c : s ) {
    std::pair<double,double> p = model.getProbability(c);
    double range = high - low;
    high = low + range * p.second;
    low = low + range * p.first; 
  }
  return low + (high-low)/2;
}

std::string decompress(double message)
{
  std::string result;
  double high = 1.0;
  double low = 0.0;
  for ( ; ; ) 
  {
    double range = high - low;
    char c = model.getSymbol((message - low)/range);
    result += c;
    if ( c == 'Z' )
      return result;
    std::pair<double,double> p = model.getProbability(c);
    high = low + range * p.second;
    low = low + range * p.first; 
  }
}

int main(int, char **)
{
  double val = compress("WXYZ");
  //double val = compress("abcdWXYZ"); //this will be more interesting
  std::cout << "Compressed value: " 
            << std::setprecision(15)
            << val << "\n";
  std::string result = decompress(val);
  std::cout << "Decompressed result: " << result << "\n";
  return 0;
}
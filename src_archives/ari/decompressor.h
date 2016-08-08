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
#ifndef DECOMPESSOR_DOT_H
#define DECOMPESSOR_DOT_H

#ifdef LOG
#include <iostream>
#include <iomanip>
#endif

#include "byteio.h"
#include "bitio.h"

//
// The arithmetic decompressor is a general purpose decompressor that
// is parameterized on the types of the input, output, and 
// model objects, in an attempt to make it as flexible as
// possible. It is easiest to use by calling the compress()
// convenience function found at the bottom of this header file
//
// The INPUT class is expected to provide a get_bit() function,
// while the output function is expected to provider a put_byte()
// function. Both of these functions should throw exceptions on
// errors. We expect the EOF to be embedded in the compressed
// stream, so it needs to be extracted by the decoder. If the
// compression goes awry, the get_bit() function will be 
// repeatedly called on EOF(), in which case it would be good
// for it to return an error.
//
template<typename INPUT, typename OUTPUT, typename MODEL>
class decompressor
{
  typedef typename MODEL::CODE_VALUE CODE_VALUE;
  typedef typename MODEL::prob prob;
public :
  decompressor(INPUT &input, OUTPUT &output, MODEL &model ) 
  : m_input(input),
    m_output(output),
    m_model(model)
  {
  }
  int operator()()
  {
#ifdef LOG
    std::ofstream log("decompressor.log");
    log << std::hex;
#endif
    CODE_VALUE high = MODEL::MAX_CODE;
    CODE_VALUE low = 0;
    CODE_VALUE value = 0;
    for ( int i = 0 ; i < MODEL::CODE_VALUE_BITS ; i++ ) {
      value <<= 1;
      value += m_input.get_bit() ? 1 : 0;
    }
    for ( ; ; ) {
      CODE_VALUE range = high - low + 1;
      CODE_VALUE scaled_value =  ((value - low + 1) * m_model.getCount() - 1 ) / range;
      int c;
      prob p = m_model.getChar( scaled_value, c );
      if ( c == 256 )
        break;
      m_output.putByte(c);
#ifdef LOG
      log << std::hex << "0x" << std::setw(2) << std::setfill('0') << c;
      if ( c > 0x20 && c <= 0x7f )
        log << "(" << char(c) << ")";
      log << " 0x" << low << " 0x" << high << " => ";
#endif 
      high = low + (range*p.high)/p.count -1;
      low = low + (range*p.low)/p.count;
#ifdef LOG
      log << "0x" << low << " 0x" << high << "\n";
#endif
      for( ; ; ) {
        if ( high < MODEL::ONE_HALF ) {
          //do nothing, bit is a zero
        } else if ( low >= MODEL::ONE_HALF ) {
          value -= MODEL::ONE_HALF;  //subtract one half from all three code values
          low -= MODEL::ONE_HALF;
          high -= MODEL::ONE_HALF;
        } else if ( low >= MODEL::ONE_FOURTH && high < MODEL::THREE_FOURTHS ) {
          value -= MODEL::ONE_FOURTH;
          low -= MODEL::ONE_FOURTH;
          high -= MODEL::ONE_FOURTH;
        } else
          break;
        low <<= 1;
        high <<= 1;
        high++;
        value <<= 1;
        value += m_input.get_bit() ? 1 : 0;
      }
    }
#ifdef LOG
      log << std::hex << "0x" << std::setw(2) << std::setfill('0') << 256;
      log << " 0x" << low << " 0x" << high << "\n";
#endif 
    return 0;
  }
private :
  OUTPUT &m_output;
  INPUT &m_input;
  MODEL &m_model;
};

//
// This convenience function takes care of
// constructing the decompressor and the
// input and output objects, then calling
// the decompressor.
//
template<typename INPUT, typename OUTPUT, typename MODEL>
int decompress(INPUT &source, OUTPUT &target, MODEL &model)
{
  input_bits<INPUT> in(source,MODEL::CODE_VALUE_BITS);
  output_bytes<OUTPUT> out(target);
  decompressor<input_bits<INPUT>, output_bytes<OUTPUT>, MODEL> d(in,out, model);
  return d();
}

#endif //#ifndef DECOMPESSOR_DOT_H

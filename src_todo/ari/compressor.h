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
#ifndef COMPRESSOR_DOT_H
#define COMPRESSOR_DOT_H

#include <stdexcept>

#include "byteio.h"
#include "bitio.h"

#ifdef LOG
#include <iomanip>
#include <iostream>
#endif

//
// The arithmetic compressor is a general purpose compressor that
// is parameterized on the types of the input, output, and 
// model objects, in an attempt to make it as flexible as
// possible. It is easiest to use by calling the compress()
// convenience function found at the bottom of this header file
//

template<typename INPUT, typename OUTPUT, typename MODEL>
class compressor
{
  typedef typename MODEL::CODE_VALUE CODE_VALUE;
  typedef typename MODEL::prob prob;
public :
  compressor(INPUT &input, OUTPUT &output, MODEL &model ) 
  : m_input(input),
    m_output(output),
    m_model(model)
  {
  }
  int operator()()
  {
#ifdef LOG
    std::ofstream log("compressor.log");
    log << std::hex;
#endif
    int pending_bits = 0;
    CODE_VALUE low = 0;
    CODE_VALUE high = MODEL::MAX_CODE;
    for ( ; ; ) {
      int c = m_input.getByte();
      if ( c == -1 )
        c = 256;
#ifdef LOG
      log << std::hex << "0x" << std::setw(2) << std::setfill('0') << c;
      if ( c > 0x20 && c <= 0x7f )
        log << "(" << char(c) << ")";
      log << " 0x" << low << " 0x" << high << " => ";
#endif
      prob p = m_model.getProbability( c );
      CODE_VALUE range = high - low + 1;
      high = low + (range * p.high / p.count) - 1;
      low = low + (range * p.low / p.count);
#ifdef LOG
      log << "0x" << low << " 0x" << high << "\n";
#endif
      //
      // On each pass there are six possible configurations of high/low,
      // each of which has its own set of actions. When high or low
      // is converging, we output their MSB and upshift high and low.
      // When they are in a near-convergent state, we upshift over the
      // next-to-MSB, increment the pending count, leave the MSB intact,
      // and don't output anything. If we are not converging, we do
      // no shifting and no output.
      // high: 0xxx, low anything : converging (output 0)
      // low: 1xxx, high anything : converging (output 1)
      // high: 10xxx, low: 01xxx : near converging
      // high: 11xxx, low: 01xxx : not converging
      // high: 11xxx, low: 00xxx : not converging
      // high: 10xxx, low: 00xxx : not converging
      //
      for ( ; ; ) {
        if ( high < MODEL::ONE_HALF )
          put_bit_plus_pending(0, pending_bits);
        else if ( low >= MODEL::ONE_HALF )
          put_bit_plus_pending(1, pending_bits);
        else if ( low >= MODEL::ONE_FOURTH && high < MODEL::THREE_FOURTHS ) {
          pending_bits++;
          low -= MODEL::ONE_FOURTH;  
          high -= MODEL::ONE_FOURTH;  
        } else
          break;
        high <<= 1;
        high++;
        low <<= 1;
        high &= MODEL::MAX_CODE;
        low &= MODEL::MAX_CODE;
      }
      if ( c == 256 ) //256 is the special EOF code
        break;
    }
    pending_bits++;
    if ( low < MODEL::ONE_FOURTH )
      put_bit_plus_pending(0, pending_bits);
    else
      put_bit_plus_pending(1, pending_bits);
#ifdef LOG
    log.close();
#endif
    return 0;
  }

  inline void put_bit_plus_pending(bool bit, int &pending_bits)
  {
    m_output.put_bit(bit);
    for ( int i = 0 ; i < pending_bits ; i++ )
      m_output.put_bit(!bit);
    pending_bits = 0;
  }
private :
  OUTPUT &m_output;
  INPUT &m_input;
  MODEL &m_model;
};

//
// This convenience function takes care of
// constructing the compressor and the
// input and output objects, then calling
// the compressor. Letting the user of the class
// call a template function instead of instantating
// the template class object eases syntax
// requirements a bit.
//
template<typename INPUT, typename OUTPUT, typename MODEL>
int compress(INPUT &source, OUTPUT &target, MODEL &model)
{
  input_bytes<INPUT> in(source);
  output_bits<OUTPUT> out(target);
  compressor<input_bytes<INPUT>, output_bits<OUTPUT>, MODEL> c(in,out, model);
  return c();
}

#endif //#ifndef COMPRESSOR_DOT_H

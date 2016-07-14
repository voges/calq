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
#ifndef MODEL_A_DOT_H
#define MODEL_A_DOT_H

#include <iostream>
#include <stdexcept>
#include "model_metrics.h"

template<typename CODE_VALUE_ = unsigned int, 
         int CODE_VALUE_BITS_ = (std::numeric_limits<CODE_VALUE_>::digits + 3) / 2,
         int FREQUENCY_BITS_ = std::numeric_limits<CODE_VALUE_>::digits - CODE_VALUE_BITS_>
struct modelA : public model_metrics<CODE_VALUE_, CODE_VALUE_BITS_, FREQUENCY_BITS_> 
{
  typedef model_metrics<CODE_VALUE_, CODE_VALUE_BITS_, FREQUENCY_BITS_> metrics;
  typedef typename metrics::prob prob;
  typedef CODE_VALUE_ CODE_VALUE;
  using metrics::MAX_CODE;
  using metrics::MAX_FREQ;
  using metrics::CODE_VALUE_BITS;
  using metrics::ONE_FOURTH;
  using metrics::ONE_HALF;
  using metrics::THREE_FOURTHS;
  //
  // variables used by the model
  //
  CODE_VALUE cumulative_frequency[258]; //Character a is defined by the range cumulative_frequency[a],
                                        //cumulative_frequency[a+1], with cumulative_frequency[257]
                                        //containing the total count for the model. Note that entry
                                        //256 is the EOF.
  unsigned long long m_bytesProcessed;
  static_assert( MAX_FREQ > 257, "Not enough code bits to represent the needed symbol library" );

  modelA()
  {
    reset();
  }
  void reset()
  {
    for ( int i = 0 ; i < 258 ; i++ )
      cumulative_frequency[i] = i;
    m_bytesProcessed = 0;
    m_frozen = false;
  }
  virtual inline void pacify()
  {
    if ( (++m_bytesProcessed % 1000) == 0 )
      std::cout << m_bytesProcessed << "\r";
  }
  virtual void frozen()
  {
    std::cout << "Frozen at: " << m_bytesProcessed << "\n";
  }
  void inline update(int c)
  {
    for ( int i = c + 1 ; i < 258 ; i++ )
      cumulative_frequency[i]++;
    if ( cumulative_frequency[257] >= MAX_FREQ ) {
      m_frozen = true;
      frozen();
    }
  }
  prob getProbability(int c)
  {
    prob p = { cumulative_frequency[c], cumulative_frequency[c+1], cumulative_frequency[257] };
    if ( !m_frozen ) 
      update(c);
    pacify();
    return p;
  }
  prob getChar(CODE_VALUE scaled_value, int &c)
  {
    pacify();
    for ( int i = 0 ; i < 257 ; i++ )
      if ( scaled_value < cumulative_frequency[i+1] ) {
        c = i;
        prob p = {cumulative_frequency[i], cumulative_frequency[i+1],cumulative_frequency[257]};
        if ( !m_frozen)
          update(c);
        return p;
      }
      throw std::logic_error("error");
  }
  CODE_VALUE getCount()
  {
    return cumulative_frequency[257];
  }
  bool m_frozen;

};

#endif //#ifndef MODEL_A_DOT_H

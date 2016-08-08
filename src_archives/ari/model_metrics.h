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
#ifndef MODEL_METRICS_DOT_H
#define MODEL_METRICS_DOT_H

#include <typeinfo>
#include <limits>
#include <stdint.h>

//
// By default we set this whole thing up so that 
// all math is done IN CODE_VALUE integers, which
// are usually 32 bit ints or longs. In these
// you might set things up so that your max
// frequency or your max code value fit into
// a smaller size integer, say 16 bits. You can
// possibly get some efficiency then by storing
// counts in a smaller type of integer, but at this
// time we are not going to try to implement that.
// If we were going to do so automatically, this 
// template taken from this stackoverflow convo:
// http://stackoverflow.com/questions/12082571/how-to-figure-out-the-smallest-integral-type-that-can-represent-a-number-in-com
// 

template<typename CODE_VALUE, int CODE_VALUE_BITS_, int FREQUENCY_BITS_>
struct model_metrics {

  static const int PRECISION            = std::numeric_limits<CODE_VALUE>::digits;
  static const int CODE_VALUE_BITS      = CODE_VALUE_BITS_;
  static const int FREQUENCY_BITS       = FREQUENCY_BITS_;
  static const CODE_VALUE MAX_CODE      = (CODE_VALUE(1) << CODE_VALUE_BITS) - 1;
  static const CODE_VALUE MAX_FREQ      = (CODE_VALUE(1) << FREQUENCY_BITS) - 1;
  static const CODE_VALUE ONE_FOURTH    = CODE_VALUE(1) << (CODE_VALUE_BITS - 2);;
  static const CODE_VALUE ONE_HALF      = 2 * ONE_FOURTH;
  static const CODE_VALUE THREE_FOURTHS = 3 * ONE_FOURTH;

  static_assert( std::numeric_limits<CODE_VALUE>::digits >= CODE_VALUE_BITS, 
                 "CODE_VALUE_BITS is too large to fit in a CODE_VALUE type" );
  static_assert( FREQUENCY_BITS <= (CODE_VALUE_BITS + 2 ),
                 "FREQUENCY_BITS can be no greater than CODE_VALUE_BITS - 2" );
  static_assert( (CODE_VALUE_BITS + FREQUENCY_BITS) <= PRECISION,
                 "CODE_VALUE_BITS + FREQUENCY_BITS cannot exceed precision of CODE_VALUE" );

  template<typename STRING, typename STREAM>
  static void dump(const STRING &name, STREAM &s)
  {
    s << "Model " << name << " created with:\n"
      << "CODE_VALUE of type "   << typeid(CODE_VALUE).name() << " with " << PRECISION << " bits\n"
      << "CODE_VALUE_BITS " << CODE_VALUE_BITS << " bits giving MAX_CODE of "      << MAX_CODE << "\n"
      << "FREQUENCY_BITS "  << FREQUENCY_BITS  << " bits giving MAX_FREQUENCY of " << MAX_FREQ << "\n"
      << "MAX_CODE: "       << MAX_CODE      << " (0x" << std::hex << MAX_CODE      << std::dec << ")\n"
      << "MAX_FREQ: "       << MAX_FREQ      << " (0x" << std::hex << MAX_FREQ      << std::dec << ")\n"
      << "ONE_FOURTH: "     << ONE_FOURTH    << " (0x" << std::hex << ONE_FOURTH    << std::dec << ")\n"
      << "ONE_HALF: "       << ONE_HALF      << " (0x" << std::hex << ONE_HALF      << std::dec << ")\n"
      << "THREE_FOURTHS: "  << THREE_FOURTHS << " (0x" << std::hex << THREE_FOURTHS << std::dec << ")\n";
  }
  struct prob { 
    CODE_VALUE low;
    CODE_VALUE high; 
    CODE_VALUE count;
  };
};

#endif //#ifndef MODEL_METRICS_DOT_H

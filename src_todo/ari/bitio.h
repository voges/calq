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
#ifndef BITIO_DOT_H
#define BITIO_DOT_H

#include "byteio.h"

template<typename OUTPUT>
class output_bits
{
public :
  output_bits(OUTPUT &output)
  : m_Output(output),
    m_NextByte(0),
    m_Mask(0x80) 
  {
  }
  ~output_bits()
  {
    if ( m_Mask != 0x80 )
      m_Output.putByte(m_NextByte);
  }
  void put_bit( bool val ) {
    if ( val )
      m_NextByte |= m_Mask;
    m_Mask >>= 1;
    if ( !m_Mask) {
      m_Output.putByte(m_NextByte);
      m_Mask = 0x80;
      m_NextByte = 0;
    }
  }
private :
  output_bytes<OUTPUT> m_Output;
  char m_NextByte;
  unsigned char m_Mask;

};

template<typename INPUT>
class input_bits
{
public :
  input_bits(INPUT &input, int code_value_bits)
  : m_Input(input),
    m_CurrentByte(0),
    m_LastMask(1),
    m_CodeValueBits(code_value_bits) 
  {
  }
  bool get_bit() {
    if ( m_LastMask == 1 ) {
      m_CurrentByte = m_Input.getByte();
      if ( m_CurrentByte < 0 ) {
        if ( m_CodeValueBits <= 0 ) 
          throw std::logic_error("EOF on input");
        else
          m_CodeValueBits -= 8;
      }
      m_LastMask = 0x80;
    } else
      m_LastMask >>= 1;
    return (m_CurrentByte & m_LastMask) != 0;
  }

private :
  input_bytes<INPUT> m_Input;
  int m_CurrentByte;
  unsigned char m_LastMask;
  int m_CodeValueBits;
};

#endif //#ifndef BITIO_DOT_H


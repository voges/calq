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

#include "bitstream.h"
#include "definitions.h"


/* Why is the bitstream.h class not sufficient for the use in this project and why do we need this wrapper class:
 * obistream starts to construct a byte from the LSB to the MSB. The ARI-algorithm works the exact opposite way and expects 
 * the MSB to be read first. This is were ibitstream starts to fail and
 *
 */
template<typename OUTPUT>
class output_bits
{
public :
  output_bits(obitstream &output)
  : m_NextByte(0),
    m_Mask(0x80) 
  {
      m_Output = &output;
  }
  ~output_bits()
  {
    if ( m_Mask != 0x80 )
      m_Output->writeByte(m_NextByte);
  }
    
  void put_bit( bool val ) {
    if ( val )
      m_NextByte |= m_Mask;
    m_Mask >>= 1;
    if ( !m_Mask) {
      m_Output->writeByte(m_NextByte);
      m_Mask = 0x80;
      m_NextByte = 0;
    }
  }
    /*
     * Added by Philipp Schelske
     * The following functions are added, so we can access basic functionality of the streamclass.
     * The output_bytes are only allowd to be derived from std::ostream.
     */
    int finishBlock(){
        if ( m_Mask != 0x80 ){
            m_Output->writeByte(m_NextByte);
            m_Mask = 0x80;
            m_NextByte = 0x0;
            return 1;
        }
        m_NextByte = 0x0;
        m_Mask = 0x80;

        return 0;
    }
    std::streampos tellp(){
        return m_Output->tellp();
    }
    void seekp(std::streampos pos, std::ios_base::seekdir way){
        m_Output->seekp(pos, way);
    }
    void seekp(std::streampos pos){
        m_Output->seekp(pos);
    }
    void writeUint64(const uint64_t x){
        m_Output->writeUint64(x);
    }
    void writeByte(const BYTE byte){
        m_Output->writeByte(byte);
    }
    
private :
  obitstream *m_Output;
  char m_NextByte;
  unsigned char m_Mask;

};

template<typename INPUT>
class input_bits
{
public :
  input_bits(ibitstream &input, int code_value_bits)
  : m_CurrentByte(0),
    m_LastMask(1),
    m_CodeValueBits(code_value_bits) 
  {
      m_Input = &input;
  }
    
  bool get_bit() {
    if ( m_LastMask == 1 ) {
      m_Input->readByte(m_CureB);
        m_CurrentByte = m_CureB;
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
    /*
     * Added by Philipp Schelske
     * Following functions are implemented to make the ibitstream functions accessable.
     */
    bool eof(){
        return m_Input->eof();
    }
    void seekg(std::streampos pos, std::ios_base::seekdir way){
        m_Input->seekg(pos, way);
    }
    void seekg(std::streampos pos){
        m_Input->seekg(pos);
    }

    int readUint64(uint64_t &x){
       return m_Input->readUint64(x);
    }
    int readByte(BYTE &byte){
        return m_Input->readByte(byte);
    }
    std::streampos tellg(){
        return m_Input->tellg();
    }
    
    void reset(){
        m_CurrentByte = 0x0;
        m_LastMask = 1;
        m_CureB = 0x0;
    }


private :
  ibitstream *m_Input;
    BYTE m_CureB;
  int m_CurrentByte;
  unsigned char m_LastMask;
  int m_CodeValueBits;
};

#endif //#ifndef BITIO_DOT_H


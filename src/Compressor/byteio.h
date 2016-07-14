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
#ifndef BYTEIO_DOT_H
#define BYTEIO_DOT_H

#include <type_traits>
#include <string>
#include "bitstream.h"
#include "definitions.h"

//
// Definition of a family of template classes used 
// for byte oriented output. Speicialized classes
// need to implement putByte()
//

template<typename T, typename Enable = void>
class output_bytes
{
public :
  //
  // If you try to instantiate an output_bytes<T>
  // object for a type that doesn't have a specialization,
  // you will get an error.
  //
  output_bytes(...) {
    static_assert( !std::is_void<Enable>::value, "Instantiating output_bytes<> without a specialization" );
  };
public :
  void putByte(char){}
};

template<typename T>
class output_bytes<T,typename std::enable_if<std::is_base_of<std::ostream, T>::value>::type>
{
public :
  output_bytes(T &stream) : m_stream(stream){}
  void putByte(char c)
  {
    m_stream.put(c);
  }
    /*
     * Added by Philipp Schelske
     * The following functions are added, so we can access basic functionality of the streamclass.
     */
    std::streampos tellp(){
        return m_stream.tellp();
    }
    void seekp(std::streampos pos, std::ios_base::seekdir way){
        m_stream.seekp(pos, way);
    }
    void seekp(std::streampos pos){
        m_stream.seekp(pos);
    }
    void writeUint64(const uint64_t x){
        m_stream.writeUint64(x);
    }
    void writeByte(const BYTE byte){
        m_stream.writeByte(byte);
    }
private :
  T &m_stream;
};

//
// Specialization of output_bytes for FILE *
//
/*
 * Uncommented by Philipp Schelske
 * Unlikely to be used in the Project.
 */
/*
template<>
class output_bytes<FILE *>
{
public :
  output_bytes(FILE *pFile)
  : m_pFile(pFile)
  {}
  void putByte(char c)
  {
    putc(c,m_pFile);
  }
private :
  FILE *m_pFile;
};
*/
//
// Definition of a family of template classes used 
// for byte oriented input. Speicialized classes
// need to implement getByte()
//

template<typename T,typename Enable = void>
class input_bytes
{
public :
  //
  // If you try to instantiate an input_bytes<T>
  // object for a type that doesn't have a specialization,
  // you will get an error indicating that you are 
  // trying to use this private constructor. 
  //
  input_bytes(...) {
    static_assert( !std::is_void<Enable>::value, "Instantiating input_bytes<> without a specialization" );
  };
public :
  int getByte();
    /*
     * Added by Philipp Schelske
     * This function will be used, so that the throw exception is not triggered anymore.
     *
     */
    bool eof();
};

//
// Specialization of input_bytes for class istream
//
template<typename T>
class input_bytes<T,typename std::enable_if<std::is_base_of<std::istream, T>::value>::type>
{
public :
  input_bytes(T &stream) 
  : m_stream(stream)
  {
  }
  int getByte()
  {
    return m_stream.get();
  }
    
    /*
     * Added by Philipp Schelske
     * The following functions are added, so we can access basic functionality of the streamclass.
     */

    bool eof(){
        return m_stream.eof();
    }
    void seekg(std::streampos pos, std::ios_base::seekdir way){
        m_stream.seekg(pos, way);
    }
    void seekg(std::streampos pos){
        m_stream.seekg(pos);
    }
    int readUint64(uint64_t &x){
        return m_stream.readUint64(x);
    }
    int readByte(BYTE &byte){
        return m_stream.readByte(byte);
    }
    
    std::streampos tellg(){
        return m_stream.tellg();
    }
    
private :
  T &m_stream;
};



//
// Specialization of input_bytes for FILE *
//
/*
 * Uncommented by Philipp Schelske
 * File does not support a simple EOF function. Curretnly I am to lazy to add this functionality.
 * Therefore I removed the support of this type. Unlikely to be used in the project we are working on.
 */
/*
template<>
class input_bytes<FILE *>
{
public :
  input_bytes(FILE *pFile) 
  : m_pFile( pFile )
  {}
  int getByte() {
      eof()
    return getc(m_pFile);
  }
    
private :
    FILE *m_pFile;
};*/


#endif //BYTEIO_DOT_H

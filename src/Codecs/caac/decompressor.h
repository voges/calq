/** @file decompressor.h
 *  @brief This file contains the definition of the decompressor template class.
 *  @author Philipp Schelske (schelske)
 *  @bug No known bugs
 */

#ifndef DECOMPRESSOR_H
#define DECOMPRESSOR_H

#include "bitstream.h"
#include <vector>


/* @brief class Decompressor
 *
 * The decompressor provides the functionality, to decompress a binary file, compressed with the compressor class.
 * Due to the fact, that neither the decompressor nor the compressor are adaptive, both need to use the exact same model.
 * 
 */

template<typename MODEL, typename INPUT>
class Decompressor
{
    typedef typename MODEL::CODE_VALUE CODE_VALUE;
    typedef typename MODEL::prob prob;
    public :
    
    Decompressor(INPUT &input,  MODEL &model )
    : m_input(input),
    m_model(model)
    {
    }
    
    void startBlock(){
        m_model.reset();
        m_input.reset();
        
        beginBlock = m_input.tellg();
        m_input.readUint64(uncompressedSize);
        m_input.readUint64(compressedSize);

        sizeCounter=0;
        high = MODEL::MAX_CODE;
        low = 0;
        value = 0;
    }
    
    int decode(std::vector<int> &output)
    {
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
            
            output.push_back(c);
            sizeCounter++;
            
            if(sizeCounter == uncompressedSize)
                break;
            
            high = low + (range*p.high)/p.count -1;
            low = low + (range*p.low)/p.count;
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
                
                if ((std::streampos)(beginBlock+(std::streampos)16+(std::streampos)compressedSize) <= m_input.tellg())  break;
                value += m_input.get_bit() ? 1 : 0;
            }
            
        }
        return 0;
    }
    
    std::streampos tellg(){
        return m_input.tellg();
    }
    int eof(){
        return m_input.eof();
    }
    
    private :
    INPUT &m_input;
    MODEL &m_model;
    
    CODE_VALUE high = MODEL::MAX_CODE;
    CODE_VALUE low = 0;
    CODE_VALUE value = 0;
    std::streampos beginBlock;
    uint64_t uncompressedSize;
    uint64_t compressedSize;
    uint64_t sizeCounter;
};


#endif /* COMPRESSOR:_H */

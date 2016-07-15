/** @file compressor.h
 *  @brief This file contains the definition of the artihmetic compressor
 *  @author Philipp Schelske (schelske)
 *  @bug No known bugs
 */



#ifndef COMPRESSOR_H
#define COMPRESSOR_H

#include <stdio.h>
#include <iostream>
#include <fstream>
#include "bitio.h"


#include <stdexcept>
#include "Compressors/bitstream.h"
#include "Compressors/caac/model_metrics.h"
#include "Compressors/caac/modelA.h"


/* @brief class Compressor
 * This class can be used to encode character per character with an arithmetic compressor.
 * It is important to use the same model, the decompressor will use.
 *
 */



template<typename MODEL, typename OUTPUT>
class Compressor {
    typedef typename MODEL::CODE_VALUE CODE_VALUE;
    typedef typename MODEL::prob prob;

public:
    Compressor<MODEL, OUTPUT>(OUTPUT &stream, MODEL &model);
    ~Compressor<MODEL, OUTPUT>();
    
    void encodeSymbol(int c);
    void finishBlock();
    void startBlock();

private:
    void put_bit_plus_pending(bool bit, int &pending_bits);
    OUTPUT &m_output;
    MODEL &m_model;
    int pending_bits = 0;
    CODE_VALUE low = 0;
    CODE_VALUE high = MODEL::MAX_CODE;
    
    void incrementCompressedSize();
    void incrementUncompressedSize();
    
    std::streampos blockBegin = 0x00;

    //compressed and uncompressed size are both in BYTE
    uint64_t compressedSize = 0;
    uint64_t uncompressedSize = 0;
    uint8_t compressedBitCounter = 0;
    uint8_t uncompressedBitCounter = 0;
    
};


/*
 * Added by Philipp Schelske
 * This function will start a new block in the outputstream. We reserve 128 Bit to save the uncompressed size
 * as well as the compressed size in byte.
 */
template<typename MODEL, typename OUTPUT>
void Compressor<MODEL, OUTPUT>::startBlock(){
    m_model.reset();
    blockBegin = m_output.tellp();
    std::streamoff offset = 2*sizeof(uint64_t);
    m_output.seekp(blockBegin+offset);

    pending_bits = 0;
    low = 0;
    high = MODEL::MAX_CODE;
    
    compressedSize = 0;
    uncompressedSize = 0;
    compressedBitCounter = 0;
    uncompressedBitCounter = 0;


}

template<typename MODEL, typename OUTPUT>
Compressor<MODEL, OUTPUT>::Compressor(OUTPUT &stream, MODEL &model)
:m_output(stream),
m_model(model)
{
    
}

template<typename MODEL, typename OUTPUT>
Compressor<MODEL, OUTPUT>::~Compressor()
{

}

template<typename MODEL, typename OUTPUT>
void Compressor<MODEL, OUTPUT>::encodeSymbol(int c){
    if ( c == -1 )
        c = 256;
    else
        uncompressedSize++;
    
    prob p = m_model.getProbability( c );
    CODE_VALUE range = high - low + 1;
    high = low + (range * p.high / p.count) - 1;
    low = low + (range * p.low / p.count);

    // On each pass there are six possible configurations of high/low,
    // each of which has its own    set of actions. When high or low
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
}

template<typename MODEL, typename OUTPUT>
inline void Compressor<MODEL, OUTPUT>::put_bit_plus_pending(bool bit, int &pending_bits){
    m_output.put_bit(bit);
    incrementCompressedSize();
    for ( int i = 0 ; i < pending_bits ; i++ )
    {
        m_output.put_bit(!bit);
        incrementCompressedSize();
    }
    pending_bits = 0;
    
}

/*
 * @author Philipp Schelske (schelske)
 * This function will start a new block in the outputstream. We reserve  128 Bit to save the uncompressed size
 * as well as the compressed size in byte.
 */
template<typename MODEL, typename OUTPUT>
void Compressor<MODEL, OUTPUT>::finishBlock(){
    pending_bits++;
    if ( low < MODEL::ONE_FOURTH )
        put_bit_plus_pending(0, pending_bits);
    else
        put_bit_plus_pending(1, pending_bits);

    compressedSize += m_output.finishBlock();
    
    std::streampos pos = m_output.tellp();
    m_output.seekp(blockBegin);
    m_output.writeUint64(uncompressedSize);
    m_output.writeUint64(compressedSize);
    m_output.seekp(pos);
}

template<typename MODEL, typename OUTPUT>
inline void Compressor<MODEL, OUTPUT>::incrementCompressedSize(){
    compressedBitCounter++;
    
    if(compressedBitCounter == 8){
        compressedBitCounter = 0;
        compressedSize++;
    }
}

#endif //COMPRESSOR_H
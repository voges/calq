/** @file QualDecoderWrapper.cc
 *  @brief C wrappers for the QualDecoder class
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#include "QualDecoderWrapper.h"
#include "QualDecoder.hpp"
#include <stdio.h>

extern "C" {

qualdecoder_t * qualdecoder_new(void) 
{
    cq::QualDecoder *qualDecoder = new cq::QualDecoder();
    return (qualdecoder_t *)qualDecoder;
}

void qualdecoder_delete(qualdecoder_t *qualdecoder) 
{
    cq::QualDecoder *qualDecoder = (cq::QualDecoder *)qualdecoder;
    delete qualDecoder;
    qualdecoder = NULL;
}

void qualdecoder_decode_block(qualdecoder_t * const qualdecoder, FILE *fp, str_t* qual[], size_t n)
{
    cq::QualDecoder *qualDecoder = (cq::QualDecoder *)qualdecoder;
    qualDecoder->decodeBlock(fp, qual, n);
}

}


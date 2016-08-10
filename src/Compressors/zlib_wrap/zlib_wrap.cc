/** @file zlib_wrap.c
 *  @brief This file contains the implementation of the zlib wrapper.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "zlib_wrap.h"
#include "Common/Exceptions.h"
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

unsigned char * zlib_compress(unsigned char *in,
                              const size_t  in_sz,
                              size_t        *out_sz)
{
    compressBound(in_sz);
    *out_sz = compressBound(in_sz) + 1;
    Byte *out = (Byte *)calloc((uInt)*out_sz, 1);
    int err = compress(out, out_sz, (const Bytef *)in, (uLong)in_sz);
    if (err != Z_OK) {
        throwErrorException("zlib failed to compress");
    }
    return out;
}

unsigned char * zlib_decompress(unsigned char *in,
                                const size_t  in_sz,
                                const size_t  out_sz)
{
    Bytef *out = (Bytef *)malloc(out_sz);
    if (out == NULL) {
        throwErrorException("malloc failed");
    }
    int err = uncompress(out, (uLongf *)&out_sz, (const Bytef *)in, (uLong)in_sz);
    if (err != Z_OK) {
        throwErrorException("zlib failed to uncompress");
    }
    return out;
}


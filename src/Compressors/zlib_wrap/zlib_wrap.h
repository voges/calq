/** @file zlib_wrap.h
 *  @brief This file contains the interface to the zlib wrapper.
 *
 *  Note it is up to the calling code to ensure that no overruns on input and
 *  output buffers occur!
 *
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef ZLIB_WRAP_H
#define ZLIB_WRAP_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>

unsigned char * zlib_compress(unsigned char *in,
                              const size_t  in_sz,
                              size_t        *out_sz);
unsigned char * zlib_decompress(unsigned char *in,
                                const size_t  in_sz,
                                const size_t  out_sz);

#ifdef __cplusplus
}
#endif

#endif // ZLIB_WRAP_H


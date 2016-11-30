/** @file rice.h
 *  @brief This file contains the interface to the Rice codec.
 *
 *  Note it is up to the calling code to ensure that no overruns on input and
 *  output buffers occur!
 *
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#ifndef RICE_H_
#define RICE_H_

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

unsigned char * rice_compress(unsigned char *in,
                              const size_t  in_sz,
                              size_t        *out_sz);
unsigned char * rice_decompress(unsigned char *in,
                                const size_t  in_sz,
                                size_t        *out_sz);

#ifdef __cplusplus
}
#endif

#endif // RICE_H_


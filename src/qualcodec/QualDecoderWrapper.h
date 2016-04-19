/** @file QualDecoderWrapper.h
 *  @brief C wrappers for the QualDecoder class
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#ifndef CQ_QUALDECODERWRAPPER_H
#define CQ_QUALDECODERWRAPPER_H

#include "misc/str.h"
#include "sam/samrec.h"
#include <stdio.h>

typedef void qualdecoder_t;

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Allocates a new qualdecoder instance
 *  @return Returns a pointer to the allocated object
 */
qualdecoder_t * qualdecoder_new(void);

/** @brief Deletes (and frees) a qualencoder object; after freeing, the pointer
 *         is set to NULL.
 */
void qualdecoder_delete(qualdecoder_t *qualdecoder);

/** @brief Decodes a block of encoded quality score from the given stream and
 *         writes the decoded quality scores to qual.
 *  @param qualdecoder The qualdecoder instance
 *  @param fp The file pointer to read from; has to be positioned at the
 *         beginning of an encoded quality score block
 *  @param qual An preallocated array of empty strings
 *  @param n Number of quality score vectors to decode
 */
void qualdecoder_decode_block(qualdecoder_t * const qualdecoder, FILE *fp, str_t* qual[], size_t n);

#ifdef __cplusplus
}
#endif

#endif // CQ_QUALDECODERWRAPPER_H


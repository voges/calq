/** @file QualEncoderWrapper.h
 *  @brief C wrappers for the QualEncoder class
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#ifndef QUALENCODERWRAPPER_H
#define QUALENCODERWRAPPER_H

#include "sam/samrec.h"
#include <stdio.h>

typedef void qualencoder_t;

#ifdef __cplusplus
extern "C" {
#endif

/** @brief Constructor: qualencoder_new
 *
 *  Allocates a new qualencoder instance.
 *
 *  @return Returns a pointer to the allocated object.
 */
qualencoder_t * qualencoder_new(void);

/** @brief Destructor: qualencoder_delete
 *
 *  Deletes (and frees) a qualencoder object; after freeing, the pointer
 *  is set to NULL.
 *
 *  @param qualencoder The qualencoder instance to delete.
 *  @return Void.
 */
void qualencoder_delete(qualencoder_t *qualencoder);

/** @brief Encodes the quality scores in a SAM record.
 * 
 *  Extracts the quality scores (plus some other information) from a SAM record
 *  and encodes them. The encoded quality score are kept in a buffer which can
 *  be written to a file.
 *
 *  @param qualencoder The qualencoder instance to operate with.
 *  @param samrec The SAM record to be encoded.
 *  @return Void.
 */
void qualencoder_encode_record(qualencoder_t * const qualencoder, const samrec_t * const samrec);

/** @brief Write out encoded quality scores and reset the qualencoder.
 *
 *  Writes encoded quality score to the given stream.
 *
 *  @param qualencoder The qualencoder instance
 *  @param fp Pointer to the output file
 *  @return Returns the number of written bytes
 */
size_t qualencoder_finish_block(qualencoder_t * const qualencoder, FILE *fp);

#ifdef __cplusplus
}
#endif

#endif // QUALENCODERWRAPPER_H


/** @file QualEncoderWrapper.cc
 *  @brief C wrappers for the QualEncoder class
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#include "QualEncoderWrapper.h"
#include "QualEncoder.hpp"
#include <stdio.h>

extern "C" {

qualencoder_t * qualencoder_new(void) 
{
    QualEncoder *qualEncoder = new QualEncoder();
    return (qualencoder_t *)qualEncoder;
}

void qualencoder_delete(qualencoder_t *qualencoder) 
{
    QualEncoder *qualEncoder = (QualEncoder *)qualencoder;
    delete qualEncoder;
    qualencoder = NULL;
}

void qualencoder_encode_record(qualencoder_t * const qualencoder, const samrec_t * const samrec)
{
    QualEncoder *qualEncoder = (QualEncoder *)qualencoder;
    qualEncoder->encodeRecord(samrec);
}

size_t qualencoder_finish_block(qualencoder_t * const qualencoder, FILE *fp)
{
    QualEncoder *qualEncoder = (QualEncoder *)qualencoder;
    return qualEncoder->finishBlock(fp);
}

}


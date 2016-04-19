#ifndef CQ_CQCODEC_H
#define CQ_CQCODEC_H

#include "qualcodec/QualDecoderWrapper.h"
#include "qualcodec/QualEncoderWrapper.h"
#include "sam/samparser.h"
#include <stdio.h>

typedef struct cqcodec_t_ {
    FILE          *ifp;
    FILE          *ofp;
    size_t        block_size;
    samparser_t   *samparser;
    qualencoder_t *qualencoder;
    qualencoder_t *qualdecoder;
} cqcodec_t;

cqcodec_t * cqcodec_new(FILE *ifp, FILE *ofp, const size_t block_size);
void cqcodec_delete(cqcodec_t *cqcodec);
int cqcodec_encode(cqcodec_t * const cqcodec);
int cqcodec_decode(cqcodec_t * const cqcodec);
int cqcodec_info(cqcodec_t * const cqcodec);

#endif // CQ_CQCODEC_H


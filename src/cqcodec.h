#ifndef CQ_CQCODEC_H
#define CQ_CQCODEC_H

#include "samparser.h"
#include "qualcodec.h"
#include <stdio.h>

typedef struct cqcodec_t_ {
    FILE        *ifp;
    FILE        *ofp;
    size_t      blk_sz;
    samparser_t *samparser;
    qualcodec_t *qualcodec;
} cqcodec_t;

cqcodec_t * cqcodec_new(FILE *ifp, FILE *ofp, size_t blocksz);
void cqcodec_free(cqcodec_t *cqcodec);
void cqcodec_encode(cqcodec_t *cqcodec);
void cqcodec_decode(cqcodec_t *cqcodec);
void cqcodec_info(cqcodec_t *cqcodec);

#endif // CQ_CQCODEC_H


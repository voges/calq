#ifndef CQ_SAMPARSER_H
#define CQ_SAMPARSER_H

#include "misc/str.h"
#include "sam/samrec.h"
#include <stdbool.h>
#include <stdio.h>

typedef struct samparser_t_ {
    FILE     *fp;   // file pointer
    str_t    *head; // SAM header
    samrec_t curr;  // current SAM record
} samparser_t;

samparser_t * samparser_new(FILE * const fp);
void samparser_delete(samparser_t *samparser);
int samparser_head(samparser_t * const samparser);
bool samparser_next(samparser_t * const samparser);

#endif // CQ_SAMPARSER_H


#ifndef CQ_QUALCODEC_H
#define CQ_QUALCODEC_H

#include "str.h"
#include <stdbool.h>
#include <stdio.h>

typedef struct qualcodec_t_ {
    size_t record_cnt; // number of records processed in the current block
    size_t pos_min;
    size_t pos_max;
    size_t *depths;
    size_t depths_len;
} qualcodec_t;

qualcodec_t * qualcodec_new(void);
void qualcodec_free(qualcodec_t *qualcodec);

// encoder methods 
bool qualcodec_add_record(qualcodec_t *qualcodec, const uint32_t pos, const char *cigar, const char *seq, const char *qual);
size_t qualcodec_finish_block(qualcodec_t *qualcodec, FILE *fp);

// decoder methods
size_t qualcodec_decode_block(qualcodec_t *qualcodec, FILE *fp, str_t **qual);

#endif // CQ_QUALCODEC_H


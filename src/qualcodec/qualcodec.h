#ifndef CQ_QUALCODEC_H
#define CQ_QUALCODEC_H

#include "str.h"
#include <stdio.h>

typedef struct qualcodec_t_ {
    size_t record_cnt; // number of records processed in the current block
} qualcodec_t;

qualcodec_t * qualcodec_new(void);
void qualcodec_free(qualcodec_t *qualcodec);

// Encoder methods
void qualcodec_add_record(qualcodec_t *qualcodec, uint32_t pos, const char *seq, const char *qual);
size_t qualcodec_write_block(qualcodec_t *qualcodec, FILE *fp);

// Decoder methods
size_t qualcodec_decode_block(qualcodec_t *qualcodec, FILE *fp, str_t **qual);

#endif // CQ_QUALCODEC_H


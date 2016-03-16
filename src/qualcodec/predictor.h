#ifndef CQ_PREDICTOR_H
#define CQ_PREDICTOR_H

#include <stdint.h>
#include <stdlib.h>

typedef struct predictor_t_ {
    size_t memory_sz;
    size_t alphabet_sz;
    uint32_t *table;
} predictor_t;

predictor_t * predictor_new(const size_t memory_sz, const size_t alphabet_sz);
void predictor_free(predictor_t *predictor);

void predictor_update(predictor_t *predictor, int *memory, int symbol);
int predictor_predict(predictor_t *predictor, int *memory);
void predictor_reset(predictor_t *predictor);

#endif // CQ_PREDICTOR_H
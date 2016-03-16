#include "predictor.h"
#include "cqlib.h"

static void init(predictor_t *predictor, const size_t memory_sz, const size_t alphabet_sz)
{
    predictor->memory_sz = memory_sz;
    predictor->alphabet_sz = alphabet_sz;
}

predictor_t * predictor_new(const size_t memory_sz, const size_t alphabet_sz)
{
    predictor_t *predictor = (predictor_t *)cq_malloc(sizeof(predictor_t));
    predictor->table = (uint32_t *)cq_calloc(memory_sz * alphabet_sz, sizeof(uint32_t));
    init(predictor, memory_sz, alphabet_sz);
    return predictor;
}

void predictor_free(predictor_t *predictor)
{
    if (predictor != NULL) {
        cq_free(predictor->table);
        free(predictor);
        predictor = NULL;
    } else {
        cq_error("Tried to free null pointer\n");
    }
}

void predictor_update(predictor_t *predictor, int *memory, int symbol)
{
    
}

int predictor_predict(predictor_t *predictor, int *memory)
{
    
}

void predictor_reset(predictor_t *predictor)
{
    
}


#include "qualcodec.h"
#include "common.h"
#include "cqlib.h"

static void qualcodec_init(qualcodec_t *qualcodec)
{
    qualcodec->record_cnt = 0;
}

qualcodec_t * qualcodec_new(void)
{
    qualcodec_t *qualcodec = (qualcodec_t *)cq_malloc(sizeof(qualcodec_t));
    qualcodec_init(qualcodec);
    return qualcodec;
}

void qualcodec_free(qualcodec_t *qualcodec)
{
    if (qualcodec != NULL) {
        free(qualcodec);
        qualcodec = NULL;
    } else {
        cq_error("Tried to free null pointer\n");
    }
}

void qualcodec_add_record(qualcodec_t *qualcodec, uint32_t pos, const char *seq, const char *qual)
{
    DEBUG("%d\n", pos);
    DEBUG("%s\n", seq);
    DEBUG("%s\n", qual);
}

size_t qualcodec_write_block(qualcodec_t *qualcodec, FILE * fp)
{
    size_t ret = 0;
    qualcodec_init(qualcodec);
    return ret;
}

size_t qualcodec_decode_block(qualcodec_t *qualcodec, FILE *fp, str_t **qual)
{
    size_t ret = 0;
    qualcodec_init(qualcodec);
    return ret;
}


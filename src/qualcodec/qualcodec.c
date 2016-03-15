#include "qualcodec.h"
#include "common.h"
#include "cqlib.h"
#include <ctype.h>
#include <string.h>

static void init(qualcodec_t *qualcodec)
{
    qualcodec->record_cnt = 0;
    qualcodec->pos_min = 1;
    qualcodec->pos_max = 1;
    qualcodec->depths_len = 1;
}

static void reset(qualcodec_t *qualcodec)
{
    qualcodec->record_cnt = 0;
    qualcodec->pos_min = 1;
    qualcodec->pos_max = 1;
    qualcodec->depths_len = 1;
    cq_free(qualcodec->depths);
    qualcodec->depths = (size_t *)cq_calloc(1, sizeof(size_t));
}

qualcodec_t * qualcodec_new(void)
{
    qualcodec_t *qualcodec = (qualcodec_t *)cq_malloc(sizeof(qualcodec_t));
    qualcodec->depths = (size_t *)cq_calloc(1, sizeof(size_t));
    init(qualcodec);
    return qualcodec;
}

void qualcodec_free(qualcodec_t *qualcodec)
{
    if (qualcodec != NULL) {
        cq_free(qualcodec->depths);
        free(qualcodec);
        qualcodec = NULL;
    } else {
        cq_error("Tried to free null pointer\n");
    }
}

static void expand(str_t *exp, const char *cigar, const char *seq)
{
    size_t cigar_idx = 0;
    size_t cigar_len = strlen(cigar);
    size_t op_len = 0; // length of current CIGAR operation
    size_t seq_idx = 0;

    for (cigar_idx = 0; cigar_idx < cigar_len; cigar_idx++) {
        if (isdigit(cigar[cigar_idx])) {
            op_len = op_len * 10 + (size_t)cigar[cigar_idx] - (size_t)'0';
            continue;
        }

        size_t i = 0;
        switch (cigar[cigar_idx]) {
        case 'M':
        case '=':
        case 'X':
            // add matching part to expanded sequence
            str_append_cstrn(exp, &seq[seq_idx], op_len);
            seq_idx += op_len;
            break;
        case 'I':
        case 'S':
            seq_idx += op_len; // skip inserted part
            break;
        case 'D':
        case 'N':
            // inflate expanded sequence
            for (i = 0; i < op_len; i++) str_append_char(exp, 'D');
            break;
        case 'H':
        case 'P':
            break; // these have been clipped
        default: cq_error("Bad CIGAR string: %s\n", cigar);
        }

        op_len = 0;
    }
}

bool qualcodec_add_record(qualcodec_t *qualcodec, const uint32_t pos, const char *cigar, const char *seq, const char *qual)
{
    qualcodec->record_cnt++;

    // check if this alignment is complete
    if (   (pos == 0)
        || (strlen(cigar) == 0 || (cigar[0] == '*' && cigar[1] == '\0'))
        || (strlen(seq) == 0 || (seq[0] == '*' && seq[1] == '\0'))
        || (strlen(qual) == 0 || (qual[0] == '*' && qual[1] == '\0'))) {
        cq_log("Alignment #%d is incomplete, skipping it\n", qualcodec->record_cnt);
        return false;
    }

    // expand current sequence
    str_t *exp = str_new();
    expand(exp, cigar, seq);

    // if this is the first record in a new block, simply allocate depths
    // vector; otherwise reallocate depths vector to cover the new region
    if (qualcodec->record_cnt == 1) {
        qualcodec->pos_min = pos;
        qualcodec->pos_max = pos + exp->len - 1;
        qualcodec->depths_len = qualcodec->pos_max - qualcodec->pos_min + 1;
        qualcodec->depths = (size_t *)cq_realloc(qualcodec->depths, sizeof(size_t) * qualcodec->depths_len);
        memset(qualcodec->depths, 0x00, sizeof((qualcodec->depths)) * qualcodec->depths_len);
    } else {
        size_t pos_max_new = pos + exp->len - 1;
        if (pos_max_new > qualcodec->pos_max) {
            size_t pos_max_old = qualcodec->pos_max;
            qualcodec->pos_max = pos_max_new;
            qualcodec->depths_len = qualcodec->pos_max-qualcodec->pos_min+1;
            qualcodec->depths = (size_t *)cq_realloc(qualcodec->depths, sizeof(size_t) * qualcodec->depths_len);
            memset(&(qualcodec->depths[pos_max_old - qualcodec->pos_min + 1]), 0x00, sizeof((qualcodec->depths)) * (pos_max_new-pos_max_old));
        }
    }

    // accumulate depths
    size_t i = 0;
    for (i = 0; i < exp->len; i++) {
        if (exp->s[i] == 'D') continue; // skip deletions
        qualcodec->depths[pos - qualcodec->pos_min + i]++;
    }

    // TODO: design/update quantizers for genomic columns
    // TODO: quantize the added quality score vector
    // TODO: perform Markov-chain prediction on current vector
    // TODO: pass vector to arithmetic coder

    str_free(exp);
    return true;
}

size_t qualcodec_finish_block(qualcodec_t *qualcodec, FILE *fp)
{
    size_t ret = 0;

    // DUMMY BEGIN
    ret += cq_fwrite_uint64(fp, qualcodec->record_cnt);
    size_t i = 0;
    for (i = 0; i < qualcodec->record_cnt; i++) {
        ret += cq_fwrite_byte(fp, 'q');
    }
    // DUMMY END

    reset(qualcodec);
    return ret;
}

size_t qualcodec_decode_block(qualcodec_t *qualcodec, FILE *fp, str_t **qual)
{
    size_t ret = 0;

    // DUMMY BEGIN
    unsigned char byte = 0x00;
    ret += cq_fread_uint64(fp, &(qualcodec->record_cnt));
    size_t i = 0;
    for (i = 0; i < qualcodec->record_cnt; i++) {
        ret += cq_fread_byte(fp, &byte);
        str_append_char(qual[i], (char)byte);
    }
    // DUMMY END

    reset(qualcodec);
    return ret;
}


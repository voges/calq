#include "sam/samparser.h"
#include "cqlib.h"
#include <string.h>

static void init(samparser_t * const samparser, FILE * const fp)
{
    samparser->fp = fp;
}

samparser_t * samparser_new(FILE * const fp)
{
    samparser_t *samparser = (samparser_t *)cq_malloc(sizeof(samparser_t));
    samparser->head = str_new();
    init(samparser, fp);
    return samparser;
}

void samparser_delete(samparser_t *samparser)
{
    if (samparser != NULL) {
        str_free(samparser->head);
        free(samparser);
        samparser = NULL;
    } else {
        cq_err("Tried to free null pointer\n");
		exit(EXIT_FAILURE);
    }
}

int samparser_head(samparser_t * const samparser)
{
    bool samheader = false;

    while (fgets(samparser->curr.line, sizeof(samparser->curr.line), samparser->fp)) {
        if (*(samparser->curr.line) == '@') {
            str_append_cstr(samparser->head, samparser->curr.line);
            samheader = true;
        } else {
            if (!samheader) {
                cq_err("SAM header missing\n");
                return CQ_FAILURE;
            }
            size_t offset = -strlen(samparser->curr.line);
            fseek(samparser->fp, (long)offset, SEEK_CUR);
            break;
        }
    }

    return CQ_SUCCESS;
}

static void parse(samparser_t * const samparser)
{
    size_t l = strlen(samparser->curr.line) - 1;

    while (l && (samparser->curr.line[l] == '\r' || samparser->curr.line[l] == '\n'))
        samparser->curr.line[l--] = '\0';

    char *c = samparser->curr.qname = samparser->curr.line;
    int f = 1;

    while (*c) {
        if (*c == '\t') {
            if (f ==  1) samparser->curr.flag  = (uint16_t)atoi(c + 1);
            if (f ==  2) samparser->curr.rname = c + 1;
            if (f ==  3) samparser->curr.pos   = (uint32_t)atoi(c + 1);
            if (f ==  4) samparser->curr.mapq  = (uint8_t)atoi(c + 1);
            if (f ==  5) samparser->curr.cigar = c + 1;
            if (f ==  6) samparser->curr.rnext = c + 1;
            if (f ==  7) samparser->curr.pnext = (uint32_t)atoi(c + 1);
            if (f ==  8) samparser->curr.tlen  = (int64_t)atoi(c + 1);
            if (f ==  9) samparser->curr.seq   = c + 1;
            if (f == 10) samparser->curr.qual  = c + 1;
            if (f == 11) samparser->curr.opt   = c + 1;
            f++;
            *c = '\0';
            if (f == 12) break;
        }
        c++;
    }

    if (f == 11) samparser->curr.opt = c;
}

bool samparser_next(samparser_t * const samparser)
{
    // try to read and parse next line
    if (fgets(samparser->curr.line, sizeof(samparser->curr.line), samparser->fp)) {
        parse(samparser);
    } else {
        return false;
    }

    return true;
}


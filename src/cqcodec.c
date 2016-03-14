#include "cqcodec.h"
#include "common.h"
#include "cqconfig.h"
#include "cqlib.h"
#include <string.h>
#include <sys/time.h>

static void cqcodec_init(cqcodec_t *cqcodec, FILE *ifp, FILE *ofp, size_t blk_sz)
{
    cqcodec->ifp = ifp;
    cqcodec->ofp = ofp;
    cqcodec->blk_sz = blk_sz;
}

cqcodec_t * cqcodec_new(FILE *ifp, FILE *ofp, size_t blk_sz)
{
    cqcodec_t *cqcodec = (cqcodec_t *)cq_malloc(sizeof(cqcodec_t));
    cqcodec->samparser = samparser_new(ifp);
    cqcodec->qualcodec = qualcodec_new();
    cqcodec_init(cqcodec, ifp, ofp, blk_sz);
    return cqcodec;
}

void cqcodec_free(cqcodec_t *cqcodec)
{
    if (cqcodec != NULL) {
        samparser_free(cqcodec->samparser);
        qualcodec_free(cqcodec->qualcodec);
        free(cqcodec);
        cqcodec = NULL;
    } else {
        cq_error("Tried to free null pointer\n");
    }
}

void cqcodec_encode(cqcodec_t *cqcodec)
{
    size_t qual_sz = 0; // accumulated QUAL size
    size_t cq_sz = 0; // compressed CALQ file size
    long fpos_prev = 0; // fp offset of the previous block
    size_t blk_n = 0; // block counter
    size_t rec_n = 0; // record counter
    size_t rec_cnt = 0; // record counter for current block
    size_t rec_max = cqcodec->blk_sz; // number of records per block

    struct timeval tv0, tv1;
    gettimeofday(&tv0, NULL);

    // write calq file header, skip 2x8 bytes for block and record counter
    unsigned char magic[5] = "calq";
    unsigned char version_major = CQ_VERSION_MAJOR + 48;
    unsigned char version_minor = CQ_VERSION_MINOR + 48;
    cq_sz += cq_fwrite_buf(cqcodec->ofp, magic, sizeof(magic));
    cq_sz += cq_fwrite_byte(cqcodec->ofp, version_major);
    cq_sz += cq_fwrite_byte(cqcodec->ofp, version_minor);
    fseek(cqcodec->ofp, 2*sizeof(uint64_t), SEEK_CUR);

    // parse (and seek past) SAM header
    samparser_head(cqcodec->samparser);

    samrec_t *samrec = &(cqcodec->samparser->curr);
    while (samparser_next(cqcodec->samparser)) {
        if (rec_cnt >= rec_max) {
            rec_cnt = 0;
            
            // store the file pointer offset of this block in the previous
            // block header
            long fpos_curr = ftell(cqcodec->ofp);
            if (blk_n > 0) {
                fseek(cqcodec->ofp, fpos_prev+4, SEEK_SET);
                cq_fwrite_uint64(cqcodec->ofp, (uint64_t)fpos_curr);
                fseek(cqcodec->ofp, fpos_curr, SEEK_SET);
            }
            fpos_prev = fpos_curr;

            // write last block
            cq_fwrite_uint64(cqcodec->ofp, (uint64_t)ftell(cqcodec->ofp));
            cq_fwrite_uint64(cqcodec->ofp, 0); // last block header has a zero here
            cq_sz += qualcodec_write_block(cqcodec->qualcodec, cqcodec->ofp);
            blk_n++;
        }

        // add new record
        qualcodec_add_record(cqcodec->qualcodec, samrec->pos, samrec->seq, samrec->qual);
        qual_sz += strlen(samrec->qual);
        rec_cnt++;
        rec_n++;
    }

    // store the file pointer offset of this block in the previous
    // block header
    long fpos_curr = ftell(cqcodec->ofp);
    if (blk_n > 0) {
        fseek(cqcodec->ofp, fpos_prev+4, SEEK_SET);
        cq_fwrite_uint64(cqcodec->ofp, (uint64_t)fpos_curr);
        fseek(cqcodec->ofp, fpos_curr, SEEK_SET);
    }
    fpos_prev = fpos_curr;

    // write last block
    cq_fwrite_uint64(cqcodec->ofp, (uint64_t)ftell(cqcodec->ofp));
    cq_fwrite_uint64(cqcodec->ofp, 0); // last block header has a zero here
    cq_sz += qualcodec_write_block(cqcodec->qualcodec, cqcodec->ofp);
    blk_n++;

    // finish file header
    // TODO

    // print summary
    gettimeofday(&tv1, NULL);
    cq_log("Took %ld us ~= %.2f s\n", tvdiff(tv0, tv1), (double)tvdiff(tv0, tv1)/1000000);
    cq_log("Compressed %zu record(s)\n", rec_n);
    cq_log("Wrote %zu block(s)\n", blk_n);
    cq_log("Uncompressed QUAL size: %zu\n", qual_sz);
    cq_log("Compressed CQ size: %zu\n", cq_sz);
    cq_log("CR: %.2f%%\n", qual_sz/cq_sz);
    cq_log("CF: %.2f%%\n", cq_sz/qual_sz);
}

void cqcodec_decode(cqcodec_t *cqcodec)
{
    // TODO
    cq_log("cqcodec_decode not yet implemented\n");
}


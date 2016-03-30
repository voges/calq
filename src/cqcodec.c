#include "cqcodec.h"
#include "common.h"
#include "cqconfig.h"
#include "cqlib.h"
#include <string.h>
#include <sys/time.h>

static void init(cqcodec_t *cqcodec, FILE *ifp, FILE *ofp, size_t blk_sz)
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
    init(cqcodec, ifp, ofp, blk_sz);
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
        cq_err("Tried to free null pointer\n");
        exit(EXIT_FAILURE);
    }
}

int cqcodec_encode(cqcodec_t *cqcodec)
{
    struct timeval tv0, tv1;
    gettimeofday(&tv0, NULL);

    size_t qual_sz = 0; // accumulated QUAL size
    size_t cq_sz = 0; // compressed CALQ file size
    long fpos_prev = 0; // fp offset of the previous block
    size_t blk_n = 0; // block counter
    size_t rec_n = 0; // record counter
    size_t rec_cnt = 0; // record counter for current block
    size_t rec_max = cqcodec->blk_sz; // number of records per block

    // write calq file header, skip 2x8 bytes for block and record counters
    unsigned char magic[5] = "calq";
    unsigned char version_major = CQ_VERSION_MAJOR + '0';
    unsigned char version_minor = CQ_VERSION_MINOR + '0';
    unsigned char version_patch = CQ_VERSION_PATCH + '0';
    cq_sz += cq_fwrite_buf(cqcodec->ofp, magic, sizeof(magic));
    cq_sz += cq_fwrite_byte(cqcodec->ofp, version_major);
    cq_sz += cq_fwrite_byte(cqcodec->ofp, version_minor);
    cq_sz += cq_fwrite_byte(cqcodec->ofp, version_patch);
    cq_sz += cq_fwrite_uint64(cqcodec->ofp, (uint64_t)rec_max);
    fseek(cqcodec->ofp, sizeof(uint64_t), SEEK_CUR); // space for blk_n
    fseek(cqcodec->ofp, sizeof(uint64_t), SEEK_CUR); // space for rec_n

    // parse (and seek past) SAM header
    if (CQ_SUCCESS != samparser_head(cqcodec->samparser)) {
        cq_err("Failed to parse SAM header!\n");
        return CQ_FAILURE;
    }

    samrec_t *samrec = &(cqcodec->samparser->curr);
    while (samparser_next(cqcodec->samparser)) {
        if (rec_cnt >= rec_max) {
            rec_cnt = 0;

            // store the file pointer offset of this block in the previous
            // block header
            long fpos_curr = ftell(cqcodec->ofp);
            if (blk_n > 0) {
                fseek(cqcodec->ofp, fpos_prev+(long)sizeof(uint64_t), SEEK_SET);
                cq_sz += cq_fwrite_uint64(cqcodec->ofp, (uint64_t)fpos_curr);
                fseek(cqcodec->ofp, fpos_curr, SEEK_SET);
            }
            fpos_prev = fpos_curr;

            // write block: fpos, fpos_next, data
            cq_sz += cq_fwrite_uint64(cqcodec->ofp, (uint64_t)ftell(cqcodec->ofp));
            fseek(cqcodec->ofp, sizeof(uint64_t), SEEK_CUR); // space for fpos_next
            cq_sz += cq_fwrite_uint64(cqcodec->ofp, (uint64_t)rec_max);
            cq_sz += qualcodec_finish_block(cqcodec->qualcodec, cqcodec->ofp);
            blk_n++;
        }

        // add new record
        qualcodec_add_record(cqcodec->qualcodec, samrec->pos, samrec->cigar, samrec->seq, samrec->qual);
        qual_sz += strlen(samrec->qual);
        rec_cnt++;
        rec_n++;
    }

    // store the file pointer offset of this block in the previous
    // block header
    long fpos_curr = ftell(cqcodec->ofp);
    if (blk_n > 0) {
        fseek(cqcodec->ofp, fpos_prev+(long)sizeof(uint64_t), SEEK_SET);
        cq_sz += cq_fwrite_uint64(cqcodec->ofp, (uint64_t)fpos_curr);
        fseek(cqcodec->ofp, fpos_curr, SEEK_SET);
    }
    fpos_prev = fpos_curr;

    // write last block: fpos, fpos_next=0, rec_cnt, data
    cq_sz += cq_fwrite_uint64(cqcodec->ofp, (uint64_t)ftell(cqcodec->ofp));
    cq_sz += cq_fwrite_uint64(cqcodec->ofp, 0); // last block header has a zero here
    cq_sz += cq_fwrite_uint64(cqcodec->ofp, (uint64_t)rec_cnt);
    cq_sz += qualcodec_finish_block(cqcodec->qualcodec, cqcodec->ofp);
    blk_n++;

    // finish file header
    size_t off = sizeof(magic) + sizeof(version_major) + sizeof(version_minor) + sizeof(version_patch) + sizeof(rec_max);
    fseek(cqcodec->ofp, (long)off, SEEK_SET);
    cq_sz += cq_fwrite_uint64(cqcodec->ofp, (uint64_t)blk_n);
    cq_sz += cq_fwrite_uint64(cqcodec->ofp, (uint64_t)rec_n);
    fseek(cqcodec->ofp, 0, SEEK_END);

    // print summary
    gettimeofday(&tv1, NULL);
    cq_out("Took %ld us ~= %.2f s\n", tvdiff(tv0, tv1), (double)tvdiff(tv0, tv1)/1000000);
    cq_out("Compressed %zu record(s)\n", rec_n);
    cq_out("Wrote %zu block(s)\n", blk_n);
    cq_out("Uncompressed QUAL size: %zu\n", qual_sz);
    cq_out("Compressed CQ size: %zu\n", cq_sz);
    cq_out("Compression Ratio (CR): %.2f%%\n", (double)qual_sz/(double)cq_sz*100);
    cq_out("Compression Factor (CF): %.2f%%\n", (double)cq_sz/(double)qual_sz*100);

    return CQ_SUCCESS;
}

int cqcodec_decode(cqcodec_t *cqcodec)
{
    struct timeval tv0, tv1;
    gettimeofday(&tv0, NULL);

    // read and check file header
    unsigned char magic[5];
    unsigned char version_major = 0;
    unsigned char version_minor = 0;
    unsigned char version_patch = 0;
    uint64_t rec_max = 0;
    uint64_t blk_n = 0;
    uint64_t rec_n = 0;

    cq_fread_buf(cqcodec->ifp, magic, sizeof(magic));
    cq_fread_byte(cqcodec->ifp, &version_major);
    cq_fread_byte(cqcodec->ifp, &version_minor);
    cq_fread_byte(cqcodec->ifp, &version_patch);
    cq_fread_uint64(cqcodec->ifp, &rec_max);
    cq_fread_uint64(cqcodec->ifp, &blk_n);
    cq_fread_uint64(cqcodec->ifp, &rec_n);

    if (version_major-'0' != CQ_VERSION_MAJOR || version_minor-'0' != CQ_VERSION_MINOR) {
        cq_err("Program version: %d.%d.%d\n", CQ_VERSION_MAJOR, CQ_VERSION_MINOR, CQ_VERSION_PATCH);
        cq_err("File version: %c.%c.%c\n", version_major, version_minor, version_patch);
        cq_err("Program version does not match file version\n");
        return CQ_FAILURE;
    }

    size_t b = 0, r = 0;
    for (b = 0; b < blk_n; b++) {
        // read number of records in current block
        uint64_t rec_cnt = 0;
        fseek(cqcodec->ifp, 2*sizeof(uint64_t), SEEK_CUR);
        cq_fread_uint64(cqcodec->ifp, &rec_cnt);

        // allocate memory for decoded quality scores
        str_t **qual =(str_t **)cq_malloc(sizeof(str_t *) * rec_max);
        for (r = 0; r < rec_max; r++) qual[r] = str_new();

        // decode block
        qualcodec_decode_block(cqcodec->qualcodec, cqcodec->ifp, qual);

        // write decoded quality score to output
        for (r = 0; r < rec_cnt; r++) {
            cq_fwrite_buf(cqcodec->ofp, (unsigned char *)qual[r]->s, qual[r]->len);
            cq_fwrite_byte(cqcodec->ofp, '\n');
            str_free(qual[r]);
        }

        free(qual);
    }

    // print summary
    gettimeofday(&tv1, NULL);
    cq_out("Took %ld us ~= %.2f s\n", tvdiff(tv0, tv1), (double)tvdiff(tv0, tv1)/1000000);
    cq_out("Decoded %zu record(s) in %zu block(s)\n", rec_n, blk_n);

    return CQ_SUCCESS;
}

void cqcodec_info(cqcodec_t *cqcodec)
{
    // read and check file header
    unsigned char magic[5];
    unsigned char version_major = 0;
    unsigned char version_minor = 0;
    unsigned char version_patch = 0;
    uint64_t rec_max = 0;
    uint64_t blk_n = 0;
    uint64_t rec_n = 0;

    cq_fread_buf(cqcodec->ifp, magic, sizeof(magic));
    cq_fread_byte(cqcodec->ifp, &version_major);
    cq_fread_byte(cqcodec->ifp, &version_minor);
    cq_fread_byte(cqcodec->ifp, &version_patch);
    cq_fread_uint64(cqcodec->ifp, &rec_max);
    cq_fread_uint64(cqcodec->ifp, &blk_n);
    cq_fread_uint64(cqcodec->ifp, &rec_n);

    cq_out("magic: %s\n", magic);
    cq_out("version: %c.%c.%c\n", version_major, version_minor, version_patch);
    cq_out("block size: %"PRIu64"\n", rec_max);
    cq_out("blocks: %"PRIu64"\n", blk_n);
    cq_out("records: %"PRIu64"\n", rec_n);

    // read and print block headers
    printf("\n        fpos     fpos_next       rec_cnt\n");

    while (1) {
        uint64_t fpos = 0;
        uint64_t fpos_next = 0;
        uint64_t rec_cnt = 0;
        cq_fread_uint64(cqcodec->ifp, &fpos);
        cq_fread_uint64(cqcodec->ifp, &fpos_next);
        cq_fread_uint64(cqcodec->ifp, &rec_cnt);

        printf("%12"PRIu64"  ", fpos);
        printf("%12"PRIu64"  ", fpos_next);
        printf("%12"PRIu64"  ", rec_cnt);
        printf("\n");

        if (fpos_next)
            fseek(cqcodec->ifp, (long)fpos_next, SEEK_SET);
        else
            break; // last block has zeros here
    }
    printf("\n");
}


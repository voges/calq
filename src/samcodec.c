/*
 * The copyright in this software is being made available under the TNT
 * License, included below. This software may be subject to other third party
 * and contributor rights, including patent rights, and no such rights are
 * granted under this license.
 *
 * Copyright (c) 2015, Leibniz Universitaet Hannover, Institut fuer
 * Informationsverarbeitung (TNT)
 * Contact: <voges@tnt.uni-hannover.de>
 * All rights reserved.
 *
 * * Redistribution in source or binary form is not permitted.
 *
 * * Use in source or binary form is only permitted in the context of scientific
 *   research.
 *
 * * Commercial use without specific prior written permission is prohibited.
 *   Neither the name of the TNT nor the names of its contributors may be used
 *   to endorse or promote products derived from this software without specific
 *   prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "samcodec.h"
#include "common/common.h"
#include "osro.h"
#include "tscfmt.h"
#include "tsclib.h"
#include "version.h"
#include <inttypes.h>
#include <string.h>
#include <sys/time.h>

// Indices for SAM statistics
enum {
    SAM_QNAME,
    SAM_FLAG,
    SAM_RNAME,
    SAM_POS,
    SAM_MAPQ,
    SAM_CIGAR,
    SAM_RNEXT,
    SAM_PNEXT,
    SAM_TLEN,
    SAM_SEQ,
    SAM_QUAL,
    SAM_OPT,
    SAM_CTRL, // Count \t's and \n's here
    SAM_HEAD  // SAM header
};

// Indices for tsc statistics
enum {
    TSC_FH,     // File header
    TSC_SH,     // SAM header
    TSC_BH,     // Block header(s)
    TSC_AUX,    // Total no. of bytes written by auxcodec
    TSC_ID,     // Total no. of bytes written by idcodec
    TSC_NUC,    // Total no. of bytes written by nuccodec
    TSC_PAIR,   // Total no. of bytes written by paircodec
    TSC_QUAL    // Total no. of bytes written by qualcodec
};

// Indices for timing statistics
enum {
    ET_TOT, // Total time elapsed
    ET_AUX,
    ET_ID,
    ET_NUC,
    ET_PAIR,
    ET_QUAL,
    ET_REM  // Remaining (I/O, statistics, etc.)
};

static void samcodec_init(samcodec_t   *samcodec,
                          FILE         *ifp,
                          FILE         *ofp,
                          unsigned int blk_sz)
{
    samcodec->ifp = ifp;
    samcodec->ofp = ofp;
    samcodec->blk_sz = blk_sz;
}

samcodec_t * samcodec_new(FILE *ifp, FILE *ofp, unsigned int blk_sz)
{
    samcodec_t *samcodec = (samcodec_t *)osro_malloc(sizeof(samcodec_t));
    samcodec->samparser = samparser_new(ifp);
    samcodec->auxcodec = auxcodec_new();
    samcodec->idcodec = idcodec_new();
    samcodec->nuccodec = nuccodec_new();
    samcodec->paircodec = paircodec_new();
    samcodec->qualcodec = qualcodec_new();
    samcodec_init(samcodec, ifp, ofp, blk_sz);
    return samcodec;
}

void samcodec_free(samcodec_t *samcodec)
{
    if (samcodec != NULL) {
        samparser_free(samcodec->samparser);
        auxcodec_free(samcodec->auxcodec);
        idcodec_free(samcodec->idcodec);
        nuccodec_free(samcodec->nuccodec);
        paircodec_free(samcodec->paircodec);
        qualcodec_free(samcodec->qualcodec);
        free(samcodec);
        samcodec = NULL;
    } else {
        tsc_error("Tried to free null pointer\n");
    }
}

static void samcodec_print_stats(const size_t  *sam_sz,
                                 const size_t  *tsc_sz,
                                 const tscfh_t *tscfh,
                                 const long    *et)
{
    size_t sam_total_sz = sam_sz[SAM_QNAME]
                        + sam_sz[SAM_FLAG]
                        + sam_sz[SAM_RNAME]
                        + sam_sz[SAM_POS]
                        + sam_sz[SAM_MAPQ]
                        + sam_sz[SAM_CIGAR]
                        + sam_sz[SAM_RNEXT]
                        + sam_sz[SAM_PNEXT]
                        + sam_sz[SAM_TLEN]
                        + sam_sz[SAM_SEQ]
                        + sam_sz[SAM_QUAL]
                        + sam_sz[SAM_OPT]
                        + sam_sz[SAM_CTRL]  // \t's and \n's
                        + sam_sz[SAM_HEAD]; // Don't forget the SAM header
    size_t sam_aux_sz = sam_sz[SAM_FLAG]
                      + sam_sz[SAM_MAPQ]
                      + sam_sz[SAM_OPT];
    size_t sam_id_sz = sam_sz[SAM_QNAME];
    size_t sam_nuc_sz = sam_sz[SAM_RNAME]
                      + sam_sz[SAM_POS]
                      + sam_sz[SAM_CIGAR]
                      + sam_sz[SAM_SEQ];
    size_t sam_pair_sz = sam_sz[SAM_RNEXT]
                         + sam_sz[SAM_PNEXT]
                         + sam_sz[SAM_TLEN];
    size_t sam_qual_sz = sam_sz[SAM_QUAL];
    size_t tsc_total_sz = tsc_sz[TSC_FH]
                        + tsc_sz[TSC_SH]
                        + tsc_sz[TSC_BH]
                        + tsc_sz[TSC_AUX]
                        + tsc_sz[TSC_ID]
                        + tsc_sz[TSC_NUC]
                        + tsc_sz[TSC_PAIR]
                        + tsc_sz[TSC_QUAL];

    fprintf(stdout,
            "\n"
            "\tStatistics:\n"
            "\t-----------\n"
            "\tNumber of records   : %12"PRIu64"\n"
            "\tNumber of blocks    : %12"PRIu64"\n"
            "\n"
            "\tSAM file size       : %12zu (%6.2f%%)\n"
            "\t  QNAME             : %12zu (%6.2f%%)\n"
            "\t  FLAG              : %12zu (%6.2f%%)\n"
            "\t  RNAME             : %12zu (%6.2f%%)\n"
            "\t  POS               : %12zu (%6.2f%%)\n"
            "\t  MAPQ              : %12zu (%6.2f%%)\n"
            "\t  CIGAR             : %12zu (%6.2f%%)\n"
            "\t  RNEXT             : %12zu (%6.2f%%)\n"
            "\t  PNEXT             : %12zu (%6.2f%%)\n"
            "\t  TLEN              : %12zu (%6.2f%%)\n"
            "\t  SEQ               : %12zu (%6.2f%%)\n"
            "\t  QUAL              : %12zu (%6.2f%%)\n"
            "\t  OPT               : %12zu (%6.2f%%)\n"
            "\t  CTRL (\\t, \\n)     : %12zu (%6.2f%%)\n"
            "\t  HEAD (SAM header) : %12zu (%6.2f%%)\n"
            "\n"
            "\tTsc file size       : %12zu (%6.2f%%)\n"
            "\t  File header       : %12zu (%6.2f%%)\n"
            "\t  SAM header        : %12zu (%6.2f%%)\n"
            "\t  Block header(s)   : %12zu (%6.2f%%)\n"
            "\t  Aux               : %12zu (%6.2f%%)\n"
            "\t  Id                : %12zu (%6.2f%%)\n"
            "\t  Nuc               : %12zu (%6.2f%%)\n"
            "\t  Pair              : %12zu (%6.2f%%)\n"
            "\t  Qual              : %12zu (%6.2f%%)\n"
            "\n"
            "\tCompression ratios             SAM /          tsc\n"
            "\t  Total             : %12zu / %12zu (%6.2f%%)\n"
            "\t  Aux               : %12zu / %12zu (%6.2f%%)\n"
            "\t  Id                : %12zu / %12zu (%6.2f%%)\n"
            "\t  Nuc               : %12zu / %12zu (%6.2f%%)\n"
            "\t  Pair              : %12zu / %12zu (%6.2f%%)\n"
            "\t  Qual              : %12zu / %12zu (%6.2f%%)\n"
            "\n"
            "\tTiming\n"
            "\t  Total time elapsed: %12ld us ~= %12.2f s (%6.2f%%)\n"
            "\t  Aux               : %12ld us ~= %12.2f s (%6.2f%%)\n"
            "\t  Id                : %12ld us ~= %12.2f s (%6.2f%%)\n"
            "\t  Nuc               : %12ld us ~= %12.2f s (%6.2f%%)\n"
            "\t  Pair              : %12ld us ~= %12.2f s (%6.2f%%)\n"
            "\t  Qual              : %12ld us ~= %12.2f s (%6.2f%%)\n"
            "\t  Remaining         : %12ld us ~= %12.2f s (%6.2f%%)\n"
            "\n"
            "\tSpeed\n"
            "\t  Total             : %12.2f MB/s\n"
            "\t  Aux               : %12.2f MB/s\n"
            "\t  Id                : %12.2f MB/s\n"
            "\t  Nuc               : %12.2f MB/s\n"
            "\t  Pair              : %12.2f MB/s\n"
            "\t  Qual              : %12.2f MB/s\n"
            "\n",
            tscfh->rec_n,
            tscfh->blk_n,

            sam_total_sz,
            (100 * (double)sam_total_sz / (double)sam_total_sz),
            sam_sz[SAM_QNAME],
            (100 * (double)sam_sz[SAM_QNAME] / (double)sam_total_sz),
            sam_sz[SAM_FLAG],
            (100 * (double)sam_sz[SAM_FLAG] / (double)sam_total_sz),
            sam_sz[SAM_RNAME],
            (100 * (double)sam_sz[SAM_RNAME] / (double)sam_total_sz),
            sam_sz[SAM_POS],
            (100 * (double)sam_sz[SAM_POS] / (double)sam_total_sz),
            sam_sz[SAM_MAPQ],
            (100 * (double)sam_sz[SAM_MAPQ] / (double)sam_total_sz),
            sam_sz[SAM_CIGAR],
            (100 * (double)sam_sz[SAM_CIGAR] / (double)sam_total_sz),
            sam_sz[SAM_RNEXT],
            (100 * (double)sam_sz[SAM_RNEXT] / (double)sam_total_sz),
            sam_sz[SAM_PNEXT],
            (100 * (double)sam_sz[SAM_PNEXT] / (double)sam_total_sz),
            sam_sz[SAM_TLEN],
            (100 * (double)sam_sz[SAM_TLEN] / (double)sam_total_sz),
            sam_sz[SAM_SEQ],
            (100 * (double)sam_sz[SAM_SEQ] / (double)sam_total_sz),
            sam_sz[SAM_QUAL],
            (100 * (double)sam_sz[SAM_QUAL] / (double)sam_total_sz),
            sam_sz[SAM_OPT],
            (100 * (double)sam_sz[SAM_OPT] / (double)sam_total_sz),
            sam_sz[SAM_CTRL],
            (100 * (double)sam_sz[SAM_CTRL] / (double)sam_total_sz),
            sam_sz[SAM_HEAD],
            (100 * (double)sam_sz[SAM_HEAD] / (double)sam_total_sz),

            tsc_total_sz,
            (100 * (double)tsc_total_sz / (double)tsc_total_sz),
            tsc_sz[TSC_FH],
            (100 * (double)tsc_sz[TSC_FH] / (double)tsc_total_sz),
            tsc_sz[TSC_SH],
            (100 * (double)tsc_sz[TSC_SH] / (double)tsc_total_sz),
            tsc_sz[TSC_BH],
            (100 * (double)tsc_sz[TSC_BH] / (double)tsc_total_sz),
            tsc_sz[TSC_AUX],
            (100 * (double)tsc_sz[TSC_AUX] / (double)tsc_total_sz),
            tsc_sz[TSC_ID],
            (100 * (double)tsc_sz[TSC_ID] / (double)tsc_total_sz),
            tsc_sz[TSC_NUC],
            (100 * (double)tsc_sz[TSC_NUC] / (double)tsc_total_sz),
            tsc_sz[TSC_PAIR],
            (100 * (double)tsc_sz[TSC_PAIR] / (double)tsc_total_sz),
            tsc_sz[TSC_QUAL],
            (100 * (double)tsc_sz[TSC_QUAL] / (double)tsc_total_sz),

            sam_total_sz,
            tsc_total_sz,
            (100*(double)tsc_total_sz/(double)sam_total_sz),
            sam_aux_sz,
            tsc_sz[TSC_AUX],
            (100*(double)tsc_sz[TSC_AUX]/(double)sam_aux_sz),
            sam_id_sz,
            tsc_sz[TSC_ID],
            (100*(double)tsc_sz[TSC_ID]/(double)sam_id_sz),
            sam_nuc_sz,
            tsc_sz[TSC_NUC],
            (100*(double)tsc_sz[TSC_NUC]/(double)sam_nuc_sz),
            sam_pair_sz,
            tsc_sz[TSC_PAIR],
            (100*(double)tsc_sz[TSC_PAIR]/(double)sam_pair_sz),
            sam_qual_sz,
            tsc_sz[TSC_QUAL],
            (100*(double)tsc_sz[TSC_QUAL]/(double)sam_qual_sz),

            et[ET_TOT],
            (double)et[ET_TOT] / (double)1000000,
            (100*(double)et[ET_TOT]/(double)et[ET_TOT]),
            et[ET_AUX],
            (double)et[ET_AUX] / (double)1000000,
            (100*(double)et[ET_AUX]/(double)et[ET_TOT]),
            et[ET_ID],
            (double)et[ET_ID] / (double)1000000,
            (100*(double)et[ET_ID]/(double)et[ET_TOT]),
            et[ET_NUC],
            (double)et[ET_NUC] / (double)1000000,
            (100*(double)et[ET_NUC]/(double)et[ET_TOT]),
            et[ET_PAIR],
            (double)et[ET_PAIR] / (double)1000000,
            (100*(double)et[ET_PAIR]/(double)et[ET_TOT]),
            et[ET_QUAL],
            (double)et[ET_QUAL] / (double)1000000,
            (100*(double)et[ET_QUAL]/(double)et[ET_TOT]),
            et[ET_REM],
            (double)et[ET_REM] / (double)1000000,
            (100*(double)et[ET_REM] / (double)et[ET_TOT]),

            ((double)sam_total_sz / MB) / ((double)et[ET_TOT] / (double)1000000),
            ((double)sam_aux_sz / MB) / ((double)et[ET_AUX] / (double)1000000),
            ((double)sam_id_sz / MB) / ((double)et[ET_ID] / (double)1000000),
            ((double)sam_nuc_sz / MB) / ((double)et[ET_NUC] / (double)1000000),
            ((double)sam_pair_sz / MB) / ((double)et[ET_PAIR] / (double)1000000),
            ((double)sam_qual_sz / MB) / ((double)et[ET_QUAL] / (double)1000000));
}

void samcodec_encode(samcodec_t *samcodec)
{
    size_t   sam_sz[14] = {0}; // SAM statistics
    size_t   tsc_sz[8]  = {0}; // Tsc statistics
    long     et[7]      = {0}; // Timing statistics
    long     fpos_prev  = 0;   // fp offset of the previous block
    struct timeval tt0, tt1, tv0, tv1;
    gettimeofday(&tt0, NULL);

    tscfh_t *tscfh = tscfh_new();
    tscsh_t *tscsh = tscsh_new();
    tscbh_t *tscbh = tscbh_new();

    // Set up tsc file header, then seek past it for now
    tscfh->flags = 0x1; // LSB signals SAM
    tscfh->sblk_n = 3; // aux, nux, qual
    fseek(samcodec->ofp, (long)tscfh_size(tscfh), SEEK_SET);

    // Copy SAM header to tsc file
    samparser_head(samcodec->samparser);
    sam_sz[SAM_HEAD] += tscsh->data_sz = (uint64_t)samcodec->samparser->head->len;
    tscsh->data = (unsigned char *)samcodec->samparser->head->s;
    tsc_sz[TSC_SH] = tscsh_write(tscsh, samcodec->ofp);
    tscsh->data = NULL; // Need this before freeing tscsh

    // Set up block header
    tscbh->rec_max = samcodec->blk_sz;

    samrec_t *samrec = &(samcodec->samparser->curr);
    while (samparser_next(samcodec->samparser)) {
        if (tscbh->rec_cnt >= tscbh->rec_max) {
            // Store the file pointer offset of this block in the previous
            // block header
            long fpos_curr = ftell(samcodec->ofp);
            if (tscbh->blk_cnt > 0) {
                fseek(samcodec->ofp, fpos_prev, SEEK_SET);
                fseek(samcodec->ofp, (long)sizeof(tscbh->fpos), SEEK_CUR);
                osro_fwrite_uint64(samcodec->ofp, (uint64_t)fpos_curr);
                fseek(samcodec->ofp, fpos_curr, SEEK_SET);
            }
            fpos_prev = fpos_curr;

            // Write block header
            tscbh->fpos = (uint64_t)ftell(samcodec->ofp);
            tsc_sz[TSC_BH] += tscbh_write(tscbh, samcodec->ofp);
            tscbh->blk_cnt++;
            tscbh->rec_cnt = 0;

            // Write sub-blocks
            gettimeofday(&tv0, NULL);
            tsc_sz[TSC_AUX] += auxcodec_write_block(samcodec->auxcodec, samcodec->ofp);
            gettimeofday(&tv1, NULL);
            et[ET_AUX] += tvdiff(tv0, tv1);

            gettimeofday(&tv0, NULL);
            tsc_sz[TSC_ID] += idcodec_write_block(samcodec->idcodec, samcodec->ofp);
            gettimeofday(&tv1, NULL);
            et[ET_ID] += tvdiff(tv0, tv1);

            gettimeofday(&tv0, NULL);
            tsc_sz[TSC_NUC] += nuccodec_write_block(samcodec->nuccodec, samcodec->ofp);
            gettimeofday(&tv1, NULL);
            et[ET_NUC] += tvdiff(tv0, tv1);

            gettimeofday(&tv0, NULL);
            tsc_sz[TSC_PAIR] += paircodec_write_block(samcodec->paircodec, samcodec->ofp);
            gettimeofday(&tv1, NULL);
            et[ET_PAIR] += tvdiff(tv0, tv1);

            gettimeofday(&tv0, NULL);
            tsc_sz[TSC_QUAL] += qualcodec_write_block(samcodec->qualcodec, samcodec->ofp);
            gettimeofday(&tv1, NULL);
            et[ET_QUAL] += tvdiff(tv0, tv1);
        }

        // Add record to encoders
        gettimeofday(&tv0, NULL);
        auxcodec_add_record(samcodec->auxcodec,
                            samrec->flag,
                            samrec->mapq,
                            samrec->opt);
        gettimeofday(&tv1, NULL);
        et[ET_AUX] += tvdiff(tv0, tv1);

        gettimeofday(&tv0, NULL);
        idcodec_add_record(samcodec->idcodec,
                         samrec->qname);
        gettimeofday(&tv1, NULL);
        et[ET_ID] += tvdiff(tv0, tv1);

        gettimeofday(&tv0, NULL);
        nuccodec_add_record(samcodec->nuccodec,
                            samrec->rname,
                            samrec->pos,
                            samrec->cigar,
                            samrec->seq);
        gettimeofday(&tv1, NULL);
        et[ET_NUC] += tvdiff(tv0, tv1);

        gettimeofday(&tv0, NULL);
        paircodec_add_record(samcodec->paircodec,
                           samrec->rnext,
                           samrec->pnext,
                           samrec->tlen);
        gettimeofday(&tv1, NULL);
        et[ET_QUAL] += tvdiff(tv0, tv1);

        gettimeofday(&tv0, NULL);
        qualcodec_add_record(samcodec->qualcodec,
                             samrec->qual);
        gettimeofday(&tv1, NULL);
        et[ET_QUAL] += tvdiff(tv0, tv1);

        // Accumulate input statistics
        sam_sz[SAM_QNAME] += strlen(samrec->qname);
        sam_sz[SAM_FLAG]  += ndigits(samrec->flag);
        sam_sz[SAM_RNAME] += strlen(samrec->rname);
        sam_sz[SAM_POS]   += ndigits(samrec->pos);
        sam_sz[SAM_MAPQ]  += ndigits(samrec->mapq);
        sam_sz[SAM_CIGAR] += strlen(samrec->cigar);
        sam_sz[SAM_RNEXT] += strlen(samrec->rnext);
        sam_sz[SAM_PNEXT] += ndigits(samrec->pnext);
        sam_sz[SAM_TLEN]  += ndigits(samrec->tlen);
        sam_sz[SAM_SEQ]   += strlen(samrec->seq);
        sam_sz[SAM_QUAL]  += strlen(samrec->qual);
        sam_sz[SAM_OPT]   += strlen(samrec->opt);
        sam_sz[SAM_CTRL]  += 11 /* \t's */ + 1 /* \n's */;

        tscbh->rec_cnt++;
        tscfh->rec_n++;
    }

    // Store the file pointer offset of this block in the previous
    // block header
    long fpos_curr = ftell(samcodec->ofp);
    if (tscbh->blk_cnt > 0) {
        fseek(samcodec->ofp, fpos_prev, SEEK_SET);
        fseek(samcodec->ofp, (long)sizeof(tscbh->fpos), SEEK_CUR);
        osro_fwrite_uint64(samcodec->ofp, (uint64_t)fpos_curr);
        fseek(samcodec->ofp, fpos_curr, SEEK_SET);
    }
    fpos_prev = fpos_curr;

    // Write -last- block header
    tscbh->fpos = (uint64_t)ftell(samcodec->ofp);
    tscbh->fpos_nxt = 0; // Last block header has a zero here
    tsc_sz[TSC_BH] += tscbh_write(tscbh, samcodec->ofp);
    tscbh->blk_cnt++;
    tscbh->rec_cnt = 0;

    // Write -last- sub-blocks
    gettimeofday(&tv0, NULL);
    tsc_sz[TSC_AUX] += auxcodec_write_block(samcodec->auxcodec, samcodec->ofp);
    gettimeofday(&tv1, NULL);
    et[ET_AUX] += tvdiff(tv0, tv1);

    gettimeofday(&tv0, NULL);
    tsc_sz[TSC_ID] += idcodec_write_block(samcodec->idcodec, samcodec->ofp);
    gettimeofday(&tv1, NULL);
    et[ET_ID] += tvdiff(tv0, tv1);

    gettimeofday(&tv0, NULL);
    tsc_sz[TSC_NUC] += nuccodec_write_block(samcodec->nuccodec, samcodec->ofp);
    gettimeofday(&tv1, NULL);
    et[ET_NUC] += tvdiff(tv0, tv1);

    gettimeofday(&tv0, NULL);
    tsc_sz[TSC_PAIR] += paircodec_write_block(samcodec->paircodec, samcodec->ofp);
    gettimeofday(&tv1, NULL);
    et[ET_PAIR] += tvdiff(tv0, tv1);

    gettimeofday(&tv0, NULL);
    tsc_sz[TSC_QUAL] += qualcodec_write_block(samcodec->qualcodec, samcodec->ofp);
    gettimeofday(&tv1, NULL);
    et[ET_QUAL] += tvdiff(tv0, tv1);

    // Write file header
    tscfh->blk_n = tscbh->blk_cnt;
    rewind(samcodec->ofp);
    tsc_sz[TSC_FH] = tscfh_write(tscfh, samcodec->ofp);
    fseek(samcodec->ofp, (long)0, SEEK_END);

    // Print summary
    gettimeofday(&tt1, NULL);
    et[ET_TOT] = tvdiff(tt0, tt1);
    et[ET_REM] = et[ET_TOT] - et[ET_AUX] - et[ET_ID]
               - et[ET_NUC] - et[ET_PAIR] - et[ET_QUAL];
    tsc_log("Compressed %zu record(s)\n", tscfh->rec_n);
    tsc_log("Wrote %zu block(s)\n", tscfh->blk_n);
    tsc_log("Took %ld us ~= %.2f s\n", et[ET_TOT], (double)et[ET_TOT]/1000000);

    // If selected, print detailed statistics
    if (tsc_stats) samcodec_print_stats(sam_sz, tsc_sz, tscfh, et);

    tscfh_free(tscfh);
    tscsh_free(tscsh);
    tscbh_free(tscbh);
}

void samcodec_decode(samcodec_t *samcodec)
{
    size_t   sam_sz[14] = {0}; // SAM statistics
    size_t   tsc_sz[8]  = {0}; // Tsc statistics
    long     et[7]      = {0}; // Timing statistics
    struct timeval tt0, tt1, tv0, tv1;
    gettimeofday(&tt0, NULL);

    tscfh_t *tscfh = tscfh_new();
    tscsh_t *tscsh = tscsh_new();
    tscbh_t *tscbh = tscbh_new();

    // Read and check file header
    tsc_sz[TSC_FH] = tscfh_read(tscfh, samcodec->ifp);

    // Read SAM header from tsc file and write it to SAM file
    tsc_sz[TSC_SH] = tscsh_read(tscsh, samcodec->ifp);
    sam_sz[SAM_HEAD] += osro_fwrite_buf(samcodec->ofp, tscsh->data, tscsh->data_sz);

    size_t b = 0;
    for (b = 0; b < tscfh->blk_n; b++) {
        tsc_sz[TSC_BH] += tscbh_read(tscbh, samcodec->ifp);

        // Allocate memory to prepare decoding of the sub-blocks
        str_t   **qname=(str_t **)osro_malloc(sizeof(str_t *)*tscbh->rec_cnt);
        uint16_t *flag =(uint16_t *)osro_malloc(sizeof(uint16_t)*tscbh->rec_cnt);
        str_t   **rname=(str_t **)osro_malloc(sizeof(str_t *)*tscbh->rec_cnt);
        uint32_t *pos  =(uint32_t *)osro_malloc(sizeof(uint32_t)*tscbh->rec_cnt);
        uint8_t *mapq  =(uint8_t *)osro_malloc(sizeof(uint8_t)*tscbh->rec_cnt);
        str_t   **cigar=(str_t **)osro_malloc(sizeof(str_t *)*tscbh->rec_cnt);
        str_t   **rnext=(str_t **)osro_malloc(sizeof(str_t *)*tscbh->rec_cnt);
        uint32_t *pnext=(uint32_t *)osro_malloc(sizeof(uint32_t)*tscbh->rec_cnt);
        int64_t *tlen  =(int64_t *)osro_malloc(sizeof(int64_t)*tscbh->rec_cnt);
        str_t   **seq  =(str_t **)osro_malloc(sizeof(str_t *)*tscbh->rec_cnt);
        str_t   **qual =(str_t **)osro_malloc(sizeof(str_t *)*tscbh->rec_cnt);
        str_t   **opt  =(str_t **)osro_malloc(sizeof(str_t *)*tscbh->rec_cnt);

        size_t r = 0;
        for (r = 0; r < tscbh->rec_cnt; r++) {
            qname[r] = str_new();
            flag[r]  = 0;
            rname[r] = str_new();
            pos[r]   = 0;
            mapq[r]  = 0;
            cigar[r] = str_new();
            rnext[r] = str_new();
            pnext[r] = 0;
            tlen[r]  = 0;
            seq[r]   = str_new();
            qual[r]  = str_new();
            opt[r]   = str_new();
        }

        // Decode sub-blocks
        gettimeofday(&tv0, NULL);
        tsc_sz[TSC_AUX] += auxcodec_decode_block(samcodec->auxcodec,
                                                 samcodec->ifp,
                                                 flag,
                                                 mapq,
                                                 opt);
        gettimeofday(&tv1, NULL);
        et[ET_AUX] += tvdiff(tv0, tv1);

        gettimeofday(&tv0, NULL);
        tsc_sz[TSC_ID] += idcodec_decode_block(samcodec->idcodec,
                                               samcodec->ifp,
                                               qname);
        gettimeofday(&tv1, NULL);
        et[ET_ID] += tvdiff(tv0, tv1);

        gettimeofday(&tv0, NULL);
        tsc_sz[TSC_NUC] += nuccodec_decode_block(samcodec->nuccodec,
                                                 samcodec->ifp,
                                                 rname,
                                                 pos,
                                                 cigar,
                                                 seq);
        gettimeofday(&tv1, NULL);
        et[ET_NUC] += tvdiff(tv0, tv1);

        gettimeofday(&tv0, NULL);
        tsc_sz[TSC_PAIR] += paircodec_decode_block(samcodec->paircodec,
                                                   samcodec->ifp,
                                                   rnext,
                                                   pnext,
                                                   tlen);
        gettimeofday(&tv1, NULL);
        et[ET_PAIR] += tvdiff(tv0, tv1);

        gettimeofday(&tv0, NULL);
        tsc_sz[TSC_QUAL] += qualcodec_decode_block(samcodec->qualcodec,
                                                   samcodec->ifp,
                                                   qual);
        gettimeofday(&tv1, NULL);
        et[ET_QUAL] += tvdiff(tv0, tv1);

        // These are dummies for testing
        //int i = 0;
        //for (i = 0; i < tscbh->rec_cnt; i++) {
            //str_append_cstr(qname[i], "qname");
            //flag[i] = 2146;
            //str_append_cstr(rname[i], "rname");
            //pos[i] = 905;
            //mapq[i] = 3490;
            //str_append_cstr(cigar[i], "cigar");
            //str_append_cstr(rnext[i], "rnext");
            //pnext[i] = 68307;
            //tlen[i] = 7138;
            //str_append_cstr(seq[i], "seq");
            //str_append_cstr(qual[i], "qual");
            //str_append_cstr(opt[i], "opt");
        //}

        // Write decoded sub-blocks in correct order to outfile
        for (r = 0; r < tscbh->rec_cnt; r++) {
            // Convert int-fields to C-strings (101 bytes should be enough)
            char flag_cstr[101];
            char pos_cstr[101];
            char mapq_cstr[101];
            char pnext_cstr[101];
            char tlen_cstr[101];

            snprintf(flag_cstr, sizeof(flag_cstr), "%"PRIu16, flag[r]);
            snprintf(pos_cstr, sizeof(pos_cstr), "%"PRIu32, pos[r]);
            snprintf(mapq_cstr, sizeof(mapq_cstr), "%"PRIu8, mapq[r]);
            snprintf(pnext_cstr, sizeof(pnext_cstr), "%"PRIu32, pnext[r]);
            snprintf(tlen_cstr, sizeof(tlen_cstr), "%"PRId64, tlen[r]);

            char *sam_fields[12];
            sam_fields[SAM_QNAME] = qname[r]->s;
            sam_fields[SAM_FLAG]  = flag_cstr;
            sam_fields[SAM_RNAME] = rname[r]->s;
            sam_fields[SAM_POS]   = pos_cstr;
            sam_fields[SAM_MAPQ]  = mapq_cstr;
            sam_fields[SAM_CIGAR] = cigar[r]->s;
            sam_fields[SAM_RNEXT] = rnext[r]->s;
            sam_fields[SAM_PNEXT] = pnext_cstr;
            sam_fields[SAM_TLEN]  = tlen_cstr;
            sam_fields[SAM_SEQ]   = seq[r]->s;
            sam_fields[SAM_QUAL]  = qual[r]->s;
            sam_fields[SAM_OPT]   = opt[r]->s;

            // Write data to file
            size_t f = 0;
            for (f = 0; f < 12; f++) {
                sam_sz[f] += osro_fwrite_buf(samcodec->ofp, (unsigned char *)sam_fields[f], strlen(sam_fields[f]));
                if (f != 11 && strlen(sam_fields[f + 1]))
                    sam_sz[SAM_CTRL] += osro_fwrite_byte(samcodec->ofp, '\t');
            }
            sam_sz[SAM_CTRL] += osro_fwrite_byte(samcodec->ofp, '\n');

            // Free the memory used for the current line
            str_free(qname[r]);
            str_free(rname[r]);
            str_free(cigar[r]);
            str_free(rnext[r]);
            str_free(seq[r]);
            str_free(qual[r]);
            str_free(opt[r]);
        }

        // Free memory allocated for this block
        free(qname);
        free(flag);
        free(rname);
        free(pos);
        free(mapq);
        free(cigar);
        free(rnext);
        free(pnext);
        free(tlen);
        free(seq);
        free(qual);
        free(opt);
    }

    // Print summary
    gettimeofday(&tt1, NULL);
    et[ET_TOT] = tvdiff(tt0, tt1);
    et[ET_REM] = et[ET_TOT] - et[ET_AUX] - et[ET_ID]
               - et[ET_NUC] - et[ET_PAIR] - et[ET_QUAL];
    tsc_log("Decompressed %zu record(s)\n", tscfh->rec_n);
    tsc_log("Read %zu block(s)\n", tscfh->blk_n);
    tsc_log("Took %ld us ~= %.2f s\n", et[ET_TOT], (double)et[ET_TOT]/1000000);

    // If selected, print detailed statistics
    if (tsc_stats) samcodec_print_stats(sam_sz, tsc_sz, tscfh, et);

    tscfh_free(tscfh);
    tscsh_free(tscsh);
    tscbh_free(tscbh);
}

void samcodec_info(samcodec_t *samcodec)
{
    tscfh_t *tscfh = tscfh_new();
    tscsh_t *tscsh = tscsh_new();
    tscbh_t *tscbh = tscbh_new();

    // Read and check file header
    tscfh_read(tscfh, samcodec->ifp);

    // Skip SAM header
    tscsh_read(tscsh, samcodec->ifp);

    // Read and print block headers
    printf("\n"
           "\t        fpos      fpos_nxt       blk_cnt       rec_cnt"
           "       rec_max         rname       pos_min       pos_max\n");

    while (1) {
        tscbh_read(tscbh, samcodec->ifp);

        printf("\t");
        printf("%12"PRIu64"  ", tscbh->fpos);
        printf("%12"PRIu64"  ", tscbh->fpos_nxt);
        printf("%12"PRIu64"  ", tscbh->blk_cnt);
        printf("%12"PRIu64"  ", tscbh->rec_cnt);
        printf("%12"PRIu64"  ", tscbh->rec_max);
        printf("%12"PRIu64"  ", tscbh->rname);
        printf("%12"PRIu64"  ", tscbh->pos_min);
        printf("%12"PRIu64"  ", tscbh->pos_max);
        printf("\n");

        if (tscbh->fpos_nxt)
            fseek(samcodec->ifp, (long)tscbh->fpos_nxt, SEEK_SET);
        else
            break; // Last block has zeros here
    }
    printf("\n");

    tscfh_free(tscfh);
    tscsh_free(tscsh);
    tscbh_free(tscbh);
}


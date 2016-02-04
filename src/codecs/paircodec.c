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

//
// Pair block format:
//   unsigned char id[8]          : "pair---" + '\0'
//   uint64_t      record_cnt     : No. of lines in block
//   uint64_t      uncompressed_sz: Size of uncompressed data
//   uint64_t      compressed_sz  : Compressed data size
//   uint64_t      compressed_crc : CRC64 of compressed data
//   unsigned char *compressed    : Compressed data
//

#include "paircodec.h"
#include "tsclib.h"
#include "osro.h"
#include "zlib_wrap.h"
#include <inttypes.h>
#include <string.h>

static void paircodec_init(paircodec_t *paircodec)
{
    paircodec->record_cnt = 0;
    str_clear(paircodec->uncompressed);
    if (paircodec->compressed != NULL) {
        free(paircodec->compressed);
        paircodec->compressed = NULL;
    }
    paircodec->compressed_sz = 0;
}

paircodec_t * paircodec_new(void)
{
    paircodec_t *paircodec = (paircodec_t *)osro_malloc(sizeof(paircodec_t));
    paircodec->uncompressed = str_new();
    paircodec->compressed = NULL;
    paircodec_init(paircodec);
    return paircodec;
}

void paircodec_free(paircodec_t *paircodec)
{
    if (paircodec != NULL) {
        str_free(paircodec->uncompressed);
        if (paircodec->compressed != NULL) {
            free(paircodec->compressed);
            paircodec->compressed = NULL;
        }
        free(paircodec);
        paircodec = NULL;
    } else {
        tsc_error("Tried to free null pointer\n");
    }
}

// Encoder methods
// -----------------------------------------------------------------------------

void paircodec_add_record(paircodec_t    *paircodec,
                          const char     *rnext,
                          const uint32_t pnext,
                          const int64_t  tlen)
{
    paircodec->record_cnt++;

    char pnext_cstr[101];
    char tlen_cstr[101];

    snprintf(pnext_cstr, sizeof(pnext_cstr), "%"PRIu32, pnext);
    snprintf(tlen_cstr, sizeof(tlen_cstr), "%"PRId64, tlen);

    str_append_cstr(paircodec->uncompressed, rnext);
    str_append_cstr(paircodec->uncompressed, "\t");
    str_append_cstr(paircodec->uncompressed, pnext_cstr);
    str_append_cstr(paircodec->uncompressed, "\t");
    str_append_cstr(paircodec->uncompressed, tlen_cstr);
    str_append_cstr(paircodec->uncompressed, "\n");
}

size_t paircodec_write_block(paircodec_t *paircodec, FILE *fp)
{
    size_t ret = 0;

    // Compress block
    unsigned char *uncompressed = (unsigned char *)paircodec->uncompressed->s;
    size_t uncompressed_sz = paircodec->uncompressed->len;
    size_t compressed_sz = 0;
    unsigned char *compressed = zlib_compress(uncompressed, uncompressed_sz, &compressed_sz);

    // Compute CRC64
    uint64_t compressed_crc = osro_crc64(compressed, compressed_sz);

    // Write compressed block
    unsigned char id[8] = "pair---"; id[7] = '\0';
    ret += osro_fwrite_buf(fp, id, sizeof(id));
    ret += osro_fwrite_uint64(fp, (uint64_t)paircodec->record_cnt);
    ret += osro_fwrite_uint64(fp, (uint64_t)uncompressed_sz);
    ret += osro_fwrite_uint64(fp, (uint64_t)compressed_sz);
    ret += osro_fwrite_uint64(fp, (uint64_t)compressed_crc);
    ret += osro_fwrite_buf(fp, compressed, compressed_sz);

    // Free memory allocated by zlib_compress
    free(compressed);

    paircodec_init(paircodec);

    return ret;
}

// Decoder methods
// -----------------------------------------------------------------------------

static size_t paircodec_decode(unsigned char *tmp,
                               size_t        tmp_sz,
                               str_t         **rnext,
                               uint32_t      *pnext,
                               int64_t       *tlen)
{
    size_t ret = 0;
    size_t i = 0;
    size_t line = 0;
    unsigned char *cstr = tmp;
    unsigned int idx = 0;

    for (i = 0; i < tmp_sz; i++) {
        if (tmp[i] == '\n') {
            tmp[i] = '\0';
            tlen[line] = strtoll((const char *)cstr, NULL, 10);
            ret += sizeof(*tlen);
            cstr = &tmp[i+1];
            idx = 0;
            line++;
            continue;
        }

        if (tmp[i] == '\t') {
            tmp[i] = '\0';
            switch (idx++) {
            case 0:
                str_append_cstr(rnext[line], (const char *)cstr);
                ret += strlen((const char *)cstr);
                break;
            case 1:
                pnext[line] = (uint32_t)strtoul((const char *)cstr, NULL, 10);
                ret += sizeof(*pnext);
                break;
            default:
                tsc_error("Wrong pair bit stream\n");
            }
            cstr = &tmp[i+1];
        }
    }

    return ret;
}

size_t paircodec_decode_block(paircodec_t *paircodec,
                              FILE        *fp,
                              str_t       **rnext,
                              uint32_t    *pnext,
                              int64_t     *tlen)
{
    size_t ret = 0;

    unsigned char id[8];
    uint64_t      record_cnt;
    uint64_t      uncompressed_sz;
    uint64_t      compressed_sz;
    uint64_t      compressed_crc;
    unsigned char *compressed;

    // Read block
    ret += osro_fread_buf(fp, id, sizeof(id));
    ret += osro_fread_uint64(fp, &record_cnt);
    ret += osro_fread_uint64(fp, &uncompressed_sz);
    ret += osro_fread_uint64(fp, &compressed_sz);
    ret += osro_fread_uint64(fp, &compressed_crc);
    compressed = (unsigned char *)osro_malloc((size_t)compressed_sz);
    ret += osro_fread_buf(fp, compressed, compressed_sz);

    // Check CRC64
    if (osro_crc64(compressed, compressed_sz) != compressed_crc)
        tsc_error("CRC64 check failed for pair block\n");

    // Decompress block
    unsigned char *uncompressed = zlib_decompress(compressed, compressed_sz, uncompressed_sz);
    free(compressed);

    // Decode block
    paircodec_decode(uncompressed, uncompressed_sz, rnext, pnext, tlen);
    free(uncompressed); // Free memory used for decoded bitstream

    paircodec_init(paircodec);

    return ret;
}


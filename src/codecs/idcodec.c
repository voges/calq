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
// Id block format:
//   unsigned char id[8]          : "id-----" + '\0'
//   uint64_t      record_cnt     : No. of lines in block
//   uint64_t      uncompressed_sz: Size of uncompressed data
//   uint64_t      compressed_sz  : Compressed data size
//   uint64_t      compressed_crc : CRC64 of compressed data
//   unsigned char *compressed    : Compressed data
//

#include "idcodec.h"
#include "tsclib.h"
#include "osro.h"
#include "zlib_wrap.h"
#include <string.h>

static void idcodec_init(idcodec_t *idcodec)
{
    idcodec->record_cnt = 0;
    str_clear(idcodec->uncompressed);
    if (idcodec->compressed != NULL) {
        free(idcodec->compressed);
        idcodec->compressed = NULL;
    }
    idcodec->compressed_sz = 0;
}

idcodec_t * idcodec_new(void)
{
    idcodec_t *idcodec = (idcodec_t *)osro_malloc(sizeof(idcodec_t));
    idcodec->uncompressed = str_new();
    idcodec->compressed = NULL;
    idcodec_init(idcodec);
    return idcodec;
}

void idcodec_free(idcodec_t *idcodec)
{
    if (idcodec != NULL) {
        str_free(idcodec->uncompressed);
        if (idcodec->compressed != NULL) {
            free(idcodec->compressed);
            idcodec->compressed = NULL;
        }
        free(idcodec);
        idcodec = NULL;
    } else {
        tsc_error("Tried to free null pointer\n");
    }
}

// Encoder methods
// -----------------------------------------------------------------------------

void idcodec_add_record(idcodec_t *idcodec, const char *qname)
{
    idcodec->record_cnt++;

    str_append_cstr(idcodec->uncompressed, qname);
    str_append_cstr(idcodec->uncompressed, "\n");
}

size_t idcodec_write_block(idcodec_t *idcodec, FILE * fp)
{
    size_t ret = 0;

    // Compress block
    unsigned char *uncompressed = (unsigned char *)idcodec->uncompressed->s;
    size_t uncompressed_sz = idcodec->uncompressed->len;
    size_t compressed_sz = 0;
    unsigned char *compressed = zlib_compress(uncompressed, uncompressed_sz, &compressed_sz);

    // Compute CRC64
    uint64_t compressed_crc = osro_crc64(compressed, compressed_sz);

    // Write compressed block
    unsigned char id[8] = "id-----"; id[7] = '\0';
    ret += osro_fwrite_buf(fp, id, sizeof(id));
    ret += osro_fwrite_uint64(fp, (uint64_t)idcodec->record_cnt);
    ret += osro_fwrite_uint64(fp, (uint64_t)uncompressed_sz);
    ret += osro_fwrite_uint64(fp, (uint64_t)compressed_sz);
    ret += osro_fwrite_uint64(fp, (uint64_t)compressed_crc);
    ret += osro_fwrite_buf(fp, compressed, compressed_sz);

    // Free memory allocated by zlib_compress
    free(compressed);

    idcodec_init(idcodec);

    return ret;
}

// Decoder methods
// -----------------------------------------------------------------------------

static size_t idcodec_decode(unsigned char *tmp,
                             size_t        tmp_sz,
                             str_t**       qname)
{
    size_t ret = 0;
    size_t i = 0;
    size_t rec = 0;

    for (i = 0; i < tmp_sz; i++) {
        if (tmp[i] != '\n') {
            str_append_char(qname[rec], (const char)tmp[i]);
            ret++;
        } else {
            rec++;
        }
    }

    return ret;
}

size_t idcodec_decode_block(idcodec_t *idcodec, FILE *fp, str_t **qname)
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
        tsc_error("CRC64 check failed for id block\n");

    // Decompress block
    unsigned char *uncompressed = zlib_decompress(compressed, compressed_sz, uncompressed_sz);
    free(compressed);

    // Decode block
    idcodec_decode(uncompressed, uncompressed_sz, qname);
    free(uncompressed); // Free memory used for decoded bitstream

    idcodec_init(idcodec);

    return ret;
}


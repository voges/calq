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

#ifndef TSC_NUCCODEC_O1_H
#define TSC_NUCCODEC_O1_H

#include "common/cbufint64.h"
#include "common/cbufstr.h"
#include "common/str.h"
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

#define TSC_NUCCODEC_O1_WINDOW_SZ 10
#define TSC_NUCCODEC_O1_ALPHA 0.9

typedef struct nuccodec_t_ {
    // Only used in encoder
    size_t record_cnt;  // No. of records processed in the current block
    bool   first;       // 'false', if first line has not been processed yet
    str_t  *rname_prev; // Holding current RNAME
    str_t  *ctrl;
    str_t  *poff;
    str_t  *stogy;
    str_t  *mod;
    str_t  *trail;

    // Circular buffers
    cbufint64_t *neo_cbuf;
    cbufint64_t *pos_cbuf;
    cbufstr_t   *exs_cbuf;
} nuccodec_t;

nuccodec_t * nuccodec_new(void);
void nuccodec_free(nuccodec_t *nuccodec);

// Encoder
// -----------------------------------------------------------------------------

void nuccodec_add_record(nuccodec_t     *nuccodec,
                         const char     *rname,
                         const uint32_t pos,
                         const char     *cigar,
                         const char     *seq);
size_t nuccodec_write_block(nuccodec_t *nuccodec, FILE *fp);

// Decoder methods
// -----------------------------------------------------------------------------

size_t nuccodec_decode_block(nuccodec_t *nuccodec,
                             FILE       *fp,
                             str_t      **rname,
                             uint32_t   *pos,
                             str_t      **cigar,
                             str_t      **seq);

#endif // TSC_NUCCODEC_O1_H


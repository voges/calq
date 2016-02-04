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

#include "cbufint64.h"
#include <stdio.h>
#include <string.h>

static void cbufint64_init(cbufint64_t *cbufint64, const size_t sz)
{
    cbufint64->sz = sz;
    cbufint64->nxt = 0;
    cbufint64->n = 0;
}

cbufint64_t * cbufint64_new(const size_t sz)
{
    cbufint64_t *cbufint64 = (cbufint64_t *)malloc(sizeof(cbufint64_t));
    if (!cbufint64) abort();
    cbufint64->buf = (int64_t *)malloc(sizeof(int64_t) * sz);
    if (!cbufint64->buf) abort();
    memset(cbufint64->buf, 0x00, sz * sizeof(int64_t));
    cbufint64_init(cbufint64, sz);
    return cbufint64;
}

void cbufint64_free(cbufint64_t *cbufint64)
{
    if (cbufint64 != NULL) {
        free(cbufint64->buf);
        free(cbufint64);
        cbufint64 = NULL;
    } else {
        fprintf(stderr, "Error: Tried to free null pointer\n");
        exit(EXIT_FAILURE);
    }
}

void cbufint64_clear(cbufint64_t *cbufint64)
{
    memset(cbufint64->buf, 0x00, cbufint64->sz * sizeof(int64_t));
    cbufint64->nxt = 0;
    cbufint64->n = 0;
}

void cbufint64_push(cbufint64_t *cbufint64, int64_t x)
{
    cbufint64->buf[cbufint64->nxt++] = x;
    if (cbufint64->nxt == cbufint64->sz)
        cbufint64->nxt = 0;
    if (cbufint64->n < cbufint64->sz)
        cbufint64->n++;
}

int64_t cbufint64_top(cbufint64_t *cbufint64)
{
    if (cbufint64->n == 0) {
        fprintf(stderr, "Error: Tried to access empty cbufint64\n");
        exit(EXIT_FAILURE);
    }

    size_t nxt = cbufint64->nxt;
    size_t last = 0;

    if (nxt == 0)
        last = cbufint64->sz - 1;
    else
        last = cbufint64->nxt - 1;

    return cbufint64->buf[last];
}

int64_t cbufint64_get(const cbufint64_t *cbufint64, size_t pos)
{
    if (cbufint64->n == 0) {
        fprintf(stderr, "Error: Tried to access empty cbufint64\n");
        exit(EXIT_FAILURE);
    }
    if (pos > (cbufint64->n - 1)) {
        fprintf(stderr, "Error: Not enough elements in cbufint64\n");
        exit(EXIT_FAILURE);
    }

    return cbufint64->buf[pos];
}


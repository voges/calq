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

#include "cbufstr.h"
#include <stdio.h>

static void cbufstr_init(cbufstr_t *cbufstr, const size_t sz)
{
    cbufstr->sz = sz;
    cbufstr->nxt = 0;
    cbufstr->n = 0;
}

cbufstr_t * cbufstr_new(const size_t sz)
{
    cbufstr_t *cbufstr = (cbufstr_t *)malloc(sizeof(cbufstr_t));
    if (!cbufstr) abort();
    cbufstr->buf = (str_t **)malloc(sizeof(str_t*) * sz);
    if (!cbufstr->buf) abort();
    size_t i = 0;
    for (i = 0; i < sz; i++)
        cbufstr->buf[i] = str_new();
    cbufstr_init(cbufstr, sz);
    return cbufstr;
}

void cbufstr_free(cbufstr_t *cbufstr)
{
    if (cbufstr != NULL) {
        size_t i = 0;
        for (i = 0; i < cbufstr->sz; i++)
            str_free(cbufstr->buf[i]);
        free(cbufstr->buf);
        free(cbufstr);
        cbufstr = NULL;
    } else {
        fprintf(stderr, "Error: Tried to free null pointer\n");
        exit(EXIT_FAILURE);
    }
}

void cbufstr_clear(cbufstr_t *cbufstr)
{
    size_t i = 0;
    for (i = 0; i < cbufstr->sz; i++)
        str_clear(cbufstr->buf[i]);
    cbufstr->nxt = 0;
    cbufstr->n = 0;
}

void cbufstr_push(cbufstr_t *cbufstr, const char *s)
{
    str_copy_cstr(cbufstr->buf[cbufstr->nxt++], s);
    if (cbufstr->nxt == cbufstr->sz)
        cbufstr->nxt = 0;
    if (cbufstr->n < cbufstr->sz)
        cbufstr->n++;
}

str_t * cbufstr_top(cbufstr_t *cbufstr)
{
    if (cbufstr->n == 0) {
        fprintf(stderr, "Error: Tried to access empty cbufstr\n");
        exit(EXIT_FAILURE);
    }

    size_t nxt = cbufstr->nxt;
    size_t last = 0;

    if (nxt == 0)
        last = cbufstr->sz - 1;
    else
        last = cbufstr->nxt - 1;

    return cbufstr->buf[last];
}

str_t * cbufstr_get(const cbufstr_t *cbufstr, size_t pos)
{
    if (cbufstr->n == 0) {
        fprintf(stderr, "Error: Tried to access empty cbufstr\n");
        exit(EXIT_FAILURE);
    }
    if (pos > (cbufstr->n - 1)) {
        fprintf(stderr, "Error: Not enough elements in cbufstr\n");
        exit(EXIT_FAILURE);
    }

    return cbufstr->buf[pos];
}


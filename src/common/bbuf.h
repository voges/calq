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

#ifndef BBUF_H
#define BBUF_H

#include <stdint.h>
#include <stdlib.h>

typedef struct bbuf_t_ {
    unsigned char *bytes;
    size_t        sz;
} bbuf_t;

bbuf_t * bbuf_new(void);
void bbuf_free(bbuf_t *bbuf);
void bbuf_clear(bbuf_t *bbuf);
void bbuf_reserve(bbuf_t *bbuf, const size_t sz);
void bbuf_extend(bbuf_t *bbuf, const size_t ex);
void bbuf_trunc(bbuf_t *bbuf, const size_t tr);
void bbuf_append_bbuf(bbuf_t *bbuf, const bbuf_t *app);
void bbuf_append_byte(bbuf_t *bbuf, const unsigned char byte);
void bbuf_append_uint64(bbuf_t *bbuf, const uint64_t x);
void bbuf_append_buf(bbuf_t *bbuf, const unsigned char *buf, const size_t n);

#endif // BBUF_H


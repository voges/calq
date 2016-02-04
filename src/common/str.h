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

#ifndef STR_H
#define STR_H

#include <inttypes.h>
#include <stdint.h>
#include <stdlib.h>

typedef struct str_t_ {
    char   *s;  // Null-terminated string
    size_t len; // Length of s
    size_t sz;  // Bytes allocated for s
} str_t;

str_t * str_new(void);
void str_free(str_t *str);
void str_clear(str_t *str);
void str_reserve(str_t *str, const size_t sz);
void str_extend(str_t *str, const size_t ex);
void str_trunc(str_t *str, const size_t tr);
void str_append_str(str_t *str, const str_t *app);
void str_append_cstr(str_t *str , const char *cstr);
void str_append_cstrn(str_t *str, const char *cstr, const size_t len);
void str_append_char(str_t *str, const char c);
void str_append_int(str_t *str, const int64_t num);
void str_append_double2(str_t *str, const double dbl);
void str_copy_str(str_t *dest, const str_t *src);
void str_copy_cstr(str_t *dest, const char *src);

#endif // STR_H


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

#include "str.h"
#include <stdio.h>
#include <string.h>

static void str_init(str_t *str)
{
    str->s = NULL;
    str->len = 0;
    str->sz = 0;

    // Init string with terminating null byte
    str_reserve(str, 1);
    str->s[str->len] = '\0';
}

str_t * str_new(void)
{
    str_t *str = (str_t *)malloc(sizeof(str_t));
    if (!str) abort();
    str_init(str);
    return str;
}

void str_free(str_t *str)
{
    if (str != NULL) {
        if (str->s != NULL) {
            free(str->s);
            str->s = NULL;
        }
        free(str);
        str = NULL;
    } else {
        fprintf(stderr, "Error: Tried to free null pointer\n");
        exit(EXIT_FAILURE);
    }
}

void str_clear(str_t *str)
{
    if (str->s != NULL) {
        free(str->s);
        str->s = NULL;
    }
    str_init(str);
}

void str_reserve(str_t *str, const size_t sz)
{
    str->sz = sz;
    str->s = (char *)realloc(str->s, str->sz * sizeof(char));
    if (!(str->s)) abort();
}

void str_extend(str_t *str, const size_t ex)
{
    str_reserve(str, str->sz + ex);
}

void str_trunc(str_t *str, const size_t tr)
{
    str->len -= tr;
    str->sz -= tr;
    str_reserve(str, str->sz);
    str->s[str->len] = '\0';
}

void str_append_str(str_t *str, const str_t *app)
{
    str_extend(str, app->len);
    memcpy(str->s + str->len, app->s, app->len);
    str->len += app->len;
    str->s[str->len] = '\0';
}

void str_append_cstr(str_t *str, const char *cstr)
{
    size_t len = strlen(cstr);
    str_extend(str, len);
    memcpy(str->s + str->len, cstr, len);
    str->len += len;
    str->s[str->len] = '\0';
}

void str_append_cstrn(str_t *str, const char *cstr, const size_t len)
{
    if (len > strlen(cstr)) {
        fprintf(stderr, "Error: Failed to append C-string\n");
        exit(EXIT_FAILURE);
    }
    str_extend(str, len);
    memcpy(str->s + str->len, cstr, len);
    str->len += len;
    str->s[str->len] = '\0';
}

void str_append_char(str_t *str, const char c)
{
    str_extend(str, 1);
    str->s[str->len++] = c;
    str->s[str->len] = '\0';
}

void str_append_int(str_t *str, const int64_t num)
{
    char num_cstr[101];
    snprintf(num_cstr, sizeof(num_cstr), "%"PRId64, num);
    str_append_cstr(str, num_cstr);
}

void str_append_double2(str_t *str, const double dbl)
{
    char dbl_cstr[101];
    snprintf(dbl_cstr, sizeof(dbl_cstr), "%.2f", dbl);
    str_append_cstr(str, dbl_cstr);
}

void str_copy_str(str_t *dest, const str_t *src)
{
    str_clear(dest);
    str_reserve(dest, src->len + 1);
    memcpy(dest->s, src->s, src->len);
    dest->len = src->len;
    dest->s[dest->len] = '\0';
}

void str_copy_cstr(str_t *str, const char *cstr)
{
    str_clear(str);
    size_t len = strlen(cstr);
    str_reserve(str, len + 1);
    memcpy(str->s, cstr, len);
    str->len = len;
    str->s[str->len] = '\0';
}


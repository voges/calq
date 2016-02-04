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

/*
 * Copyright (c) 2013, 2014, Simon Fraser University, Indiana University
 * Bloomington. All rights reserved. Redistribution and use in source and
 * binary forms, with or without modification, are permitted provided that the
 * following conditions are met:
 *
 *   * Redistributions of source code must retain the above copyright notice,
 *     this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in the
 *     documentation and/or other materials provided with the distribution.
 *   * Neither the name of the Simon Fraser University, Indiana University
 *     Bloomington nor the names of its contributors may be used to endorse
 *     or promote products derived from this software without specific prior
 *     written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS “AS IS”
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include "samparser.h"
#include "tsclib.h"
#include <string.h>

static void samparser_init(samparser_t *samparser, FILE *fp)
{
    samparser->fp = fp;
}

samparser_t * samparser_new(FILE *fp)
{
    samparser_t *samparser = (samparser_t *)osro_malloc(sizeof(samparser_t));
    samparser->head = str_new();
    samparser_init(samparser, fp);
    return samparser;
}

void samparser_free(samparser_t *samparser)
{
    if (samparser != NULL) {
        str_free(samparser->head);
        free(samparser);
        samparser = NULL;
    } else {
        tsc_error("Tried to free null pointer\n");
    }
}

void samparser_head(samparser_t *samparser)
{
    // Read the SAM header
    bool samheader = false;
    while (fgets(samparser->curr.line, sizeof(samparser->curr.line),
                 samparser->fp)) {
        if (*(samparser->curr.line) == '@') {
            str_append_cstr(samparser->head, samparser->curr.line);
            samheader = true;
        } else {
            if (!samheader) tsc_error("SAM header missing\n");
            size_t offset = -strlen(samparser->curr.line);
            fseek(samparser->fp, (long)offset, SEEK_CUR);
            break;
        }
    }
}

static void samparser_parse(samparser_t *samparser)
{
    size_t l = strlen(samparser->curr.line) - 1;

    while (l && (samparser->curr.line[l] == '\r'
               || samparser->curr.line[l] == '\n'))
        samparser->curr.line[l--] = '\0';

    char *c = samparser->curr.qname = samparser->curr.line;
    int f = 1;

    while (*c) {
        if (*c == '\t') {
            if (f ==  1) samparser->curr.flag  = (uint16_t)atoi(c + 1);
            if (f ==  2) samparser->curr.rname = c + 1;
            if (f ==  3) samparser->curr.pos   = (uint32_t)atoi(c + 1);
            if (f ==  4) samparser->curr.mapq  = (uint8_t)atoi(c + 1);
            if (f ==  5) samparser->curr.cigar = c + 1;
            if (f ==  6) samparser->curr.rnext = c + 1;
            if (f ==  7) samparser->curr.pnext = (uint32_t)atoi(c + 1);
            if (f ==  8) samparser->curr.tlen  = (int64_t)atoi(c + 1);
            if (f ==  9) samparser->curr.seq   = c + 1;
            if (f == 10) samparser->curr.qual  = c + 1;
            if (f == 11) samparser->curr.opt   = c + 1;
            f++;
            *c = '\0';
            if (f == 12) break;
        }
        c++;
    }

    if (f == 11) samparser->curr.opt = c;
}

bool samparser_next(samparser_t *samparser)
{
    // Try to read and parse next line
    if (fgets(samparser->curr.line, sizeof(samparser->curr.line),
              samparser->fp)) {
        if (*(samparser->curr.line) == '@')
            tsc_error("Tried to read SAM record but found header line\n");
        else
            samparser_parse(samparser);
    } else {
        return false;
    }
    return true;
}


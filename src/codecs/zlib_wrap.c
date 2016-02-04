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

#include "zlib.h"
#include <stdio.h>
#include <stdlib.h>

unsigned char * zlib_compress(unsigned char *in,
                              size_t        in_sz,
                              size_t        *out_sz)
{
    Byte *out;
    compressBound(in_sz);
    *out_sz = compressBound(in_sz) + 1;
    out = (Byte *)calloc((uInt)*out_sz, 1);
    int err = compress(out, out_sz, (const Bytef *)in, (uLong)in_sz);
    if (err != Z_OK) {
        fprintf(stderr, "Error: zlib failed to compress: %d\n", err);
        exit(EXIT_FAILURE);
    }
    return out;
}

unsigned char * zlib_decompress(unsigned char *in,
                                size_t        in_sz,
                                size_t        out_sz)
{
    Bytef *out = (Bytef *)malloc(out_sz);
    if (!out) abort();
    int err = uncompress(out, (uLongf *)&out_sz, (const Bytef *)in, (uLong)in_sz);
    if (err != Z_OK) {
        fprintf(stderr, "Error: zlib failed to uncompress: %d\n", err);
        exit(EXIT_FAILURE);
    }
    return out;
}


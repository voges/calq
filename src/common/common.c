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

#include "common.h"

long tvdiff(struct timeval tv0, struct timeval tv1)
{
    return (tv1.tv_sec - tv0.tv_sec) * 1000000 + tv1.tv_usec - tv0.tv_usec;
}

size_t ndigits(int64_t x)
{
    // Ugly but fast
    size_t n = 0;
    if (x < 0) n++;
    x = llabs(x);

    if (x < 10) return n+1;
    if (x < 100) return n+2;
    if (x < 1000) return n+3;
    if (x < 10000) return n+4;
    if (x < 100000) return n+5;
    if (x < 1000000) return n+6;
    if (x < 10000000) return n+7;
    if (x < 100000000) return n+8;
    if (x < 1000000000) return n+9;
    if (x < 10000000000) return n+10;
    if (x < 100000000000) return n+11;
    if (x < 1000000000000) return n+12;
    if (x < 10000000000000) return n+13;
    if (x < 100000000000000) return n+14;
    if (x < 1000000000000000) return n+15;
    if (x < 10000000000000000) return n+16;
    if (x < 100000000000000000) return n+17;
    if (x < 1000000000000000000) return n+18;
    return n+19; /* INT64_MAX: 2^63 - 1 = 9223372036854775807 */
}


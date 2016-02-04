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

#ifndef TSC_SAMCODEC_H
#define TSC_SAMCODEC_H

#include "samparser.h"
#include "codecs/auxcodec.h"
#include "codecs/idcodec.h"
//#include "codecs/nuccodec_o0.h"
#include "codecs/nuccodec_o1.h"
#include "codecs/paircodec.h"
#include "codecs/qualcodec.h"
#include "common/str.h"
#include <stdio.h>

typedef struct samcodec_t_ {
    FILE         *ifp;
    FILE         *ofp;
    unsigned int blk_sz;
    samparser_t  *samparser;
    auxcodec_t   *auxcodec;
    idcodec_t    *idcodec;
    nuccodec_t   *nuccodec;
    paircodec_t  *paircodec;
    qualcodec_t  *qualcodec;
} samcodec_t;

samcodec_t * samcodec_new(FILE *ifp, FILE *ofp, unsigned int blk_sz);
void samcodec_free(samcodec_t *samcodec);
void samcodec_encode(samcodec_t *samcodec);
void samcodec_decode(samcodec_t *samcodec);
void samcodec_info(samcodec_t *samcodec);

#endif // TSC_SAMCODEC_H


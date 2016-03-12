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

//
// Wrapper functions to safely read/write different data types from/to files (or
// streams).
// The 'osro_fwrite_uintXX' resp. 'osro_fread_uintXX' functions are compatible
// with each other. They are independent from the endianness of the system this
// code is built on.
//

#ifndef OSRO_IO_H
#define OSRO_IO_H

#include <stdint.h>
#include <stdio.h>

FILE * osro_fopen(const char *fname, const char *mode);
void osro_fclose(FILE *fp);

size_t osro_fwrite_byte(FILE *fp, const unsigned char byte);
size_t osro_fwrite_buf(FILE *fp, const unsigned char *buf, const size_t n);
size_t osro_fwrite_uint32(FILE *fp, const uint32_t dword);
size_t osro_fwrite_uint64(FILE *fp, const uint64_t qword);

size_t osro_fread_byte(FILE *fp, unsigned char *byte);
size_t osro_fread_buf(FILE *fp, unsigned char *buf, const size_t n);
size_t osro_fread_uint32(FILE *fp, uint32_t *dword);
size_t osro_fread_uint64(FILE *fp, uint64_t *qword);

#endif // OSRO_IO_H


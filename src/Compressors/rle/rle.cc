/** @file rle.cc
 *  @brief This files contains the implementation of the RLE codec.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#include "rle.h"

#include <iostream>

#include "Common/helpers.h"
#include "Common/Exceptions.h"
#include "Common/os.h"

#if defined(OS_WINDOWS)
    #include <windows.h>
#elif defined(OS_APPLE) || defined(OS_LINUX)
    #include <unistd.h> // sysconf(3)
#else
    #error "Operating system not supported"
#endif

// TODO: Let pointers returned by malloc point to addresses at the start of a
//       page and verify that realloc just extends the continuous memory
//       block

//
// Decoder specification:
// ----------------------
// A byte can have 256 values: from 0x00 = 0 up to 0xFF = 255. By definition
// the smallest run-length is 3. Therefore, with one byte, we can express
// run-lengths from 3 (= 0x00) up to 258 (= 0xFF). However, we define the
// longest possible run-length as 255.
// Assume, we only want to encode N symbols 0...(N-1). We assign the N byte
// values 0x00 - 0x** = (N-1) to our symbols. Then, we have (256-N) byte values
// left to encode run-lengths from 3 up to 255-N.
//
// Example for N = 4 symbols 0,1,2,3. We can encode 256 - 4 = 252 run-lengths
// ranging from 3 up to 254:
//   Symbol 0 is encoded as 0x00
//   Symbol 1 is encoded as 0x01
//   Symbol 2 is encoded as 0x02
//   Symbol 3 is encoded as 0x03
//   Run-length 3 is encoded as 0x04
//   Run-length 4 is encoded as 0x05
//   ...
//   Run-length 254 is encoded as 0xFF
//

static const unsigned char RLE_MIN_RUN_LENGTH = 3;
static const unsigned char RLE_MAX_RUN_LENGTH = 255;

unsigned char * rle_encode(unsigned char       *in,
                           const size_t        in_sz,
                           size_t              *out_sz,
                           const unsigned char N,
                           const unsigned char offset)
{
#if defined(OS_WINDOWS)
    SYSTEM_INFO si;
    GetSystemInfo(&si);
    const unsigned long pageSize = (unsigned long)si.dwPageSize;
#elif defined(OS_APPLE) || defined(OS_LINUX)
    const long pageSize = sysconf(_SC_PAGESIZE); // _SC_PAGE_SIZE is OK, too
#else
    #error "Operating system not supported"
#endif

    if (in_sz == 0) {
        throwErrorException("in_sz = 0, nothing to encode");
    }

    // Less than 2 symbols do not make sense, the same holds for more than
    // 255 symbols, because we need at least one byte to represent a run-length
    // of 3.
    if (N < 2/* || N > 255*/) {
        throwErrorException("N out of range");
    }

    // We have at least 2 symbols. The highest possible values for these are
    // 0xFE = 254 and 0xFF = 255, respectively. Thus, an offset greater than
    // 254 cannot make sense.
    if (offset > 254) {
        throwErrorException("offset out of range");
    }

    const unsigned char maxRunLength = RLE_MAX_RUN_LENGTH - N;

    // Allocate one page for 'out'
    size_t out_alloc_sz = pageSize;
    unsigned char *out = (unsigned char *)malloc(out_alloc_sz);
    if (out == NULL) { throwErrorException("malloc failed"); }

    size_t in_itr = 0;
    size_t out_itr = 0;

    unsigned char prevChar = in[in_itr++] - offset;
    if (prevChar > N-1) {
        throwErrorException("Symbol in stream out of range");
    }
    unsigned char n = 1;

    while (in_itr < in_sz) {
        unsigned char currChar = in[in_itr++] - offset;

        if (currChar > N-1) {
            throwErrorException("Symbol in stream out of range");
        }

        if ((currChar != prevChar) || (n == maxRunLength)) {
            // Allocate additional page if needed
            if (out_alloc_sz-out_itr-1 < RLE_MIN_RUN_LENGTH) {
                out_alloc_sz = out_itr + pageSize;
                out = (unsigned char*)realloc(out, out_alloc_sz);
                if (out == NULL) { throwErrorException("realloc failed"); }
            }

            if (n < RLE_MIN_RUN_LENGTH) {
                for (unsigned char j = 0; j < n; j++) {
                    out[out_itr++] = prevChar;
                }
            } else {
                out[out_itr++] = n + N;
                out[out_itr++] = prevChar;
            }
            n = 0;
        }

        prevChar = currChar;
        n++;
    }

    // Allocate additional page if needed
    if (out_alloc_sz-out_itr-1 < RLE_MIN_RUN_LENGTH) {
        out_alloc_sz = out_itr + pageSize;
        out = (unsigned char*)realloc(out, out_alloc_sz);
        if (out == NULL) { throwErrorException("realloc failed"); }
    }

    // Write out last run
    if (n < RLE_MIN_RUN_LENGTH) {
        for (unsigned char j = 0; j < n; j++) {
            out[out_itr++] = prevChar;
        }
    } else {
        out[out_itr++] = n + N;
        out[out_itr++] = prevChar;
    }

    // Realloc out to fit out_sz
    *out_sz = out_itr;
    out = (unsigned char *)realloc(out, *out_sz);
    if (out == NULL) { throwErrorException("realloc failed"); }

    return out;
}

unsigned char * rle_decode(unsigned char       *in,
                           const size_t        in_sz,
                           size_t              *out_sz,
                           const unsigned char N,
                           const unsigned char offset)
{
#if defined(OS_WINDOWS)
    SYSTEM_INFO si;
    GetSystemInfo(&si);
    const unsigned long pageSize = (unsigned long)si.dwPageSize;
#elif defined(OS_APPLE) || defined(OS_LINUX)
    const long pageSize = sysconf(_SC_PAGESIZE); // _SC_PAGE_SIZE is OK, too
#else
    #error "Operating system not supported"
#endif

    if (in_sz == 0) {
        throwErrorException("in_sz = 0, nothing to decode");
    }

    // Less than 2 symbols do not make sense, the same holds for more than
    // 255 symbols, because we need at least one byte to represent a run-length
    // of 3.
    if (N < 2/* || N > 255*/) {
        throwErrorException("N out of range");
    }

    // We have at least 2 symbols. The highest possible values for these are
    // 0xFE = 254 and 0xFF = 255, respectively. Thus, an offset greater than
    // 254 cannot make sense.
    if (offset > 254) {
        throwErrorException("offset out of range");
    }

    // Allocate one page for 'out'
    size_t out_alloc_sz = pageSize;
    unsigned char *out = (unsigned char *)malloc(out_alloc_sz);
    if (out == NULL) { throwErrorException("malloc failed"); }

    size_t in_itr = 0;
    size_t out_itr = 0;
    while (in_itr < in_sz) {
        // Allocate additional page if needed
        if (out_alloc_sz-out_itr-1 < RLE_MAX_RUN_LENGTH) {
            out_alloc_sz = out_itr + pageSize;
            out = (unsigned char*)realloc(out, out_alloc_sz);
            if (out == NULL) { throwErrorException("realloc failed"); }
        }

        unsigned char currChar = in[in_itr++];

        // Check, if currChar contains a run-length
        if (currChar > (unsigned char)(N-1)) {
            unsigned char runSymbol = in[in_itr++] + offset;
            unsigned char runLength = currChar - (unsigned char)N;
            for (size_t j = 0; j < runLength; j++) {
                out[out_itr++] = runSymbol;
            }
        } else {
            out[out_itr++] = currChar + offset;
        }
    }

    // Realloc out to fit out_sz
    *out_sz = out_itr;
    out = (unsigned char *)realloc(out, *out_sz);
    if (out == NULL) { throwErrorException("realloc failed"); }

    return out;
}


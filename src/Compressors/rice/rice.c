/** @file rice.c
 *  @brief This file contains the implementation of the Rice codec.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "Compressors/rice/rice.h"
#include "Common/os_config.h"
#include <stdbool.h>
#include <string.h>

#if defined(CQ_OS_WINDOWS)
    #include <windows.h>
#elif defined(CQ_OS_APPLE) || defined(CQ_OS_LINUX)
    #include <unistd.h> /* sysconf(3) */
#else
    #error "Operating system not supported"
#endif

// TODO: Let pointers returned by malloc point to addresses at the start of a
//       page and verify that realloc just extends the continuous memory
//       block

#ifndef SIZE_MAX
#define SIZE_MAX (size_t)-1 // due to portability
#endif

typedef struct ricecodec_t_ {
    unsigned char bbuf;        // bit buffer
    size_t        bbuf_filled; // bit count
    size_t        in_idx;      // current 'in' index
    size_t        out_idx;     // current 'out' index
    size_t        in_sz;       // 'in' size
    bool          eof;         // flag to indicate EOF
} ricecodec_t;

static void ricecodec_init(ricecodec_t *rc, size_t in_sz)
{
    rc->bbuf = 0x00;
    rc->bbuf_filled = 0;
    rc->in_idx = 0;
    rc->out_idx = 0;
    rc->in_sz = in_sz;
    rc->eof = false;
}

// Encoder
// -----------------------------------------------------------------------------

static int codelen(unsigned char x, int k)
{
    int m = 1 << k;
    int q = x / m;
    return (q + 1 + k);
}

static void put_bit(ricecodec_t *rc, unsigned char *out, unsigned char b)
{
    rc->bbuf = rc->bbuf | ((b & 1) << rc->bbuf_filled);

    if (rc->bbuf_filled == 7) {
        out[rc->out_idx++] = rc->bbuf;
        rc->bbuf = 0x00;
        rc->bbuf_filled = 0;
    } else {
        rc->bbuf_filled++;
    }
}

static void encode(ricecodec_t   *rc,
                   unsigned char *out,
                   unsigned char x,
                   int           k)
{
    int m = 1 << k;
    int q = x / m;
    int i = 0;

    for (i = 0; i < q; i++)
        put_bit(rc, out, 1);

    put_bit(rc, out, 0);

    for (i = (k - 1); i >= 0; i--)
        put_bit(rc, out, (x >> i) & 1);
}

unsigned char * rice_compress(unsigned char *in,
                              const size_t  in_sz,
                              size_t        *out_sz)
{
    ricecodec_t rc;
    ricecodec_init(&rc, in_sz);

    // Find best Rice parameter k
    unsigned int k = 0;
    unsigned int k_best = 0;
    size_t bit_cnt = 0;
    size_t bit_cnt_best = SIZE_MAX;
    for (k = 0; k < 8; k++) {
        size_t i = 0;
        for (i = 0; i < in_sz; i++) bit_cnt += codelen(in[i], k);
        if (bit_cnt < bit_cnt_best) {
            bit_cnt_best = bit_cnt;
            k_best = k;
        }
    }

    // Compute number of bytes needed
    *out_sz = bit_cnt_best / 8;
    unsigned int bit_rem = (bit_cnt_best % 8) + 3; // +3 bits for k in [0-7]
    *out_sz += (bit_rem / 8) + 1;

    // Allocate enough memory for 'out'
    unsigned char *out = (unsigned char *)malloc(*out_sz);
    if (!out) abort();
    memset(out, 0x00, *out_sz);

    // Output k
    put_bit(&rc, out, (k_best >> 2) & 1);
    put_bit(&rc, out, (k_best >> 1) & 1);
    put_bit(&rc, out, (k_best     ) & 1);

    // Encode
    size_t i = 0;
    for (i = 0; i < in_sz; i++)
        encode(&rc, out, in[i], k_best);

    // Flush bit buffer: fill with 1s, so decoder will run into EOF
    size_t f = 8 - rc.bbuf_filled;
    while (f--)
        put_bit(&rc, out, 1);

    return out;
}

// Decoder
// -----------------------------------------------------------------------------

static int get_bit(ricecodec_t *rc, unsigned char *in)
{
    if (!(rc->bbuf_filled)) {
        if (rc->in_idx < rc->in_sz) {
            rc->bbuf = in[rc->in_idx++];
            rc->bbuf_filled = 8;
        } else {
            rc->eof = true; // reached EOF
        }
    }

    int tmp = rc->bbuf & 0x01;
    rc->bbuf = rc->bbuf >> 1;
    rc->bbuf_filled--;

    return tmp;
}

unsigned char * rice_decompress(unsigned char *in,
                                const size_t  in_sz,
                                size_t        *out_sz)
{
#if defined(CQ_OS_WINDOWS)
    SYSTEM_INFO si;
    GetSystemInfo(&si);
    const unsigned long pageSize = (unsigned long)si.dwPageSize;
#elif defined(CQ_OS_APPLE) || defined(CQ_OS_LINUX)
    const long pageSize = sysconf(_SC_PAGESIZE); // _SC_PAGE_SIZE is OK, too
#else
    #error "Operating system not supported"
#endif

    ricecodec_t rc;
    ricecodec_init(&rc, in_sz);
    *out_sz = 0;

    // Allocate enough memory (10 is a bad estimate) for 'out'
    size_t out_alloc = 10 * in_sz;
    unsigned char *out = (unsigned char *)malloc(out_alloc);
    if (!out) abort();

    unsigned int k = (get_bit(&rc, in) << 2) |
                     (get_bit(&rc, in) << 1) |
                     (get_bit(&rc, in)     );

    while (1) {
        int m = 1 << k, q = 0, x = 0, i = 0;

        while (get_bit(&rc, in))
            q++;

        if (rc.eof) break;

        x = m * q;

        for (i = (k - 1); i >= 0; i--)
            x = x | (get_bit(&rc, in) << i);

        out[rc.out_idx++] = x;
        (*out_sz)++;

        // Allocate additional page if needed
        if (rc.out_idx == out_alloc) {
            out_alloc = *out_sz + pageSize;
            out = (unsigned char*)realloc(out, out_alloc);
            if (!out) abort();
        }
    }

    out = (unsigned char *)realloc(out, *out_sz);
    if (!out) abort();

    return out;
}


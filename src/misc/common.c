#include "misc/common.h"
#include <stdio.h>

bool yesno(void)
{
    int c = getchar();
    bool yes = c == 'y' || c == 'Y';
    while (c != '\n' && c != EOF)
        c = getchar();
    return yes;
}

long tvdiff(const struct timeval tv0, const struct timeval tv1)
{
    return (tv1.tv_sec - tv0.tv_sec) * 1000000 + tv1.tv_usec - tv0.tv_usec;
}

size_t ndigits(const int64_t x)
{
    // Ugly but fast
    size_t n = 0;
    if (x < 0) n++;
    int64_t tmp = llabs(x);

    if (tmp < 10) return n+1;
    if (tmp < 100) return n+2;
    if (tmp < 1000) return n+3;
    if (tmp < 10000) return n+4;
    if (tmp < 100000) return n+5;
    if (tmp < 1000000) return n+6;
    if (tmp < 10000000) return n+7;
    if (tmp < 100000000) return n+8;
    if (tmp < 1000000000) return n+9;
    if (tmp < 10000000000) return n+10;
    if (tmp < 100000000000) return n+11;
    if (tmp < 1000000000000) return n+12;
    if (tmp < 10000000000000) return n+13;
    if (tmp < 100000000000000) return n+14;
    if (tmp < 1000000000000000) return n+15;
    if (tmp < 10000000000000000) return n+16;
    if (tmp < 100000000000000000) return n+17;
    if (tmp < 1000000000000000000) return n+18;
    return n+19; /* INT64_MAX: 2^63 - 1 = 9223372036854775807 */
}


#include "common.h"
#include <stdio.h>

bool yesno(void)
{
    int c = getchar();
    bool yes = c == 'y' || c == 'Y';
    while (c != '\n' && c != EOF)
        c = getchar();
    return yes;
}

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


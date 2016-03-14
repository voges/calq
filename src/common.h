#ifndef CQ_COMMON_H
#define CQ_COMMON_H

#include <stdbool.h>
#include <stdlib.h>

// constants
#define KB 1000LL
#define MB (KB*1000LL)
#define GB (MB*1000LL)

// safe debug macro
#if DBG
    #define DEBUG(c,...)\
        do {\
            fprintf(stderr, "%s:%d: %s: "c, __FILE__, __LINE__, \
                    __FUNCTION__, ##__VA_ARGS__);\
        } while (false)
#else
    #define DEBUG(c,...) do { } while (false)
#endif

bool yesno(void);
long tvdiff(struct timeval tv0, struct timeval tv1);
size_t ndigits(int64_t x);

#endif // CQ_COMMON_H


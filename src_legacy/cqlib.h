#ifndef CQ_CQLIB_H
#define CQ_CQLIB_H

#define _GNU_SOURCE

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

// constants


// safe debug macro
//#if DBG
//    #define DEBUG(c,...)\
//        do {\
//            fprintf(stderr, "%s:%d: %s: "c, __FILE__, __LINE__, \
//                    __FUNCTION__, ##__VA_ARGS__);\
//        } while (false)
//#else
//    #define DEBUG(c,...) do { } while (false)
//#endif

#define CQ_SUCCESS 0
#define CQ_FAILURE 1

#define CQ_COLOR_RED     "\x1b[31m"
#define CQ_COLOR_GREEN   "\x1b[32m"
#define CQ_COLOR_YELLOW  "\x1b[33m"
#define CQ_COLOR_BLUE    "\x1b[34m"
#define CQ_COLOR_MAGENTA "\x1b[35m"
#define CQ_COLOR_CYAN    "\x1b[36m"
#define CQ_COLOR_RESET   "\x1b[0m"

typedef struct cq_opts_t_ {
    size_t blocksz;
    const char *fname_in;
    const char *fname_out;
    bool force;
    enum { 
        CQ_OPTS_MODE_COMPRESS, 
        CQ_OPTS_MODE_DECOMPRESS, 
        CQ_OPTS_MODE_INFO 
    } mode;
} cq_opts_t;

void cq_out(const char *fmt, ...);
void cq_err(const char *fmt, ...);

void * cq_malloc(const size_t size);
void * cq_calloc(const size_t nmemb, const size_t size);
void * cq_realloc(void *ptr, const size_t size);
void cq_free(void *ptr);

FILE * cq_fopen(const char *fname, const char * const mode);
void cq_fclose(FILE *fp);

size_t cq_fwrite_byte(FILE *fp, const unsigned char byte);
size_t cq_fwrite_buf(FILE *fp, const unsigned char *buf, const size_t n);
size_t cq_fwrite_uint32(FILE *fp, const uint32_t dword);
size_t cq_fwrite_uint64(FILE *fp, const uint64_t qword);

size_t cq_fread_byte(FILE *fp, unsigned char *byte);
size_t cq_fread_buf(FILE *fp, unsigned char *buf, const size_t n);
size_t cq_fread_uint32(FILE *fp, uint32_t *dword);
size_t cq_fread_uint64(FILE *fp, uint64_t *qword);

#endif // CQ_CQLIB_H


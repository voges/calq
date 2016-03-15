#ifndef CQ_CQLIB_H
#define CQ_CQLIB_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

void cq_log(const char *fmt, ...);
void cq_error(const char *fmt, ...);

void * cq_malloc(const size_t size);
void * cq_calloc(const size_t nmemb, const size_t size);
void * cq_realloc(void *ptr, const size_t size);
void cq_free(void *ptr);

FILE * cq_fopen(const char *fname, const char *mode);
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


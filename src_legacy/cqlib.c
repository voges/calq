//#define _GNU_SOURCE

#include "cqlib.h"
#include <stdarg.h>

void cq_out(const char *fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    char *msg;
    vasprintf(&msg, fmt, args);
    va_end(args);
    fprintf(stdout, "%s", msg);
    free(msg);
}

void cq_err(const char *fmt, ...)
{
    va_list args;
    va_start(args, fmt);
    char *msg;
    vasprintf(&msg, fmt, args);
    va_end(args);
    fprintf(stderr, CQ_COLOR_RED "%s" CQ_COLOR_RESET, msg);
    free(msg);
}

void * cq_malloc(const size_t size)
{
    void *ptr = malloc(size);
    if (ptr == NULL) {
        cq_err("Could not allocate %zu bytes\n", size);
        exit(EXIT_FAILURE);
    }
    return ptr;
}

void * cq_calloc(const size_t nmemb, const size_t size)
{
    void *ptr = calloc(nmemb, size);
    if (ptr == NULL) {
        cq_err("Could not allocate %zu bytes\n", nmemb*size);
        exit(EXIT_FAILURE);
    }
    return ptr;
}

void * cq_realloc(void *ptr, const size_t size)
{
    void *p = realloc(ptr, size);
    if (p == NULL) {
        cq_err("Could not allocate %zu bytes\n", size);
        exit(EXIT_FAILURE);
    }
    return p;
}

void cq_free(void *ptr)
{
    if (ptr != NULL) {
        free(ptr);
        ptr = NULL;
    } else {
        cq_err("Tried to free null pointer\n");
        exit(EXIT_FAILURE);
    }
}

FILE * cq_fopen(const char* fname, const char*const mode)
{
    FILE *fp = fopen(fname, mode);
    if (fp == NULL) {
        fclose(fp);
        cq_err("Failed to open file: %s\n", fname);
        exit(EXIT_FAILURE);
    }
    return fp;
}

void cq_fclose(FILE *fp)
{
    if (fp != NULL) {
        fclose(fp);
        fp = NULL;
    } else {
        cq_err("Failed to close file\n");
        exit(EXIT_FAILURE);
    }
}

size_t cq_fwrite_byte(FILE *fp, const unsigned char byte)
{
    if (fwrite(&byte, 1, 1, fp) != 1) {
        cq_err("Could not write byte\n");
        exit(EXIT_FAILURE);
    }
    return 1;
}

size_t cq_fwrite_buf(FILE *fp, const unsigned char *buf, const size_t n)
{
    if (fwrite(buf, 1, n, fp) != n) {
        cq_err("Could not write %zu byte(s)\n", n);
        exit(EXIT_FAILURE);
    }
    return n;
}

size_t cq_fwrite_uint32(FILE *fp, const uint32_t dword)
{
    cq_fwrite_byte(fp, (unsigned char)(dword >> 24) & 0xFF);
    cq_fwrite_byte(fp, (unsigned char)(dword >> 16) & 0xFF);
    cq_fwrite_byte(fp, (unsigned char)(dword >>  8) & 0xFF);
    cq_fwrite_byte(fp, (unsigned char)(dword      ) & 0xFF);
    return sizeof(uint32_t);
}

size_t cq_fwrite_uint64(FILE *fp, const uint64_t qword)
{
    cq_fwrite_byte(fp, (unsigned char)(qword >> 56) & 0xFF);
    cq_fwrite_byte(fp, (unsigned char)(qword >> 48) & 0xFF);
    cq_fwrite_byte(fp, (unsigned char)(qword >> 40) & 0xFF);
    cq_fwrite_byte(fp, (unsigned char)(qword >> 32) & 0xFF);
    cq_fwrite_byte(fp, (unsigned char)(qword >> 24) & 0xFF);
    cq_fwrite_byte(fp, (unsigned char)(qword >> 16) & 0xFF);
    cq_fwrite_byte(fp, (unsigned char)(qword >>  8) & 0xFF);
    cq_fwrite_byte(fp, (unsigned char)(qword      ) & 0xFF);
    return sizeof(uint64_t);
}

size_t cq_fread_byte(FILE *fp, unsigned char *byte)
{
    return fread(byte, 1, 1, fp);
}

size_t cq_fread_buf(FILE *fp, unsigned char *buf, const size_t n)
{
    return fread(buf, 1, n, fp);
}

size_t cq_fread_uint32(FILE *fp, uint32_t *dword)
{
    unsigned char *bytes = (unsigned char *)cq_malloc(sizeof(uint32_t));
    size_t ret = fread(bytes, 1, sizeof(uint32_t), fp);

    if (ret != sizeof(uint32_t)) {
        free(bytes);
        return ret;
    }

    *dword = (uint32_t)bytes[0] << 24 |
             (uint32_t)bytes[1] << 16 |
             (uint32_t)bytes[2] <<  8 |
             (uint32_t)bytes[3];

    free(bytes);
    return ret;
}

size_t cq_fread_uint64(FILE *fp, uint64_t *qword)
{
    unsigned char *bytes = (unsigned char *)cq_malloc(sizeof(uint64_t));
    size_t ret = fread(bytes, 1, sizeof(uint64_t), fp);

    if (ret != sizeof(uint64_t)) {
        free(bytes);
        return ret;
    }

    *qword = (uint64_t)bytes[0] << 56 |
             (uint64_t)bytes[1] << 48 |
             (uint64_t)bytes[2] << 40 |
             (uint64_t)bytes[3] << 32 |
             (uint64_t)bytes[4] << 24 |
             (uint64_t)bytes[5] << 16 |
             (uint64_t)bytes[6] <<  8 |
             (uint64_t)bytes[7];

    free(bytes);
    return ret;
}


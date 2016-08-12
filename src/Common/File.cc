/** @file File.cc
 *  @brief This file contains the implementation of the File class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "File.h"
#include "Common/Exceptions.h"
#include "Common/os_config.h"
#include <limits.h>
#include <stdio.h>

File::File(void)
    : fp(NULL)
    , fsize(0)
{
    // empty
}

File::File(const std::string &path, const char *mode)
    : fp(NULL)
    , fsize(0)
{
    open(path, mode);
}

File::~File(void)
{
    close();
}

void File::open(const std::string &path, const char *mode)
{
#ifdef OS_WINDOWS
    int err = fopen_s(&fp, path.c_str(), mode);
    if (err != 0) {
        throwErrorException("Failed to open file");
    }
#else
    fp = fopen(path.c_str(), mode);
    if (fp == NULL) {
        throwErrorException("Failed to open file");
    }
#endif

    // Compute file size
    fseek(fp, 0, SEEK_END);
    fsize = ftell(fp);
    fseek(fp, 0, SEEK_SET);
}

void File::close(void) 
{ 
    if (fp != NULL) {
        fclose(fp);
        fp = NULL;
    } else {
        throwErrorException("Failed to close file");
    }
}

void File::advance(const size_t &offset)
{
    int ret = fseek(fp, (long int)offset, SEEK_CUR);
    if (ret != 0) {
        throwErrorException("fseek failed");
    }
}

bool File::eof(void) 
{
    int eof = feof(fp);
    return eof != 0 ? true : false;
}

void * File::handle(void)
{
    return fp;
}

void File::seek(const size_t &pos)
{
    if (pos > LONG_MAX) {
        throwErrorException("pos out of range");
    }
    int ret = fseek(fp, (long)pos, SEEK_SET);
    if (ret != 0) {
        throwErrorException("fseek failed");
    }
}

size_t File::size(void)
{
    return fsize;
}

size_t File::tell(void)
{
    long int offset = ftell(fp);
    if (offset == -1) {
        throwErrorException("ftell failed");
    }
    return offset;
}

size_t File::read(void *buffer, const size_t &size) 
{
    return fread(buffer, 1, size, fp); 
}

size_t File::write(void *buffer, const size_t &size) 
{
    return fwrite(buffer, 1, size, fp); 
}

size_t File::readByte(unsigned char *byte)
{
    return fread(byte, 1, 1, fp);
}

size_t File::readUint8(uint8_t *byte)
{
    return readByte(byte);
}

size_t File::readUint16(uint16_t *word)
{
    unsigned char *buffer = (unsigned char *)malloc(sizeof(uint16_t));
    if (buffer == NULL) {
        throwErrorException("malloc failed");
    }

    size_t ret = read(buffer, sizeof(uint16_t));

    if (ret != sizeof(uint16_t)) {
        free(buffer);
        throwErrorException("read failed");
    }

    *word = (uint16_t)buffer[2] <<  8 | (uint16_t)buffer[3];

    free(buffer);
    return ret;
}

size_t File::readUint32(uint32_t *dword)
{
    unsigned char *buffer = (unsigned char *)malloc(sizeof(uint32_t));
    if (buffer == NULL) {
        throwErrorException("malloc failed");
    }

    size_t ret = read(buffer, sizeof(uint32_t));

    if (ret != sizeof(uint32_t)) {
        free(buffer);
        throwErrorException("read failed");
    }

    *dword = (uint32_t)buffer[0] << 24 |
             (uint32_t)buffer[1] << 16 |
             (uint32_t)buffer[2] <<  8 |
             (uint32_t)buffer[3];

    free(buffer);
    return ret;
}

size_t File::readUint64(uint64_t *qword)
{
    unsigned char *buffer = (unsigned char *)malloc(sizeof(uint64_t));
    if (buffer == NULL) {
        throwErrorException("malloc failed");
    }

    size_t ret = read(buffer, sizeof(uint64_t));

    if (ret != sizeof(uint64_t)) {
        free(buffer);
        throwErrorException("read failed");
    }

    *qword = (uint64_t)buffer[0] << 56 |
             (uint64_t)buffer[1] << 48 |
             (uint64_t)buffer[2] << 40 |
             (uint64_t)buffer[3] << 32 |
             (uint64_t)buffer[4] << 24 |
             (uint64_t)buffer[5] << 16 |
             (uint64_t)buffer[6] <<  8 |
             (uint64_t)buffer[7];

    free(buffer);
    return ret;
}

size_t File::writeByte(const unsigned char &byte)
{
    if (fwrite(&byte, 1, 1, fp) != 1) {
        throwErrorException("fwrite failed");
    }
    return 1;
}

size_t File::writeUint8(const uint8_t &byte)
{
    return writeByte(byte);
}

size_t File::writeUint16(const uint16_t &word)
{
    size_t ret = 0;
    ret += writeByte((unsigned char)(word >> 8) & 0xFF);
    ret += writeByte((unsigned char)(word     ) & 0xFF);
    return ret;
}

size_t File::writeUint32(const uint32_t &dword)
{
    size_t ret = 0;
    ret += writeByte((unsigned char)(dword >> 24) & 0xFF);
    ret += writeByte((unsigned char)(dword >> 16) & 0xFF);
    ret += writeByte((unsigned char)(dword >>  8) & 0xFF);
    ret += writeByte((unsigned char)(dword      ) & 0xFF);
    return ret;
}

size_t File::writeUint64(const uint64_t &qword)
{
    size_t ret = 0;
    ret += writeByte((unsigned char)(qword >> 56) & 0xFF);
    ret += writeByte((unsigned char)(qword >> 48) & 0xFF);
    ret += writeByte((unsigned char)(qword >> 40) & 0xFF);
    ret += writeByte((unsigned char)(qword >> 32) & 0xFF);
    ret += writeByte((unsigned char)(qword >> 24) & 0xFF);
    ret += writeByte((unsigned char)(qword >> 16) & 0xFF);
    ret += writeByte((unsigned char)(qword >>  8) & 0xFF);
    ret += writeByte((unsigned char)(qword      ) & 0xFF);
    return ret;
}


/** @file File.cc
 *  @brief This file contains the implementation of the File class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "IO/File.h"
#include "Common/Exceptions.h"
#include "Common/os_config.h"
#include <limits.h>

cq::File::File(void)
    : m_fp(NULL)
    , m_fsize(0)
    , m_isOpen(false)
    , m_mode(File::MODE_READ)
{
    // empty
}

cq::File::File(const std::string &path, const Mode &mode)
    : m_fp(NULL)
    , m_fsize(0)
    , m_isOpen(false)
    , m_mode(mode)
{
    if (path.empty() == true) {
        throwErrorException("path is empty");
    }

    open(path, mode);
}

cq::File::~File(void)
{
    close();
}

void cq::File::open(const std::string &path, const Mode &mode)
{
    if (path.empty() == true) {
        throwErrorException("path is empty");
    }
    if (m_fp != NULL) { 
        throwErrorException("File pointer already in use");
    }

    const char *m;
    if (mode == MODE_READ) {
        m = "rb";
        m_mode = mode;
    } else if (mode == MODE_WRITE) {
        m = "wb";
        m_mode = mode;
    } else {
        throwErrorException("Unkown mode");
    }

#ifdef CQ_OS_WINDOWS
    int err = fopen_s(&m_fp, path.c_str(), m);
    if (err != 0) {
        throwErrorException("Failed to open file");
    }
#else
    m_fp = fopen(path.c_str(), m);
    if (m_fp == NULL) {
        throwErrorException("Failed to open file");
    }
#endif

    // Compute file size
    fseek(m_fp, 0, SEEK_END);
    m_fsize = ftell(m_fp);
    fseek(m_fp, 0, SEEK_SET);

    m_isOpen = true;
}

void cq::File::close(void) 
{
    if (m_isOpen == true) {
        if (m_fp != NULL) {
            fclose(m_fp);
            m_fp = NULL;
        } else {
            throwErrorException("Failed to close file");
        }
    }
}

void cq::File::advance(const size_t &offset)
{
    int ret = fseek(m_fp, (long int)offset, SEEK_CUR);
    if (ret != 0) {
        throwErrorException("fseek failed");
    }
}

bool cq::File::eof(void) const
{
    int eof = feof(m_fp);
    return eof != 0 ? true : false;
}

void * cq::File::handle(void) const
{
    return m_fp;
}

void cq::File::seek(const size_t &pos)
{
    if (pos > LONG_MAX) {
        throwErrorException("pos out of range");
    }
    int ret = fseek(m_fp, (long)pos, SEEK_SET);
    if (ret != 0) {
        throwErrorException("fseek failed");
    }
}

size_t cq::File::size(void) const
{
    return m_fsize;
}

size_t cq::File::tell(void) const
{
    long int offset = ftell(m_fp);
    if (offset == -1) {
        throwErrorException("ftell failed");
    }
    return offset;
}

size_t cq::File::read(void *buffer, const size_t &size) 
{
    if (size == 0) {
        throwErrorException("Attempted to read zero bytes");
    }
    return fread(buffer, 1, size, m_fp); 
}

size_t cq::File::write(void *buffer, const size_t &size) 
{
    if (size == 0) {
        throwErrorException("Attempted to write zero bytes");
    }
    return fwrite(buffer, 1, size, m_fp); 
}

size_t cq::File::readByte(unsigned char *byte)
{
    return fread(byte, 1, 1, m_fp);
}

size_t cq::File::readUint8(uint8_t *byte)
{
    return readByte(byte);
}

size_t cq::File::readUint16(uint16_t *word)
{
    unsigned char *buffer = (unsigned char *)malloc(sizeof(uint16_t));
    if (buffer == NULL) {
        throwErrorException("malloc failed");
    }

    size_t ret = read(buffer, sizeof(uint16_t));

    if (ret != sizeof(uint16_t)) {
        free(buffer);
        throwErrorException("read failed");
    } else {
        *word = (uint16_t)buffer[2] <<  8 | (uint16_t)buffer[3];
        free(buffer);
    }

    return ret;
}

size_t cq::File::readUint32(uint32_t *dword)
{
    unsigned char *buffer = (unsigned char *)malloc(sizeof(uint32_t));
    if (buffer == NULL) {
        throwErrorException("malloc failed");
    }

    size_t ret = read(buffer, sizeof(uint32_t));

    if (ret != sizeof(uint32_t)) {
        free(buffer);
        throwErrorException("read failed");
    } else {
        *dword = (uint32_t)buffer[0] << 24 |
                 (uint32_t)buffer[1] << 16 |
                 (uint32_t)buffer[2] <<  8 |
                 (uint32_t)buffer[3];
        free(buffer);
    }

    return ret;
}

size_t cq::File::readUint64(uint64_t *qword)
{
    unsigned char *buffer = (unsigned char *)malloc(sizeof(uint64_t));
    if (buffer == NULL) {
        throwErrorException("malloc failed");
    }

    size_t ret = read(buffer, sizeof(uint64_t));

    if (ret != sizeof(uint64_t)) {
        free(buffer);
        throwErrorException("read failed");
    } else {
        *qword = (uint64_t)buffer[0] << 56 |
                 (uint64_t)buffer[1] << 48 |
                 (uint64_t)buffer[2] << 40 |
                 (uint64_t)buffer[3] << 32 |
                 (uint64_t)buffer[4] << 24 |
                 (uint64_t)buffer[5] << 16 |
                 (uint64_t)buffer[6] <<  8 |
                 (uint64_t)buffer[7];
        free(buffer);
    }

    return ret;
}

size_t cq::File::writeByte(const unsigned char &byte)
{
    if (fwrite(&byte, 1, 1, m_fp) != 1) {
        throwErrorException("fwrite failed");
    }
    return 1;
}

size_t cq::File::writeUint8(const uint8_t &byte)
{
    return writeByte(byte);
}

size_t cq::File::writeUint16(const uint16_t &word)
{
    size_t ret = 0;
    ret += writeByte((unsigned char)(word >> 8) & 0xFF);
    ret += writeByte((unsigned char)(word     ) & 0xFF);
    return ret;
}

size_t cq::File::writeUint32(const uint32_t &dword)
{
    size_t ret = 0;
    ret += writeByte((unsigned char)(dword >> 24) & 0xFF);
    ret += writeByte((unsigned char)(dword >> 16) & 0xFF);
    ret += writeByte((unsigned char)(dword >>  8) & 0xFF);
    ret += writeByte((unsigned char)(dword      ) & 0xFF);
    return ret;
}

size_t cq::File::writeUint64(const uint64_t &qword)
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


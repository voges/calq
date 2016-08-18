/** @file File.h
 *  @brief This file contains the definitions of the File class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef FILE_H
#define FILE_H

#include <string>

class File {
public:
    File(void);
    File(const std::string &path, const char *mode);
    ~File(void);

    void open(const std::string &path, const char *mode);
    void close(void);

    void advance(const size_t &offset);
    bool eof(void);
    void * handle(void);
    void seek(const size_t &pos);
    size_t size(void);
    size_t tell(void);

    size_t read(void *buffer, const size_t &size);
    size_t write(void *buffer, const size_t &size);

    size_t readByte(unsigned char *byte);
    size_t readUint8(uint8_t *byte);
    size_t readUint16(uint16_t *word);
    size_t readUint32(uint32_t *dword);
    size_t readUint64(uint64_t *qword);

    size_t writeByte(const unsigned char &byte);
    size_t writeUint8(const uint8_t &byte);
    size_t writeUint16(const uint16_t &word);
    size_t writeUint32(const uint32_t &dword);
    size_t writeUint64(const uint64_t &qword);

private:
    FILE *fp;
    size_t fsize;
};

#endif // FILE_H


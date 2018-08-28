/** @file File.h
 *  @brief This file contains the definition of the File class.
 */

// Copyright 2015-2017 Leibniz Universitaet Hannover

#ifndef CALQ_IO_FILE_H_
#define CALQ_IO_FILE_H_

#include <string>

namespace calq {

class File {
 public:
    enum Mode {
        MODE_READ = 0, MODE_WRITE = 1
    };

    File();
    File(const std::string &path, Mode mode);
    virtual ~File();

    void open(const std::string &path, Mode mode);
    void close();

    void advance(size_t offset);
    bool eof() const;
    void* handle() const;
    void seek(size_t pos);
    size_t size() const;
    size_t tell() const;

    size_t nrReadBytes() const;
    size_t nrWrittenBytes() const;

    bool isReadable() const;
    bool isWritable() const;

    size_t read(void* buffer, size_t size);
    size_t write(const void* buffer, size_t size);

    size_t readByte(unsigned char* byte);
    size_t readUint8(uint8_t* byte);
    size_t readUint16(uint16_t* word);
    size_t readUint32(uint32_t* dword);
    size_t readUint64(uint64_t* qword);

    size_t writeByte(unsigned char byte);
    size_t writeUint8(uint8_t byte);
    size_t writeUint16(uint16_t word);
    size_t writeUint32(uint32_t dword);
    size_t writeUint64(uint64_t qword);

 protected:
    FILE* fp_;
    size_t fsize_;
    bool isOpen_;
    Mode mode_;
    size_t nrReadBytes_;
    size_t nrWrittenBytes_;
};

}  // namespace calq

#endif  // CALQ_IO_FILE_H_


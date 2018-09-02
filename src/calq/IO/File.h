/** @file File.h
 *  @brief This file contains the definition of the File class.
 */

// Copyright 2015-2017 Leibniz Universitaet Hannover

#ifndef CALQ_IO_FILE_H_
#define CALQ_IO_FILE_H_

#include <fstream>
#include <string>
#include <Common/ErrorExceptionReporter.h>

namespace calq {

class File {
 public:
    enum class Mode {
        MODE_READ = 0, MODE_WRITE = 1
    };

    File();
    File(const std::string &path, Mode mode);
    virtual ~File();

    void open(const std::string &path, Mode mode);
    void close();

    void advance(size_t offset);
    bool eof() const;
    void seek(size_t pos);
    size_t size() const;
    size_t tell();

    size_t nrReadBytes() const;
    size_t nrWrittenBytes() const;

    bool readLine(char* s, std::streamsize n);

    bool isReadable() const;
    bool isWritable() const;

    template<typename T>
    size_t readValue(T* dword, size_t number = 1) {
        size_t ret = sizeof(T) * number;
        try {
            filestream.read(reinterpret_cast<char*>(dword), sizeof(T) * number);
        } catch (std::exception &e) {
            throwErrorException(std::string("Read failed: ") + e.what());
        }
        nrReadBytes_ += ret;
        return ret;
    }

    template<typename T>
    size_t writeValue(const T* dword, size_t number = 1) {
        size_t ret = sizeof(T) * number;
        try {
            filestream.write(reinterpret_cast<const char*>(dword), sizeof(T) * number);
        } catch (std::exception &e) {
            throwErrorException(std::string("Write failed: ") + e.what());
        }
        nrWrittenBytes_ += ret;
        return ret;
    }

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
    size_t fsize_;
    Mode mode_;
    size_t nrReadBytes_;
    size_t nrWrittenBytes_;
    std::fstream filestream;
};

}  // namespace calq

#endif  // CALQ_IO_FILE_H_


/** @file File.cc
 *  @brief This file contains the implementation of the File class.
 */

// Copyright 2015-2017 Leibniz Universitaet Hannover

#include "IO/File.h"

#include <climits>

#include "Common/ErrorExceptionReporter.h"
#include "File.h"


namespace calq {

File::File() : fsize_(0), mode_(File::Mode::MODE_READ), nrReadBytes_(0), nrWrittenBytes_(0) {
    filestream.exceptions(std::ifstream::failbit | std::ifstream::badbit);
}

File::File(const std::string &path, Mode mode) : fsize_(0), mode_(mode), nrReadBytes_(0), nrWrittenBytes_(0) {
    filestream.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    open(path, mode);
}

File::~File() {
    close();
}

void File::open(const std::string &path, Mode mode) {
    if (path.empty()) {
        throwErrorException("path is empty");
    }

    mode_ = mode;

    auto m = (mode_ == Mode::MODE_READ) ? std::ifstream::in : std::ifstream::out;
    try {
        filestream.open(path, m);
    } catch (std::exception &e) {
        throwErrorException(std::string("Error opening file: ") + e.what());
    }

    try {
        filestream.seekg(0, std::ifstream::end);
        fsize_ = static_cast<size_t>(filestream.tellg());
        filestream.seekg(0, std::ifstream::beg);
    } catch (std::exception &e) {
        throwErrorException(std::string("Error obtaining file size: ") + e.what());
    }
}

void File::close() {
    if (filestream.is_open()) {
        try {
            filestream.close();
        } catch (std::exception &e) {
            throwErrorException(std::string("Failed to close file: ") + e.what());
        }
    }
}

void File::advance(size_t offset) {
    try {
        filestream.seekg(offset, std::ios_base::seekdir::_S_cur);
    } catch (std::exception &e) {
        throwErrorException(std::string("Seek failed: ") + e.what());
    }
}

bool File::eof() const {
    return filestream.eof();
}

void File::seek(size_t pos) {
    if (pos > LONG_MAX) {
        throwErrorException("pos out of range");
    }
    try {
        filestream.seekg(pos);
    } catch (std::exception &e) {
        throwErrorException(std::string("Seek failed: ") + e.what());
    }
}

size_t File::size() const {
    return fsize_;
}

size_t File::tell() {
    try {
        return static_cast<size_t>(filestream.tellg());
    } catch (std::exception &e) {
        if(!eof())
            throwErrorException(std::string("Tell failed: ") + e.what());
        filestream.clear(std::ios_base::goodbit | std::ios_base::eofbit);
    }
    return 0;
}

size_t File::nrReadBytes() const {
    if (mode_ != Mode::MODE_READ) {
        throwErrorException("File is not open in read mode");
    }
    return nrReadBytes_;
}

size_t File::nrWrittenBytes() const {
    if (mode_ != Mode::MODE_WRITE) {
        throwErrorException("File is not open in write mode");
    }
    return nrWrittenBytes_;
}

bool File::isReadable() const {
    return filestream.is_open() && mode_ == Mode::MODE_READ;
}

bool File::isWritable() const {
    return filestream.is_open() && mode_ == Mode::MODE_WRITE;
}

size_t File::read(void* buffer, size_t size) {
    return readValue(reinterpret_cast<unsigned char*> (buffer), size);
}

size_t File::write(const void* buffer, size_t size) {
    return writeValue(reinterpret_cast<const unsigned char*> (buffer), size);
}

size_t File::readByte(unsigned char* byte) {
    return readValue(byte, 1);
}

size_t File::readUint8(uint8_t* byte) {
    return readValue(byte, 1);
}

size_t File::readUint16(uint16_t* word) {
    return readValue(word, 1);
}

size_t File::readUint32(uint32_t* dword) {
    return readValue(dword, 1);
}

size_t File::readUint64(uint64_t* qword) {
    return readValue(qword, 1);
}

size_t File::writeByte(unsigned char byte) {
    return writeValue(&byte, 1);
}

size_t File::writeUint8(uint8_t byte) {
    return writeValue(&byte, 1);
}

size_t File::writeUint16(uint16_t word) {
    return writeValue(&word, 1);
}

size_t File::writeUint32(uint32_t dword) {
    return writeValue(&dword, 1);
}

size_t File::writeUint64(uint64_t qword) {
    return writeValue(&qword, 1);
}

bool File::readLine(char* s, std::streamsize n) {
    try {
        filestream.getline(s, n);
    } catch (std::exception &e) {
        if(!eof())
            throwErrorException(std::string("readLine failed: ") + e.what());
        filestream.clear(std::ios_base::goodbit | std::ios_base::eofbit);
        return false;
    }
    return !eof();
}

}  // namespace calq


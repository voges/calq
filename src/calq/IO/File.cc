/** @file File.cc
 *  @brief This file contains the implementation of the File class.
 */

// Copyright 2015-2017 Leibniz Universitaet Hannover

#include "IO/File.h"

#include <climits>

#include "Common/ErrorExceptionReporter.h"

namespace calq {

File::File() : fp_(nullptr), fsize_(0), isOpen_(false), mode_(File::MODE_READ), nrReadBytes_(0), nrWrittenBytes_(0) {
}

File::File(const std::string &path, Mode mode) : fp_(nullptr), fsize_(0), isOpen_(false), mode_(mode), nrReadBytes_(0), nrWrittenBytes_(0) {
    if (path.empty()) {
        throwErrorException("path is empty");
    }

    open(path, mode);
}

File::~File() {
    close();
}

void File::open(const std::string &path, Mode mode) {
    if (path.empty()) {
        throwErrorException("path is empty");
    }
    if (fp_ != nullptr) {
        throwErrorException("File pointer already in use");
    }

    const char* m = nullptr;
    if (mode == MODE_READ) {
        m = "rb";
        mode_ = mode;
    } else if (mode == MODE_WRITE) {
        m = "wb";
        mode_ = mode;
    } else {
        throwErrorException("Unknown mode");
    }

#ifdef CQ_OS_WINDOWS
    int err = fopen_s(&fp_, path.c_str(), m);
    if (err != 0) {
        throwErrorException("Failed to open file");
    }
#else
    fp_ = fopen(path.c_str(), m);
    if (fp_ == nullptr) {
        throwErrorException("Failed to open file");
    }
#endif

    // Compute file size
    fseek(fp_, 0, SEEK_END);
    fsize_ = static_cast<size_t>(ftell(fp_));
    fseek(fp_, 0, SEEK_SET);

    isOpen_ = true;
}

void File::close() {
    if (isOpen_) {
        if (fp_ != nullptr) {
            fclose(fp_);
            fp_ = nullptr;
        } else {
            throwErrorException("Failed to close file");
        }
    }
}

void File::advance(size_t offset) {
    int ret = fseek(fp_, (long int) offset, SEEK_CUR);
    if (ret != 0) {
        throwErrorException("fseek failed");
    }
}

bool File::eof() const {
    int eof = feof(fp_);
    return eof != 0;
}

void* File::handle() const {
    return fp_;
}

void File::seek(size_t pos) {
    if (pos > LONG_MAX) {
        throwErrorException("pos out of range");
    }
    int ret = fseek(fp_, (long) pos, SEEK_SET);
    if (ret != 0) {
        throwErrorException("fseek failed");
    }
}

size_t File::size() const {
    return fsize_;
}

size_t File::tell() const {
    long int offset = ftell(fp_);
    if (offset == -1) {
        throwErrorException("ftell failed");
    }
    return static_cast<size_t>(offset);
}

size_t File::nrReadBytes() const {
    if (mode_ != MODE_READ) {
        throwErrorException("File is not open in read mode");
    }
    return nrReadBytes_;
}

size_t File::nrWrittenBytes() const {
    if (mode_ != MODE_WRITE) {
        throwErrorException("File is not open in write mode");
    }
    return nrWrittenBytes_;
}

bool File::isReadable() const {
    return isOpen_ && mode_ == MODE_READ;
}

bool File::isWritable() const {
    return isOpen_ && mode_ == MODE_WRITE;
}

size_t File::read(void* buffer, size_t size) {
    if (buffer == nullptr) {
        throwErrorException("buffer is nullptr");
    }
    if (size == 0) {
        return 0;
    }
    size_t ret = fread(buffer, 1, size, fp_);
    if (ret != size) {
        throwErrorException("fread failed");
    }
    nrReadBytes_ += ret;
    return ret;
}

size_t File::write(const void* buffer, size_t size) {
    if (buffer == nullptr) {
        throwErrorException("buffer is nullptr");
    }
    if (size == 0) {
        return 0;
    }
    size_t ret = fwrite(buffer, 1, size, fp_);
    if (ret != size) {
        throwErrorException("fwrite failed");
    }
    nrWrittenBytes_ += ret;
    return ret;
}

size_t File::readByte(unsigned char* byte) {
    size_t ret = fread(byte, 1, 1, fp_);
    if (ret != sizeof(unsigned char)) {
        throwErrorException("fread failed");
    }
    nrReadBytes_++;
    return ret;
}

size_t File::readUint8(uint8_t* byte) {
    return readByte(byte);
}

size_t File::readUint16(uint16_t* word) {
    auto* buffer = (unsigned char*) malloc(sizeof(uint16_t));
    if (buffer == nullptr) {
        throwErrorException("malloc failed");
    }

    size_t ret = read(buffer, sizeof(uint16_t));

    if (ret != sizeof(uint16_t)) {
        free(buffer);
        throwErrorException("read failed");
    } else {
        *word = (uint16_t) buffer[2] << 8 | (uint16_t) buffer[3]; //NOLINT
        free(buffer);
    }

    return ret;
}

size_t File::readUint32(uint32_t* dword) {
    auto* buffer = (unsigned char*) malloc(sizeof(uint32_t));
    if (buffer == nullptr) {
        throwErrorException("malloc failed");
    }

    size_t ret = read(buffer, sizeof(uint32_t));

    if (ret != sizeof(uint32_t)) {
        free(buffer);
        throwErrorException("read failed");
    } else {
        *dword = (uint32_t) buffer[0] << 24 | //NOLINT
                 (uint32_t) buffer[1] << 16 | //NOLINT
                 (uint32_t) buffer[2] << 8 | //NOLINT
                 (uint32_t) buffer[3];        //NOLINT
        free(buffer);
    }

    return ret;
}

size_t File::readUint64(uint64_t* qword) {
    auto* buffer = (unsigned char*) malloc(sizeof(uint64_t));
    if (buffer == nullptr) {
        throwErrorException("malloc failed");
    }

    size_t ret = read(buffer, sizeof(uint64_t));

    if (ret != sizeof(uint64_t)) {
        free(buffer);
        throwErrorException("read failed");
    } else {
        *qword = (uint64_t) buffer[0] << 56 | //NOLINT
                 (uint64_t) buffer[1] << 48 | //NOLINT
                 (uint64_t) buffer[2] << 40 | //NOLINT
                 (uint64_t) buffer[3] << 32 | //NOLINT
                 (uint64_t) buffer[4] << 24 | //NOLINT
                 (uint64_t) buffer[5] << 16 | //NOLINT
                 (uint64_t) buffer[6] << 8 | //NOLINT
                 (uint64_t) buffer[7];        //NOLINT
        free(buffer);
    }

    return ret;
}

size_t File::writeByte(unsigned char byte) {
    size_t ret = fwrite(&byte, 1, 1, fp_);
    if (ret != sizeof(unsigned char)) {
        throwErrorException("fwrite failed");
    }
    nrWrittenBytes_++;
    return ret;
}

size_t File::writeUint8(uint8_t byte) {
    return writeByte(byte);
}

size_t File::writeUint16(uint16_t word) {
    size_t ret = 0;
    ret += writeByte((unsigned char) (word >> 8) & 0xFF); //NOLINT
    ret += writeByte((unsigned char) (word) & 0xFF); //NOLINT
    return ret;
}

size_t File::writeUint32(uint32_t dword) {
    size_t ret = 0;
    ret += writeByte((unsigned char) (dword >> 24) & 0xFF); //NOLINT
    ret += writeByte((unsigned char) (dword >> 16) & 0xFF); //NOLINT
    ret += writeByte((unsigned char) (dword >> 8) & 0xFF); //NOLINT
    ret += writeByte((unsigned char) (dword) & 0xFF); //NOLINT
    return ret;
}

size_t File::writeUint64(uint64_t qword) {
    size_t ret = 0;
    ret += writeByte((unsigned char) (qword >> 56) & 0xFF); //NOLINT
    ret += writeByte((unsigned char) (qword >> 48) & 0xFF); //NOLINT
    ret += writeByte((unsigned char) (qword >> 40) & 0xFF); //NOLINT
    ret += writeByte((unsigned char) (qword >> 32) & 0xFF); //NOLINT
    ret += writeByte((unsigned char) (qword >> 24) & 0xFF); //NOLINT
    ret += writeByte((unsigned char) (qword >> 16) & 0xFF); //NOLINT
    ret += writeByte((unsigned char) (qword >> 8) & 0xFF); //NOLINT
    ret += writeByte((unsigned char) (qword) & 0xFF); //NOLINT
    return ret;
}

}  // namespace calq


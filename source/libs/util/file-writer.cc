#include "file-writer.h"
#include <climits>

namespace cip {

FileWriter::FileWriter() : fsize_(0), mode_(FileWriter::Mode::MODE_READ), nrReadBytes_(0), nrWrittenBytes_(0) {
    filestream.exceptions(std::ifstream::failbit | std::ifstream::badbit);
}

FileWriter::FileWriter(const std::string &path, Mode mode) : fsize_(0), mode_(mode), nrReadBytes_(0), nrWrittenBytes_(0) {
    filestream.exceptions(std::ifstream::failbit | std::ifstream::badbit);
    this->open(path, mode);
}

FileWriter::~FileWriter() { close(); }

void FileWriter::open(const std::string &path, Mode mode) {
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

void FileWriter::close() {
    if (filestream.is_open()) {
        try {
            filestream.close();
        } catch (std::exception &e) {
            throwErrorException(std::string("Failed to close file: ") + e.what());
        }
    }
}

bool FileWriter::eof() const { return filestream.eof(); }

void FileWriter::seek(size_t pos) {
    if (pos > LONG_MAX) {
        throwErrorException("pos out of range");
    }
    try {
        filestream.seekg(pos);
    } catch (std::exception &e) {
        throwErrorException(std::string("Seek failed: ") + e.what());
    }
}

size_t FileWriter::size() const { return fsize_; }

size_t FileWriter::tell() {
    try {
        return static_cast<size_t>(filestream.tellg());
    } catch (std::exception &e) {
        if (!eof()) throwErrorException(std::string("Tell failed: ") + e.what());
        filestream.clear(std::ios_base::goodbit | std::ios_base::eofbit);
    }
    return fsize_;
}

size_t FileWriter::read(void *buffer, size_t size) { return readValue(reinterpret_cast<unsigned char *>(buffer), size); }

size_t FileWriter::write(const void *buffer, size_t size) {
    return writeValue(reinterpret_cast<const unsigned char *>(buffer), size);
}

size_t FileWriter::readUint8(uint8_t *byte) { return readValue(byte, 1); }

size_t FileWriter::readUint32(uint32_t *dword) { return readValue(dword, 1); }

size_t FileWriter::readUint64(uint64_t *qword) { return readValue(qword, 1); }

size_t FileWriter::writeUint8(uint8_t byte) { return writeValue(&byte, 1); }

size_t FileWriter::writeUint32(uint32_t dword) { return writeValue(&dword, 1); }

size_t FileWriter::writeUint64(uint64_t qword) { return writeValue(&qword, 1); }

bool FileWriter::readLine(char *s, std::streamsize n) {
    try {
        filestream.getline(s, n);
    } catch (std::exception &e) {
        if (!eof()) {
            throwErrorException(std::string("readLine failed: ") + e.what());
        }
        filestream.clear(std::ios_base::goodbit | std::ios_base::eofbit);
        return false;
    }
    return !eof();
}

}  // namespace cip

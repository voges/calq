/**
 * @file file-writer.cc
 */

#include "file-writer.h"

namespace calq {

FileWriter::FileWriter(const std::string &path) : ofs_() {
    ofs_.open(path, std::ifstream::out | std::ifstream::binary);
    if (!ofs_.is_open()) {
        throw ErrorException("failed to open file: " + path);
    }
}

FileWriter::~FileWriter() { close_(); }

void FileWriter::close() { close_(); }

size_t FileWriter::write(void *const buffer, const size_t size) {
    return write_(reinterpret_cast<const unsigned char *>(buffer), size);
}

size_t FileWriter::writeUint8(const uint8_t byte) { return write_(&byte); }

size_t FileWriter::writeUint16(const uint16_t word) { return write_(&word); }

size_t FileWriter::writeUint32(const uint32_t dword) { return write_(&dword); }

size_t FileWriter::writeUint64(const uint64_t qword) { return write_(&qword); }

void FileWriter::close_() { ofs_.close(); }

}  // namespace calq

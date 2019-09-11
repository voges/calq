/**
 * @file file-writer.cc
 */

#include "file-writer.h"

namespace calq {

FileWriter::FileWriter(const std::string &path) : File() { open(path, Mode::WRITE); }

size_t FileWriter::write(void *const buffer, const size_t size) {
    return writeValue(reinterpret_cast<const unsigned char *>(buffer), size);
}

size_t FileWriter::writeUint8(const uint8_t byte) { return writeValue(&byte); }

size_t FileWriter::writeUint16(const uint16_t word) { return writeValue(&word); }

size_t FileWriter::writeUint32(const uint32_t dword) { return writeValue(&dword); }

size_t FileWriter::writeUint64(const uint64_t qword) { return writeValue(&qword); }

}  // namespace calq

/**
 * @file file-reader.cc
 */

#include "file-reader.h"

namespace calq {

FileReader::FileReader(const std::string &path) : File() { open(path, Mode::READ); }

size_t FileReader::read(void *const buffer, const size_t size) {
    return readValue(reinterpret_cast<unsigned char *>(buffer), size);
}

size_t FileReader::readUint8(uint8_t *const byte) { return readValue(byte); }

size_t FileReader::readUint16(uint16_t *const word) { return readValue(word); }

size_t FileReader::readUint32(uint32_t *const dword) { return readValue(dword); }

size_t FileReader::readUint64(uint64_t *const qword) { return readValue(qword); }

}  // namespace calq

#include "file-reader.h"

namespace calq {

FileReader::FileReader(const std::string &path) : ifs_() {
    ifs_.open(path, std::ifstream::in | std::ifstream::binary);
    if (!ifs_.is_open()) {
        throw ErrorException("failed to open file: " + path);
    }
}

FileReader::~FileReader() { close_(); }

void FileReader::close() { close_(); }

size_t FileReader::read(void *const buffer, const size_t size) {
    return read_(reinterpret_cast<unsigned char *>(buffer), size);
}

size_t FileReader::readUint8(uint8_t *const byte) { return read_(byte); }

size_t FileReader::readUint16(uint16_t *const word) { return read_(word); }

size_t FileReader::readUint32(uint32_t *const dword) { return read_(dword); }

size_t FileReader::readUint64(uint64_t *const qword) { return read_(qword); }

void FileReader::close_() { ifs_.close(); }

}  // namespace calq

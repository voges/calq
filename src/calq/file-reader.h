/**
 * @file file-reader.h
 */

#ifndef CALQ_FILE_READER_H_
#define CALQ_FILE_READER_H_

#include <cassert>
#include <fstream>
#include "exceptions.h"

namespace calq {

class FileReader {
   public:
    FileReader() = delete;
    explicit FileReader(const std::string &path);
    FileReader(const FileReader &) = delete;
    FileReader &operator=(const FileReader &) = delete;
    FileReader(FileReader &&) = delete;
    FileReader &operator=(FileReader &&) = delete;
    ~FileReader();

    void close();
    size_t read(void *buffer, size_t size = 1);
    size_t readUint8(uint8_t *byte);
    size_t readUint16(uint16_t *word);
    size_t readUint32(uint32_t *dword);
    size_t readUint64(uint64_t *qword);

   private:
    void close_();

    template <typename T>
    size_t read_(T *const buffer, const size_t n = 1) {
        assert(buffer != nullptr);

        // Read from file
        ifs_.read(reinterpret_cast<char *>(buffer), sizeof(T) * n);

        // Check whether we read exactly n items
        if (!ifs_.good()) {
            if (ifs_.eof() && ifs_.fail() && !ifs_.bad()) {
                // Everything okay, we just reached the EOF; only ifs_.gcount() bytes were read
                return ifs_.gcount();
            } else {
                throw ErrorException("failed to read from file");
            }
        }

        // Return the number of read bytes
        return sizeof(T) * n;
    }

    std::ifstream ifs_;
};

}  // namespace calq

#endif  // CALQ_FILE_READER_H_

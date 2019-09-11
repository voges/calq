/**
 * @file file-reader.h
 */

#ifndef CALQ_FILE_READER_H_
#define CALQ_FILE_READER_H_

#include <cassert>
#include "errors.h"
#include "file.h"

namespace calq {

class FileReader : public File {
   public:
    explicit FileReader(const std::string &path);

    template <typename T>
    size_t readValue(T *const value, const size_t n = 1) {
        assert(value != nullptr);

        size_t ret = fread(value, sizeof(T), n, fp_);

        // Check whether we read exactly n items
        if (ret != n) {
            if (eof()) {
                // Everything okay, we just reached the EOF
                return 0;
            } else if (error()) {
                throwErrorException("Error occurred while trying to read from file");
            } else {
                throwErrorException("Failed to read from file");
            }
        }

        // Return the number of read bytes
        return sizeof(T) * n;
    }

    size_t read(void *buffer, size_t size);
    size_t readUint8(uint8_t *byte);
    size_t readUint16(uint16_t *word);
    size_t readUint32(uint32_t *dword);
    size_t readUint64(uint64_t *qword);
};

}  // namespace calq

#endif  // CALQ_FILE_READER_H_

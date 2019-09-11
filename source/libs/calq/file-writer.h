/**
 * @file file-writer.h
 */

#ifndef CALQ_FILE_WRITER_H_
#define CALQ_FILE_WRITER_H_

#include "errors.h"
#include "file.h"

namespace calq {

class FileWriter : public File {
   public:
    explicit FileWriter(const std::string &path);

    template <typename T>
    size_t writeValue(const T *const value, const size_t n = 1) {
        size_t ret = fwrite(value, sizeof(T), n, fp_);

        // Check whether the write was successful
        if (ret != n) {
            throwErrorException("Failed to write to file");
        }

        // Return the number of written bytes
        return sizeof(T) * n;
    }

    size_t write(void *buffer, size_t size);
    size_t writeUint8(uint8_t byte);
    size_t writeUint16(uint16_t word);
    size_t writeUint32(uint32_t dword);
    size_t writeUint64(uint64_t qword);
};

}  // namespace calq

#endif  // CALQ_FILE_WRITER_H_

/**
 * @file file-writer.h
 */

#ifndef CALQ_FILE_WRITER_H_
#define CALQ_FILE_WRITER_H_

#include <cassert>
#include <fstream>
#include "exceptions.h"

namespace calq {

class FileWriter {
   public:
    FileWriter() = delete;
    explicit FileWriter(const std::string &path);
    FileWriter(const FileWriter &) = delete;
    FileWriter &operator=(const FileWriter &) = delete;
    FileWriter(FileWriter &&) = delete;
    FileWriter &operator=(FileWriter &&) = delete;
    ~FileWriter();

    void close();
    size_t write(void *buffer, size_t size = 1);
    size_t writeUint8(uint8_t byte);
    size_t writeUint16(uint16_t word);
    size_t writeUint32(uint32_t dword);
    size_t writeUint64(uint64_t qword);

   private:
    void close_();

    template <typename T>
    size_t write_(const T *const buffer, const size_t n = 1) {
        assert(buffer != nullptr);

        // Write to file
        ofs_.write(reinterpret_cast<const char *>(buffer), sizeof(T) * n);

        // Check whether the write was successful
        if (!ofs_.good()) {
            throwErrorException("Failed to write to file");
        }

        // Return the number of written bytes
        return sizeof(T) * n;
    }

    std::ofstream ofs_;
};

}  // namespace calq

#endif  // CALQ_FILE_WRITER_H_

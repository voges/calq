#ifndef CALQ_FILE_WRITER_H_
#define CALQ_FILE_WRITER_H_

#include <fstream>
#include <string>
// #include "cip/errors.h"

namespace calq {

class FileWriter {
   public:
    FileWriter();
    FileWriter(const std::string &path);
    virtual ~FileWriter();

    void open(const std::string &path);
    void close();

    bool eof() const;
    void seek(size_t pos);
    size_t size() const;
    size_t tell();

    bool readLine(char *s, std::streamsize n);

    template <typename T>
    size_t readValue(T *dword, size_t number = 1) {
        size_t ret = sizeof(T) * number;
        try {
            filestream.read(reinterpret_cast<char *>(dword), sizeof(T) * number);
        } catch (std::exception &e) {
            // throwErrorException(std::string("Read failed: ") + e.what());
        }
        nrReadBytes_ += ret;
        return ret;
    }

    template <typename T>
    size_t writeValue(const T *dword, size_t number = 1) {
        size_t ret = sizeof(T) * number;
        try {
            filestream.write(reinterpret_cast<const char *>(dword), sizeof(T) * number);
        } catch (std::exception &e) {
            // throwErrorException(std::string("Write failed: ") + e.what());
        }
        nrWrittenBytes_ += ret;
        return ret;
    }

    size_t read(void *buffer, size_t size);
    size_t write(const void *buffer, size_t size);

    size_t readUint8(uint8_t *byte);
    size_t readUint32(uint32_t *dword);
    size_t readUint64(uint64_t *qword);

    size_t writeUint8(uint8_t byte);
    size_t writeUint32(uint32_t dword);
    size_t writeUint64(uint64_t qword);

   protected:
    size_t fsize_;
    size_t nrReadBytes_;
    size_t nrWrittenBytes_;
    std::fstream filestream;
};

}  // namespace calq

#endif  // CALQ_FILE_WRITER_H_

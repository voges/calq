#ifndef UTIL_FILE_READER_H_
#define UTIL_FILE_READER_H_

#include <string>
#include "constants.h"

namespace util {

class FileReader {
   public:
    explicit FileReader(const std::string &path);
    virtual ~FileReader();
    void advance(int64_t offset);
    bool eof() const;
    void *handle() const;
    void readLine(std::string *line);
    void seekFromCur(int64_t offset);
    void seekFromEnd(int64_t offset);
    void seekFromSet(int64_t offset);
    size_t size() const;
    int64_t tell() const;

   protected:
    FILE *fp_;
    size_t size_;

   private:
    void close();
    void open(const std::string &path);
    void seek(int64_t offset, int whence);

    static const size_t MAX_LINE_LENGTH = sizeof(char) * (4 * KB);
    char *line_;
};

}  // namespace util

#endif  // UTIL_FILE_READER_H_

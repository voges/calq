/**
 * @file file.h
 */

#ifndef CALQ_FILE_H_
#define CALQ_FILE_H_

#include <fstream>
#include <string>

namespace calq {

class File {
   public:
    enum class Mode { READ, WRITE };

    explicit File();
    virtual ~File();
    void advance(int64_t offset);
    void close();
    bool eof() const;
    bool error() const;
    void open(const std::string &path, Mode mode);
    void seekFromCur(int64_t offset);
    void seekFromEnd(int64_t offset);
    void seekFromSet(int64_t offset);
    size_t size() const;
    int64_t tell() const;

   protected:
    FILE *fp_;
    size_t size_;

   private:
    void seek(int64_t offset, int whence);
};

}  // namespace calq

#endif  // CALQ_FILE_H_

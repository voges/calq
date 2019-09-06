#ifndef UTIL_FILE_READER_H_
#define UTIL_FILE_READER_H_

#include <string>
#include "constants.h"
#include "file.h"

namespace util {

class FileReader : public File {
   public:
    explicit FileReader(const std::string &path);
    void readLine(std::string *line);

   private:
    static const size_t MAX_LINE_LENGTH = sizeof(char) * (4 * KB);
    char *line_;
};

}  // namespace util

#endif  // UTIL_FILE_READER_H_

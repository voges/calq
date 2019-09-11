/**
 * @file file-line-reader.h
 */

#ifndef CALQ_FILE_LINE_READER_H_
#define CALQ_FILE_LINE_READER_H_

#include <string>
#include "constants.h"
#include "file.h"

namespace calq {

class FileLineReader : public File {
   public:
    explicit FileLineReader(const std::string &path);
    ~FileLineReader() override;
    void readLine(std::string *line);

   private:
    static const size_t MAX_LINE_LENGTH = sizeof(char) * (4 * constants::KB);
    char *line_;
};

}  // namespace calq

#endif  // CALQ_FILE_LINE_READER_H_

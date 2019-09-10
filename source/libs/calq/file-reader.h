#ifndef CALQ_FILE_READER_H_
#define CALQ_FILE_READER_H_

#include <string>
#include "constants.h"
#include "file.h"

namespace calq {

class FileReader : public File {
   public:
    explicit FileReader(const std::string &path);
    virtual ~FileReader();
    void readLine(std::string *line);

   private:
    static const size_t MAX_LINE_LENGTH = sizeof(char) * (4 * KB);
    char *line_;
};

}  // namespace calq

#endif  // CALQ_FILE_READER_H_

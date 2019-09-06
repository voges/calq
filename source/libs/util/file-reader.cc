#include "file-reader.h"
#include <climits>
#include <cstring>
#include <stdexcept>

namespace util {

FileReader::FileReader(const std::string &path) : File(path) {}

void FileReader::readLine(std::string *const line) {
    line->clear();

    char *rc = fgets(line_, MAX_LINE_LENGTH, fp_);

    if (rc == nullptr) {
        // Error during read or EOF. Nothing has been read. We just return to the caller and leave 'line' empty.
        return;
    }

    if (eof()) {
        // EOF was reached but contents were read into 'line_'. We proceed processing the read contents.
    }

    // Trim line
    size_t l = strlen(line_) - 1;
    while (l && ((line_[l] == '\r') || (line_[l] == '\n'))) {
        line_[l--] = '\0';
    }

    *line = line_;
}

}  // namespace util

/**
 * @file file-line-reader.cc
 */

#include "file-line-reader.h"
#include <cassert>
#include "errors.h"
#include "string-helpers.h"

namespace calq {

FileLineReader::FileLineReader(const std::string &path) : File(), line_(nullptr) {
    open(path, Mode::READ);

    line_ = reinterpret_cast<char *>(malloc(MAX_LINE_LENGTH));
    if (line_ == nullptr) {
        throwErrorException("Failed to allocate memory");
    }
}

FileLineReader::~FileLineReader() { free(line_); }

void FileLineReader::readLine(std::string *const line) {
    assert(line != nullptr);

    line->clear();

    char *rc = fgets(line_, MAX_LINE_LENGTH, fp_);

    if (rc == nullptr) {
        // Error during read or EOF. Nothing has been read. We just return to the caller and leave 'line' empty.
        return;
    }

    if (eof()) {
        // EOF was reached but contents were read into 'line_'. We proceed processing the contents.
    }

    *line = line_;
    *line = rtrim(*line);
}

}  // namespace calq

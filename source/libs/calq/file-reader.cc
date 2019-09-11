/**
 * @file file-reader.cc
 */

#include "file-reader.h"
#include <cassert>
#include "errors.h"
#include "string-helpers.h"

namespace calq {

FileReader::FileReader(const std::string &path) : File(path), line_(nullptr) {
    line_ = reinterpret_cast<char *>(malloc(MAX_LINE_LENGTH));
    if (line_ == nullptr) {
        throwErrorException("Failed to allocate memory");
    }
}

FileReader::~FileReader() { free(line_); }

void FileReader::readLine(std::string *const line) {
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

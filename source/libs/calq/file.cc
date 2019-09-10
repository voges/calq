/**
 * @file file.cc
 */

#include "file.h"
#include <cassert>
#include <climits>
#include <stdexcept>
#include "errors.h"

namespace calq {

File::File(const std::string &path) : fp_(nullptr), size_(0) {
    open(path);

    // Compute file size
    seekFromEnd(0);
    size_ = tell();
    seekFromSet(0);
}

File::~File() { close(); }

void File::advance(const int64_t offset) { seekFromCur(offset); }

bool File::eof() const { return feof(fp_) != 0; }

void File::seekFromCur(const int64_t offset) { seek(offset, SEEK_CUR); }

void File::seekFromEnd(const int64_t offset) { seek(offset, SEEK_END); }

void File::seekFromSet(const int64_t offset) { seek(offset, SEEK_SET); }

size_t File::size() const { return size_; }

int64_t File::tell() const {
    auto offset = static_cast<int64_t>(ftell(fp_));

    if (offset == -1) {
        throwErrorException("Failed to obtain the current value of the file position indicator");
    }

    return offset;
}

void File::close() {
    if (fp_ != nullptr) {
        fclose(fp_);
        fp_ = nullptr;
    } else {
        throwErrorException("Failed to close file");
    }
}

void File::open(const std::string &path) {
    assert(fp_ == nullptr);

    const char *mode = "rb";

#ifdef _WIN32
    int rc = fopen_s(&fp_, path.c_str(), mode);
    if (rc != 0) {
        throwErrorException("Failed to open file: " + path);
    }
#else
    fp_ = fopen(path.c_str(), mode);
    if (fp_ == nullptr) {
        throwErrorException("Failed to open file: " + path);
    }
#endif
}

void File::seek(const int64_t offset, const int whence) {
    if (offset < LONG_MIN || offset > LONG_MAX) {
        throwErrorException("New value for the file position indicator is out of range");
    }

    int rc = fseek(fp_, offset, whence);
    if (rc != 0) {
        throwErrorException("Failed to set the file position indicator");
    }
}

}  // namespace calq

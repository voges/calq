/**
 * @file file.cc
 */

#include "file.h"
#include <climits>
#include <stdexcept>

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
        throw std::runtime_error{"failed to get file pointer position"};
    }

    return offset;
}

void File::close() {
    if (fp_ != nullptr) {
        fclose(fp_);
        fp_ = nullptr;
    } else {
        throw std::runtime_error{"failed to close file"};
    }
}

void File::open(const std::string &path) {
    if (fp_ != nullptr) {
        throw std::runtime_error{"file pointer already in use while opening file: " + path};
    }

    const char *mode = "rb";

#ifdef _WIN32
    int rc = fopen_s(&fp_, path.c_str(), mode);
    if (rc != 0) {
        throw std::runtime_error{"failed to open file: " + path};
    }
#else
    fp_ = fopen(path.c_str(), mode);
    if (fp_ == nullptr) {
        throw std::runtime_error{"failed to open file: " + path};
    }
#endif
}

void File::seek(const int64_t offset, const int whence) {
    if (offset > LONG_MAX) {
        throw std::runtime_error{"pos out of range"};
    }

    int rc = fseek(fp_, offset, whence);
    if (rc != 0) {
        throw std::runtime_error{"failed to get file pointer position"};
    }
}

}  // namespace calq

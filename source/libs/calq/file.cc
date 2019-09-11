/**
 * @file file.cc
 */

#include "file.h"
#include <cerrno>
#include <climits>
#include <cstring>
#include "errors.h"

namespace calq {

File::File() : fp_(nullptr), size_(0) {}

File::~File() { close(); }

void File::advance(const int64_t offset) { seekFromCur(offset); }

void File::close() {
    if (fp_ != nullptr) {
        fclose(fp_);
        fp_ = nullptr;
    }
}

bool File::eof() const { return feof(fp_) != 0; }

bool File::error() const { return ferror(fp_) != 0; }

void File::open(const std::string &path, const Mode mode = Mode::READ) {
    if (fp_ != nullptr) {
        throwErrorException("Failed to open file: " + path);
    }

    const char *m = "rb";
    if (mode == Mode::WRITE) {
        m = "wb";
    }

#ifdef _WIN32
    int rc = fopen_s(&fp_, path.c_str(), m);
    if (rc != 0) {
        std::string modeStr = mode == Mode::READ ? "reading" : "writing";
        throwErrorException("Failed to open file for " + modeStr + ": " + path);
    }
#else
    fp_ = fopen(path.c_str(), m);
    if (fp_ == nullptr) {
        std::string modeStr = mode == Mode::READ ? "reading" : "writing";
        throwErrorException("Failed to open file for " + modeStr + ": " + path + " (" + strerror(errno) + ")");
    }
#endif

    // Compute file size
    seekFromEnd(0);
    size_ = tell();
    seekFromSet(0);
}

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

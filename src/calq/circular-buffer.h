/**
 * @file circular-buffer.h
 */

#ifndef CALQ_CIRCULAR_BUFFER_H_
#define CALQ_CIRCULAR_BUFFER_H_

#include <cstddef>
#include <vector>

namespace calq {

template <typename T>
class CircularBuffer {
   public:
    CircularBuffer() = delete;
    CircularBuffer(size_t size, const T& val) : pos_(0) { data_.resize(size, val); }
    CircularBuffer(const CircularBuffer&) = default;
    CircularBuffer& operator=(const CircularBuffer&) = default;
    CircularBuffer(CircularBuffer&&) = delete;
    CircularBuffer& operator=(CircularBuffer&&) = delete;
    ~CircularBuffer() = default;

    T& operator[](size_t index) { return data_[(pos_ + index) % data_.size()]; }
    const T& operator[](size_t index) const { return data_[(pos_ + index) % data_.size()]; }
    T& back() { return (*this)[pos_]; }
    const T& back() const { return (*this)[pos_]; }
    T& front() { return (*this)[(pos_ + data_.size() - 1) % data_.size()]; }
    const T& front() const { return (*this)[(pos_ + data_.size() - 1) % data_.size()]; }
    size_t size() const { return data_.size(); }
    T push(const T& val) {
        T oldVal = data_[pos_];
        data_[pos_] = val;
        pos_ = (pos_ + 1) % data_.size();
        return oldVal;
    }

   private:
    std::vector<T> data_;
    size_t pos_;
};

}  // namespace calq

#endif  // CALQ_CIRCULAR_BUFFER_H_

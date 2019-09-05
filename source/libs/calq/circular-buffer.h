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
    CircularBuffer(size_t size, const T& val) : pos_(0) { data_.resize(size, val); }

    T& operator[](size_t index) { return data_[(pos_ + index) % data_.size()]; }

    const T& operator[](size_t index) const { return data_[(pos_ + index) % data_.size()]; }

    /**
     * Leftmost value
     */
    T& back() { return (*this)[pos_]; }

    const T& back() const { return (*this)[pos_]; }

    /**
     * Rightmost value
     */
    T& front() { return (*this)[(pos_ + data_.size() - 1) % data_.size()]; }

    const T& front() const { return (*this)[(pos_ + data_.size() - 1) % data_.size()]; }

    size_t size() const { return data_.size(); }

    /**
     * Returns oldest value, deletes it and puts new value
     */
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

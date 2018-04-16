#ifndef CIRCBUFFER_H
#define CIRCBUFFER_H

#include <stddef.h>

template<typename T>
class CircularBuffer {
private:
    std::vector<T> data;
    size_t pos;
public:
    CircularBuffer(size_t size, const T& val) : pos(0) {
        data.resize(size, val);
    }

    T& operator[] (size_t index){
        return data[(pos+index)%data.size()];
    }

    const T& operator[] (size_t index) const{
        return data[(pos+index)%data.size()];
    }

    T& back (){
        return (*this)[pos];
    }

    const T& back () const{
        return (*this)[pos];
    }

    T& front (){
        return (*this)[(pos+data.size()-1)%data.size()];
    }

    const T& front () const{
        return (*this)[(pos+data.size()-1)%data.size()];
    }

    size_t size() const{
        return data.size();
    }

    T push(const T& val) {
        T oldVal = data[pos];
        data[pos] = val;
        pos = (pos+1)%data.size();
        return oldVal;
    }
};


#endif

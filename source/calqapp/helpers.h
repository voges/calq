#ifndef CALQ_HELPERS_H_
#define CALQ_HELPERS_H_

// -----------------------------------------------------------------------------

#include <functional>
#include <memory>
#include <string>
#include <utility>

// -----------------------------------------------------------------------------

namespace calqapp {

// -----------------------------------------------------------------------------

std::string currentDateAndTime();
bool fileExists(const std::string& path);
std::string fileBaseName(const std::string& path);
std::string fileNameExtension(const std::string& path);
std::string removeFileNameExtension(const std::string& path);

// ---------- make_unique (c++11) ----------------------------------------------

template<class T> struct _Unique_if
{
    typedef std::unique_ptr<T> _Single_object;
};

template<class T> struct _Unique_if<T[]>
{
    typedef std::unique_ptr<T[]> _Unknown_bound;
};

template<class T, size_t N> struct _Unique_if<T[N]>
{
    typedef void _Known_bound;
};

template<class T, class... Args>
typename _Unique_if<T>::_Single_object
make_unique(Args&& ... args){
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

template<class T>
typename _Unique_if<T>::_Unknown_bound
make_unique(size_t n){
    typedef typename std::remove_extent<T>::type U;
    return std::unique_ptr<T>(new U[n]());
}

template<class T, class... Args>
typename _Unique_if<T>::_Known_bound
make_unique(Args&& ...) = delete;

// ---------- make_unique (c++11) ----------------------------------------------

}  // namespace calq

// -----------------------------------------------------------------------------

#endif  // CALQ_HELPERS_H_

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
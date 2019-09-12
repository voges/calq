/**
 * @file exceptions.h
 */

#ifndef CALQ_EXCEPTIONS_H_
#define CALQ_EXCEPTIONS_H_

#include <exception>
#include <iostream>
#include <string>

namespace calq {

class Exception : public std::exception {
   public:
    explicit Exception(std::string msg) : msg_(std::move(msg)) {}
    Exception(const Exception& e) noexcept : msg_(e.msg_) {}
    ~Exception() noexcept override = default;
    const char* what() const noexcept override { return msg_.c_str(); }

   protected:
    std::string msg_;
};

class ErrorException : public Exception {
   public:
    explicit ErrorException(const std::string& msg) : Exception(msg) {}
};

}  // namespace calq

#endif  // CALQ_EXCEPTIONS_H_

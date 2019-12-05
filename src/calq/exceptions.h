#ifndef CALQ_EXCEPTIONS_H_
#define CALQ_EXCEPTIONS_H_

#include <exception>
#include <iostream>
#include <string>

namespace calq {

class Exception : public std::exception {
   public:
    Exception() = delete;
    explicit Exception(std::string msg) : msg_(std::move(msg)) {}
    Exception(const Exception &e) noexcept : msg_(e.msg_) {}
    Exception &operator=(const Exception &) = delete;
    Exception(Exception &&) = default;
    Exception &operator=(Exception &&) = delete;
    ~Exception() noexcept override = default;

    const char *what() const noexcept override { return msg_.c_str(); }

   protected:
    std::string msg_;
};

class ErrorException : public Exception {
   public:
    ErrorException() = delete;
    explicit ErrorException(const std::string &msg) : Exception(msg) {}
    ErrorException(const ErrorException &) = default;
    ErrorException &operator=(const ErrorException &) = delete;
    ErrorException(ErrorException &&) = default;
    ErrorException &operator=(ErrorException &&) = delete;
    ~ErrorException() override = default;
};

}  // namespace calq

#endif  // CALQ_EXCEPTIONS_H_

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
    std::string whatStr() const noexcept { return msg_; }

   protected:
    std::string msg_;
};

class ErrorException : public Exception {
   public:
    explicit ErrorException(const std::string& msg) : Exception(msg) {}
};

inline void throwErrorException(const std::string& msg) {
    std::cout.flush();
    throw calq::ErrorException(msg);
}

class ErrorExceptionReporter {
   public:
    ErrorExceptionReporter(std::string file, std::string function, const int line)
        : file_(std::move(file)), function_(std::move(function)), line_(line) {}

    void operator()(const std::string& msg) {
        std::string tmp =
            file_.substr(file_.find_last_of("/\\") + 1) + ":" + function_ + ":" + std::to_string(line_) + ": " + msg;
        // Can use the original name here, as it is still defined
        throwErrorException(tmp);
    }

   private:
    std::string file_;
    std::string function_;
    int line_;
};

}  // namespace calq

// Remove the symbol for the function, then define a new version that instead creates a stack temporary instance of
// ErrorExceptionReporter initialized with the caller.
#undef throwErrorException
#define throwErrorException calq::ErrorExceptionReporter(__FILE__, __FUNCTION__, __LINE__)

#endif  // CALQ_EXCEPTIONS_H_

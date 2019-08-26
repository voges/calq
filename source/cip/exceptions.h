#ifndef CIP_EXCEPTIONS_H_
#define CIP_EXCEPTIONS_H_

#include <exception>
#include <iostream>
#include <string>

namespace cip {

class Exception : public std::exception {
   public:
    explicit Exception(std::string msg);

    Exception(const Exception& e) noexcept;

    ~Exception() noexcept override;

    virtual std::string getMessage() const;

    const char* what() const noexcept override;

   protected:
    std::string msg_;
};

class ErrorException : public Exception {
   public:
    explicit ErrorException(const std::string& msg) : Exception(msg) {}
};

}  // namespace cip

#endif  // CIP_EXCEPTIONS_H_

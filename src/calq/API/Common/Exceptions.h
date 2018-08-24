/** @file Exceptions.h
 *  @brief This files contains the definitions of some custom exception
 *         classes.
 */

// Copyright 2015-2018 Leibniz Universitaet Hannover

#ifndef CALQ_COMMON_EXCEPTIONS_H_
#define CALQ_COMMON_EXCEPTIONS_H_

#include <exception>
#include <iostream>
#include <string>

namespace calq {

class Exception : public std::exception {
 public:
    explicit Exception(const std::string &msg);
    virtual ~Exception(void) throw();
    virtual std::string getMessage(void) const;
    virtual const char * what(void) const throw();

 protected:
    std::string msg_;
};

class ErrorException : public Exception {
 public:
    explicit ErrorException(const std::string &msg): Exception(msg) {}
};

}  // namespace calq

#endif  // CALQ_COMMON_EXCEPTIONS_H_


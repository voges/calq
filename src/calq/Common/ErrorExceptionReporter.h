/** @file ErrorExceptionReporter.h
 *  @brief This files contains the definitions of some custom exception
 *         classes.
 */

// Copyright 2015-2018 Leibniz Universitaet Hannover

#ifndef CALQ_COMMON_ERROREXCEPTIONREPORTER_H_
#define CALQ_COMMON_ERROREXCEPTIONREPORTER_H_

#include <exception>
#include <iostream>
#include <string>
#include <utility>
#include "Common/helpers.h"
#include "Common/Exceptions.h"

namespace calq {

inline void throwErrorException(const std::string &msg) {
    std::cout.flush();
    throw ErrorException(msg);
}

class ErrorExceptionReporter {
 public:
    ErrorExceptionReporter(std::string file, std::string function, const int &line) : file_(std::move(file)), function_(std::move(function)), line_(line) {}

    void operator()(const std::string &msg) {
        // std::cerr << file << ":" << function << ":" << line << ": ";
        std::string tmp = fileBaseName(file_) + ":" + function_ + ":" + std::to_string(line_) + ": " + msg;
        // Can use the original name here, as it is still defined
        throwErrorException(tmp);
    }

 private:
    std::string file_;
    std::string function_;
    int line_;
};

}  // namespace calq

// Remove the symbol for the function, then define a new version that instead
// creates a stack temporary instance of ErrorExceptionReporter initialized
// with the caller.
#undef throwErrorException
#define throwErrorException calq::ErrorExceptionReporter(__FILE__, __FUNCTION__, __LINE__)

#endif  // CALQ_COMMON_ERROREXCEPTIONREPORTER_H_


/** @file Exceptions.cc
 *  @brief This file contains the implementations of the exception classes
 *         defined in Exceptions.h.
 */

// Copyright 2015-2017 Leibniz Universitaet Hannover

#include "Exceptions.h"

namespace calq {

Exception::Exception(const std::string &msg) : msg_(msg) {
}

Exception::Exception(const Exception &e) noexcept : msg_(e.msg_) {
}

Exception::~Exception() noexcept = default;

std::string Exception::getMessage() const {
    return msg_;
}

const char* Exception::what() const noexcept {
    return msg_.c_str();
}

}  // namespace calq

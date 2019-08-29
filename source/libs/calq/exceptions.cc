/**
 * @file exceptions.cc
 */

#include "exceptions.h"

#include <utility>

namespace calq {

Exception::Exception(std::string msg) : msg_(std::move(msg)) {}

Exception::Exception(const Exception& e) noexcept : msg_(e.msg_) {}

Exception::~Exception() noexcept = default;

std::string Exception::getMessage() const { return msg_; }

const char* Exception::what() const noexcept { return msg_.c_str(); }

}  // namespace calq

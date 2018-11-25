#include "calq/exceptions.h"

namespace calq {

Exception::Exception(const std::string &msg) : msg_(msg) {}

Exception::~Exception(void) throw() {}

std::string Exception::getMessage(void) const {
    return msg_;
}

const char * Exception::what(void) const throw() {
    return msg_.c_str();
}

}  // namespace calq

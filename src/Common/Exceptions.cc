/** @file Exceptions.cc
 *  @brief This file contains the implementations of the exception classes
 *         defined in Exceptions.h.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#include "Common/Exceptions.h"

namespace calq {

Exception::Exception(const std::string &msg) : msg_(msg) {}

Exception::~Exception(void) throw () {}

std::string Exception::getMessage(void) const
{
    return msg_;
}

const char * Exception::what(void) const throw()
{
    return msg_.c_str();
}

} // namespace calq


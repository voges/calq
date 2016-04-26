#include "Exceptions.h"

Exception::Exception(std::string msg)
{
    this->msg = msg;
}

Exception::~Exception(void) throw () 
{
    // Empty
}

std::string Exception::getMessage(void) const 
{
    return msg;
}

const char * Exception::what(void) const throw() 
{
    return msg.c_str();
}


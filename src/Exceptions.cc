/** @file Exceptions.cc
 *  @brief This file contains the implementations of the exception classes
 *         defined in Exceptions.h.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: what (who)
 */

#include "Exceptions.h"

Exception::Exception(const std::string &msg)
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


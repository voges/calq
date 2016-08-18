/** @file Exceptions.cc
 *  @brief This file contains the implementations of the exception classes
 *         defined in Exceptions.h.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "Exceptions.h"

Exception::Exception(const std::string &msg): msg(msg)
{
    // empty
}

Exception::~Exception(void) throw () 
{
    // empty
}

std::string Exception::getMessage(void) const 
{
    return msg;
}

const char * Exception::what(void) const throw() 
{
    return msg.c_str();
}


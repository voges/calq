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

#include "Common/Exceptions.h"

cq::Exception::Exception(const std::string &msg): msg(msg)
{
    // empty
}

cq::Exception::~Exception(void) throw () 
{
    // empty
}

std::string cq::Exception::getMessage(void) const 
{
    return msg;
}

const char * cq::Exception::what(void) const throw() 
{
    return msg.c_str();
}


/** @file CQ.cc
 *  @brief This file contains the implementation of the CQFile class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "IO/CQ.h"
#include "Common/Exceptions.h"
#include <string.h>

cq::CQFile::CQFile(const std::string &path, const Mode &mode)
    : File(path, mode)
{
    // Check arguments
    if (path.empty() == true) {
        throwErrorException("path is empty");
    }
}

cq::CQFile::~CQFile(void)
{
    //empty 
}

size_t cq::CQFile::readHeader(void)
{
    size_t ret = 0;

    char magic[strlen(MAGIC)];
    ret += read(magic, strlen(MAGIC));
    if (strncmp(magic, MAGIC, strlen(MAGIC)) != 0) {
        throwErrorException("magic does not match");
    }

    return ret;
}

size_t cq::CQFile::writeHeader(void)
{
    size_t ret = 0;
    ret += write((char *)MAGIC, strlen(MAGIC));
    return ret;
}


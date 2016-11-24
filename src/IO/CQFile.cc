/** @file CQFile.cc
 *  @brief This file contains the implementation of the CQFile class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "IO/CQFile.h"
#include "Common/Exceptions.h"
#include <string.h>

CQFile::CQFile(const std::string &path, const CQFile::Mode &mode)
    : File(path, mode)
{
    // Check arguments
    if (path.empty() == true) {
        throwErrorException("path is empty");
    }
}

CQFile::~CQFile(void)
{
    //empty 
}

size_t CQFile::readHeader(void)
{
    size_t ret = 0;

    char magic[strlen(CQFile::MAGIC)];
    ret += read(magic, strlen(CQFile::MAGIC));
    if (strncmp(magic, CQFile::MAGIC, strlen(CQFile::MAGIC)) != 0) {
        throwErrorException("magic does not match");
    }

    return ret;
}

size_t CQFile::writeHeader(void)
{
    size_t ret = 0;
    ret += write((char *)CQFile::MAGIC, strlen(CQFile::MAGIC));
    return ret;
}


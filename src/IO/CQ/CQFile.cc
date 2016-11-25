/** @file CQFile.cc
 *  @brief This file contains the implementation of the CQFile class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "IO/CQ/CQFile.h"
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

size_t cq::CQFile::readHeader(size_t *blockSize)
{
    size_t ret = 0;

    char magic[MAGIC_LEN];
    ret += read(magic, MAGIC_LEN);
    if (strncmp(magic, MAGIC, MAGIC_LEN) != 0) {
        throwErrorException("magic does not match");
    }

    ret += readUint64((uint64_t *)blockSize);

    return ret;
}

size_t cq::CQFile::writeHeader(const size_t &blockSize)
{
    size_t ret = 0;

    ret += write((char *)MAGIC, MAGIC_LEN);
    ret += writeUint64((uint64_t)blockSize);

    return ret;
}


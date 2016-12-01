/** @file CQFile.cc
 *  @brief This file contains the implementation of the CQFile class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#include "IO/CQ/CQFile.h"

#include <math.h>
#include <string.h>

#include "Common/constants.h"
#include "Common/Exceptions.h"
#include "Common/log.h"
#include "Compressors/range/range.h"

namespace calq {

CQFile::CQFile(const std::string &path, const Mode &mode) : File(path, mode)
{
    if (path.empty() == true) {
        throwErrorException("path is empty");
    }
}

CQFile::~CQFile(void) {}

size_t CQFile::readHeader(size_t *blockSize)
{
    if (blockSize == nullptr) {
        throwErrorException("Received nullptr as argument");
    }

    size_t ret = 0;

    char magic[MAGIC_LEN];
    ret += read(magic, MAGIC_LEN);
    if (strncmp(magic, MAGIC, MAGIC_LEN) != 0) {
        throwErrorException("magic does not match");
    }

    ret += readUint64((uint64_t *)blockSize);
    CALQ_LOG("Block size: %zu", *blockSize);

    return ret;
}

size_t CQFile::writeHeader(const size_t &blockSize)
{
    if (blockSize == 0) {
        throwErrorException("blockSize must be greater than zero");
    }

    size_t ret = 0;

    ret += write((char *)MAGIC, MAGIC_LEN);
    ret += writeUint64((uint64_t)blockSize);

    return ret;
}

size_t CQFile::readBuffer(std::string *buffer)
{
    if (buffer == NULL) {
        throwErrorException("buffer is NULL");
    }
    if (buffer->empty() == false) {
        throwErrorException("buffer is not empty");
    }

    size_t ret = 0;

    uint64_t nrBlocks = 0;
    ret += readUint64(&nrBlocks);
    CALQ_LOG("Reading %zu block(s)", (size_t)nrBlocks);

    for (uint64_t i = 0; i < nrBlocks; ++i) {
        uint8_t compressed = 0;
        ret += readUint8(&compressed);
        if (compressed == 0) {
            uint32_t tmpSize = 0;
            ret += readUint32(&tmpSize);
            unsigned char *tmp = (unsigned char *)malloc(tmpSize);
            ret += read(tmp, tmpSize);
            *buffer += std::string((const char *)tmp);
            free(tmp);
            CALQ_LOG("Read uncompressed block (%u byte(s))", tmpSize);
        } else if (compressed == 1) {
            uint32_t tmpSize = 0;
            ret += readUint32(&tmpSize);
            unsigned char *tmp = (unsigned char *)malloc(tmpSize);
            ret += read(tmp, tmpSize);
            CALQ_LOG("Read compressed block (%u byte(s))", tmpSize);
            unsigned int uncompressedSize = 0;
            unsigned char *uncompressed = range_decompress_o1(tmp, &uncompressedSize);
            free(tmp);
            *buffer += std::string((const char *)uncompressed);
            free(uncompressed);
            CALQ_LOG("Uncompressed size: %u", uncompressedSize);
        } else {
            throwErrorException("Bitstream error");
        }
    }

    return ret;
}

size_t CQFile::writeBuffer(unsigned char *buffer, const size_t &bufferSize)
{
    if (buffer == NULL) {
        throwErrorException("buffer is NULL");
    }
    if (bufferSize < 1) {
        throwErrorException("bufferSize must be greater than zero");
    }

    size_t ret = 0;

    size_t nrBlocks = (size_t)ceil((double)bufferSize / (double)(1*MB));
    ret = writeUint64((uint64_t)nrBlocks);
    CALQ_LOG("Splitting buffer containing %zu byte(s) into %zu block(s)", bufferSize, nrBlocks);

    size_t encodedBytes = 0;
    while (encodedBytes < bufferSize) {
        unsigned int bytesToEncode = 0;
        if ((bufferSize - encodedBytes) > (1*MB)) {
            bytesToEncode = (1*MB);
        } else {
            bytesToEncode = bufferSize - encodedBytes;
        }

        unsigned int compressedSize = 0;
        unsigned char *compressed = range_compress_o1(buffer+encodedBytes, (unsigned int)bytesToEncode, &compressedSize);

        if (compressedSize >= bytesToEncode) {
            ret += writeUint8(0);
            ret += writeUint32(bytesToEncode);
            ret += write(buffer+encodedBytes, bytesToEncode);
        } else {
            ret += writeUint8(1);
            ret += writeUint32(compressedSize);
            ret += write(compressed, compressedSize);
        }

        encodedBytes += bytesToEncode;
        free(compressed);
    }

    return ret;
}

} // namespace calq


#include "calqapp/cq_file.h"

// -----------------------------------------------------------------------------

#include <cstring>
#include <cmath>

// -----------------------------------------------------------------------------

#include <utility>
#include <memory>
#include <vector>

// -----------------------------------------------------------------------------

#include "calqapp/range/range.h"

// -----------------------------------------------------------------------------

namespace calq {

// -----------------------------------------------------------------------------

CQFile::CQFile(const std::string& path,
               const Mode& mode
)
        : File(path, mode),
        nrReadFileFormatBytes_(0),
        nrWrittenFileFormatBytes_(0){
    if (path.empty())
    {
        throwErrorException("path is empty");
    }
}

// -----------------------------------------------------------------------------

CQFile::~CQFile() = default;

// -----------------------------------------------------------------------------

size_t CQFile::nrReadFileFormatBytes() const{
    return nrReadFileFormatBytes_;
}

// -----------------------------------------------------------------------------

size_t CQFile::nrWrittenFileFormatBytes() const{
    return nrWrittenFileFormatBytes_;
}

// -----------------------------------------------------------------------------

size_t CQFile::readHeader(size_t *blockSize){
    if (blockSize == nullptr)
    {
        throwErrorException("Received nullptr as argument");
    }

//     CALQ_LOG("Reading header");

    size_t ret = 0;

    char magic[MAGIC_LEN];
    ret += read(magic, MAGIC_LEN);
    if (strncmp(magic, MAGIC, MAGIC_LEN) != 0)
    {
        throwErrorException("magic does not match");
    }

    ret += readUint64(reinterpret_cast<uint64_t *>(blockSize));
//     CALQ_LOG("Block size: %zu", *blockSize);

    nrReadFileFormatBytes_ += ret;

    return ret;
}

// -----------------------------------------------------------------------------

size_t CQFile::readQuantizers(std::vector<std::vector<uint8_t>>* const quantizers){
    if (!quantizers->empty())
    {
        throwErrorException("quantizers is not empty");
    }

//     CALQ_LOG("Reading quantizers");

    size_t ret = 0;

    uint64_t nrQuantizers = 0;
    ret += readUint64(&nrQuantizers);

    for (uint64_t i = 0; i < nrQuantizers; ++i)
    {
        uint64_t quantizerIdx = 0;
        ret += readUint64(&quantizerIdx);

        std::vector<uint8_t> inverseLut;
        uint64_t nrInverseLutEntries = 0;
        ret += readUint64(&nrInverseLutEntries);
        for (uint64_t j = 0; j < nrInverseLutEntries; ++j)
        {
            uint8_t qualityValueIndex = 0;
            ret += readUint8(&qualityValueIndex);
            uint8_t reconstructionValue = 0;
            ret += readUint8(&reconstructionValue);
            inverseLut.push_back(reconstructionValue);
        }
        quantizers->push_back(inverseLut);
    }

    return ret;
}

// -----------------------------------------------------------------------------

size_t CQFile::readQualBlock(std::string *block){
    if (block == nullptr)
    {
        throwErrorException("block is nullptr");
    }
    else if (!block->empty())
    {
        throwErrorException("block is not empty");
    }

//     CALQ_LOG("Reading block");

    size_t ret = 0;

    uint64_t nrBlocks = 0;
    ret += readUint64(&nrBlocks);
//     CALQ_LOG("Reading %zu sub-block(s)", (size_t)nrBlocks);

    for (uint64_t i = 0; i < nrBlocks; ++i)
    {
        uint8_t compressed = 0;
        ret += readUint8(&compressed);
        if (compressed == 0)
        {
            uint32_t tmpSize = 0;
            ret += readUint32(&tmpSize);
            auto tmp = make_unique<unsigned char[]>(tmpSize);
            ret += read(tmp.get(), tmpSize);
            *block += std::string((const char *) tmp.get(), tmpSize);
//             CALQ_LOG("Read uncompressed sub-block (%u byte(s))", tmpSize);
        }
        else if (compressed == 1)
        {
            uint32_t tmpSize = 0;
            ret += readUint32(&tmpSize);
            auto tmp = make_unique<unsigned char[]>(tmpSize);
            ret += read(tmp.get(), tmpSize);
//             CALQ_LOG("Read compressed sub-block (%u byte(s))", tmpSize);
            unsigned int uncompressedSize = 0;
            unsigned char *uncompressed =
                    range_decompress_o1(tmp.get(), &uncompressedSize);
            *block +=
                    std::string((const char *) uncompressed, uncompressedSize);
            free(uncompressed);
//             CALQ_LOG("Uncompressed size was: %u", uncompressedSize);
        }
        else
        {
            throwErrorException("Bitstream error");
        }
    }

    return ret;
}

// -----------------------------------------------------------------------------

size_t CQFile::writeHeader(const size_t& blockSize){
    if (blockSize == 0)
    {
        throwErrorException("blockSize must be greater than zero");
    }

//     CALQ_LOG("Writing header");

    size_t ret = 0;

    ret += write(reinterpret_cast<const void *>(MAGIC), MAGIC_LEN);
    ret += writeUint64((uint64_t) blockSize);

    nrWrittenFileFormatBytes_ += ret;

    return ret;
}

// -----------------------------------------------------------------------------

size_t
CQFile::writeQuantizers(const std::vector<std::vector<uint8_t>>& quantizers){
    if (quantizers.empty())
    {
        throwErrorException("lut is empty");
    }

//     CALQ_LOG("Writing quantizers");

    size_t ret = 0;

    size_t nrQuantizers = quantizers.size();
    ret += writeUint64(nrQuantizers);

    for (size_t i = 0; i < quantizers.size(); ++i)
    {
        ret += writeUint64((const uint64_t&) i);
        ret += writeUint64(quantizers[i].size());
        for (size_t j = 0; j < quantizers[i].size(); ++j)
        {
            ret += writeUint8((const uint8_t&) j);
            ret += writeUint8(static_cast<const uint8_t&>(quantizers[i][j]));
        }
    }

    return ret;
}

// -----------------------------------------------------------------------------

size_t CQFile::writeQualBlock(unsigned char *block,
                              const size_t& blockSize
){
    if (block == nullptr)
    {
        throwErrorException("block is nullptr");
    }
    if (blockSize < 1)
    {
        throwErrorException("blockSize must be greater than zero");
    }

//     CALQ_LOG("Writing block");

    size_t ret = 0;

    const size_t MB = 1000000;

    auto nrBlocks = static_cast<size_t>(ceil(
            static_cast<double>(blockSize)
            / static_cast<double>(1 * MB)));
    ret = writeUint64((uint64_t) nrBlocks);
//     CALQ_LOG("Splitting block containing %zu byte(s) into %zu sub-block(s)", blockSize, nrBlocks);

    size_t encodedBytes = 0;
    while (encodedBytes < blockSize)
    {
        unsigned int bytesToEncode = 0;
        if ((blockSize - encodedBytes) > (1 * MB))
        {
            bytesToEncode = (1 * MB);
        }
        else
        {
            bytesToEncode = static_cast<unsigned int>(blockSize - encodedBytes);
        }

        unsigned int compressedSize = 0;
        unsigned char *compressed = range_compress_o1(
                block + encodedBytes,
                bytesToEncode,
                &compressedSize
        );

        if (compressedSize >= bytesToEncode)
        {
            ret += writeUint8(0);
            ret += writeUint32(bytesToEncode);
            ret += write(block + encodedBytes, bytesToEncode);
        }
        else
        {
            ret += writeUint8(1);
            ret += writeUint32(compressedSize);
            ret += write(compressed, compressedSize);
        }

        encodedBytes += bytesToEncode;
        free(compressed);
    }

    return ret;
}

// -----------------------------------------------------------------------------

}  // namespace calq

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
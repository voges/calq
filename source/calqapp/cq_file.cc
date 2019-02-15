#include "calqapp/cq_file.h"

// -----------------------------------------------------------------------------

#include <cstring>
#include <cmath>

// -----------------------------------------------------------------------------

#include <memory>

// -----------------------------------------------------------------------------

#include "calq/calq_coder.h"

// -----------------------------------------------------------------------------

#include "calqapp/range/range.h"

// -----------------------------------------------------------------------------

namespace calqapp {

// -----------------------------------------------------------------------------

CQFile::CQFile(const std::string& path,
               const Mode& mode
)
        : File(path, mode),
        nrReadFileFormatBytes_(0),
        nrWrittenFileFormatBytes_(0){
    if (path.empty()) {
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
    if (blockSize == nullptr) {
        throwErrorException("Received nullptr as argument");
    }

//     CALQ_LOG("Reading header");

    size_t ret = 0;

    char magic[MAGIC_LEN];
    ret += read(magic, MAGIC_LEN);
    if (strncmp(magic, MAGIC, MAGIC_LEN) != 0) {
        throwErrorException("magic does not match");
    }

    ret += readUint64(blockSize);
//     CALQ_LOG("Block size: %zu", *blockSize);

    nrReadFileFormatBytes_ += ret;

    return ret;
}

// -----------------------------------------------------------------------------

size_t CQFile::readQuantizers(std::vector<std::vector<uint8_t>>
                              *const quantizers
){
    if (!quantizers->empty()) {
        throwErrorException("quantizers is not empty");
    }

//     CALQ_LOG("Reading quantizers");

    size_t ret = 0;

    uint64_t nrQuantizers = 0;
    ret += readUint64(&nrQuantizers);

    for (uint64_t i = 0; i < nrQuantizers; ++i) {
        uint64_t quantizerIdx = 0;
        ret += readUint64(&quantizerIdx);

        std::vector<uint8_t> inverseLut;
        uint64_t nrInverseLutEntries = 0;
        ret += readUint64(&nrInverseLutEntries);
        for (uint64_t j = 0; j < nrInverseLutEntries; ++j) {
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

size_t CQFile::readQualBlock(std::string *block,
                             const gabac::Configuration& configuration){
    if (block == nullptr) {
        throwErrorException("block is nullptr");
    } else if (!block->empty()) {
        throwErrorException("block is not empty");
    }

//     CALQ_LOG("Reading block");
    /*
    size_t ret = 0;

    uint64_t nrBlocks = 0;
    ret += readUint64(&nrBlocks);
//     CALQ_LOG("Reading %zu sub-block(s)", (size_t)nrBlocks);

    for (uint64_t i = 0; i < nrBlocks; ++i) {
        uint8_t compressed = 0;
        ret += readUint8(&compressed);
        if (compressed == 0) {
            uint32_t tmpSize = 0;
            ret += readUint32(&tmpSize);
            auto tmp = std::unique_ptr<unsigned char[]>(
                    new unsigned char[tmpSize]
            );
            ret += read(tmp.get(), tmpSize);
            *block += std::string((const char *) tmp.get(), tmpSize);
//             CALQ_LOG("Read uncompressed sub-block (%u byte(s))", tmpSize);
        } else if (compressed == 1) {
            uint32_t tmpSize = 0;
            ret += readUint32(&tmpSize);
            auto tmp = std::unique_ptr<unsigned char[]>(
                    new unsigned char[tmpSize]
            );
            ret += read(tmp.get(), tmpSize);
//             CALQ_LOG("Read compressed sub-block (%u byte(s))", tmpSize);
            unsigned int uncompressedSize = 0;
            unsigned char *uncompressed =
                    range_decompress_o1(tmp.get(), &uncompressedSize);
            *block +=
                    std::string((const char *) uncompressed, uncompressedSize);
            free(uncompressed);
//             CALQ_LOG("Uncompressed size was: %u", uncompressedSize);
        } else {
            throwErrorException("Bitstream error");
        }
    }
    */
    size_t ret = 0;

    uint32_t tmpSize = 0;
    ret += readUint32(&tmpSize);
    auto tmp = std::unique_ptr<unsigned char[]>(new unsigned char[tmpSize]);
    ret += read(tmp.get(), tmpSize);
    std::vector<unsigned char> byteStream(tmp.get(), tmp.get() + tmpSize);

    std::vector<uint64_t> sequence;
    
    gabac::LogInfo l{&std::cout, gabac::LogInfo::LogLevel::INFO};

    gabac::decode(&byteStream, configuration, l, &sequence);

    *block = std::string(sequence.begin(), sequence.end());

    return ret;
}

// -----------------------------------------------------------------------------

size_t CQFile::writeHeader(const size_t& blockSize){
    if (blockSize == 0) {
        throwErrorException("blockSize must be greater than zero");
    }

    size_t ret = 0;

    ret += write(reinterpret_cast<const void *>(MAGIC), MAGIC_LEN);
    ret += writeUint64((uint64_t) blockSize);

    nrWrittenFileFormatBytes_ += ret;

    return ret;
}

// -----------------------------------------------------------------------------

size_t
CQFile::writeQuantizers(const std::vector<std::vector<uint8_t>>& quantizers){
    if (quantizers.empty()) {
        throwErrorException("lut is empty");
    }

//     CALQ_LOG("Writing quantizers");

    size_t ret = 0;

    size_t nrQuantizers = quantizers.size();
    ret += writeUint64(nrQuantizers);

    for (size_t i = 0; i < quantizers.size(); ++i) {
        ret += writeUint64((const uint64_t&) i);
        ret += writeUint64(quantizers[i].size());
        for (size_t j = 0; j < quantizers[i].size(); ++j) {
            ret += writeUint8((const uint8_t&) j);
            ret += writeUint8(static_cast<const uint8_t&>(quantizers[i][j]));
        }
    }

    return ret;
}

// -----------------------------------------------------------------------------

size_t CQFile::writeQualBlock(unsigned char *block,
                              const size_t& blockSize,
                              const gabac::Configuration& configuration
){
/*
if (block == nullptr) {
        throwErrorException("block is nullptr");
    }
    if (blockSize < 1) {
        throwErrorException("blockSize must be greater than zero");
    }

//     CALQ_LOG("Writing block");

    size_t ret = 0;

    const size_t MB = 1000000;

    auto nrBlocks = static_cast<size_t>(ceil(
            static_cast<double>(blockSize)
            / static_cast<double>(1 * MB)));
    ret = writeUint64((uint64_t) nrBlocks);

    size_t encodedBytes = 0;
    while (encodedBytes < blockSize) {
        unsigned int bytesToEncode = 0;
        if ((blockSize - encodedBytes) > (1 * MB)) {
            bytesToEncode = (1 * MB);
        } else {
            bytesToEncode = static_cast<unsigned int>(blockSize - encodedBytes);
        }

        unsigned int compressedSize = 0;
        unsigned char *compressed = range_compress_o1(
                block + encodedBytes,
                bytesToEncode,
                &compressedSize
        );

        if (compressedSize >= bytesToEncode) {
            ret += writeUint8(0);
            ret += writeUint32(bytesToEncode);
            ret += write(block + encodedBytes, bytesToEncode);
        } else {
            ret += writeUint8(1);
            ret += writeUint32(compressedSize);
            ret += write(compressed, compressedSize);
        }

        encodedBytes += bytesToEncode;
        free(compressed);
    }

    return ret;

    */
    if (block == nullptr) {
        throwErrorException("block is nullptr");
    }
    if (blockSize < 1) {
        throwErrorException("blockSize must be greater than zero");
    }

//     CALQ_LOG("Writing block");



    size_t ret = 0;

    unsigned int compressedSize = 0;
    std::vector<unsigned char> compressed;
    std::vector<uint64_t> toCompress;

    for(size_t i = 0; i < blockSize; i++) {
        toCompress.push_back(block[i]);
    }
        
    gabac::LogInfo l{&std::cout, gabac::LogInfo::LogLevel::INFO};

    gabac::encode(configuration, l, &toCompress, &compressed);
    compressedSize = compressed.size() * sizeof(unsigned char);
    
    ret += writeUint32(compressedSize);
    ret += write(reinterpret_cast<unsigned char*>(compressed.data()), compressedSize);

    compressed.clear();
    toCompress.clear();

    return ret;
    
}

// -----------------------------------------------------------------------------

size_t CQFile::writeBlock(const calq::EncodingOptions& opts,
                          const calq::DecodingBlock& block,
                          const calq::SideInformation& side,
                          const std::string& unmappedQualityValues_,
                          bool STREAMOUT,
                          size_t *compressedSizeMapped,
                          size_t *compressedSizeUnmapped
){
    *compressedSizeMapped = 0;
    *compressedSizeUnmapped = 0;
    // Write block parameters
    *compressedSizeMapped += this->writeUint32(side.positions[0] - 1);
    *compressedSizeMapped += this->writeUint32(
            (uint32_t) opts.qualityValueOffset
    );

    // Write inverse quantization LUTs
    *compressedSizeMapped += this->writeQuantizers(block.codeBooks);

    // Write unmapped quality values
    auto *uqv = (unsigned char *) unmappedQualityValues_.c_str();
    size_t uqvSize = unmappedQualityValues_.length();
    if (uqvSize > 0) {
        *compressedSizeUnmapped += this->writeUint8(0x01);

        if (STREAMOUT) {
            std::cerr << "unmapped qvalues:" << std::endl;

            std::cerr << unmappedQualityValues_;

            std::cerr << std::endl;
        }
        
        gabac::Configuration configuration;
        // for now simple config:
        gabac::TransformedSequenceConfiguration tsc;
        tsc.lutTransformationEnabled = 0;
        tsc.diffCodingEnabled = 0;
        tsc.binarizationId = static_cast<gabac::BinarizationId>(0);
        tsc.binarizationParameters.push_back(8);
        tsc.contextSelectionId = static_cast<gabac::ContextSelectionId>(0);
        
        configuration.transformedSequenceConfigurations.push_back(tsc);
        
        *compressedSizeUnmapped += this->writeQualBlock(uqv, uqvSize, configuration);
    } else {
        *compressedSizeUnmapped += this->writeUint8(0x00);
    }

    // Write mapped quantizer indices
    std::string mqiString;

    for (auto const& mappedQuantizerIndex : block.quantizerIndices) {
        mqiString += std::to_string(mappedQuantizerIndex);
    }

    if (STREAMOUT) {
        std::cerr << "quantizer indices:" << std::endl;

        std::cerr << mqiString;

        std::cerr << std::endl;
    }

    auto *mqi = (unsigned char *) mqiString.c_str();
    size_t mqiSize = mqiString.length();
    if (mqiSize > 0) {
        *compressedSizeMapped += this->writeUint8(0x01);

        gabac::Configuration configuration;
        // for now simple config:
        gabac::TransformedSequenceConfiguration tsc;
        tsc.lutTransformationEnabled = 0;
        tsc.diffCodingEnabled = 0;
        tsc.binarizationId = static_cast<gabac::BinarizationId>(0);
        tsc.binarizationParameters.push_back(8);
        tsc.contextSelectionId = static_cast<gabac::ContextSelectionId>(0);
        
        configuration.transformedSequenceConfigurations.push_back(tsc);

        *compressedSizeMapped += this->writeQualBlock(mqi, mqiSize, configuration);
    } else {
        *compressedSizeMapped += this->writeUint8(0x00);
    }

    // Write mapped quality value indices
    for (int i = 0; i < opts.quantizationMax - opts.quantizationMin + 1; ++i) {
        std::vector<uint8_t> mqviStream = block.stepindices[i];
        std::string mqviString;


        for (auto const& mqviInt : mqviStream) {
            mqviString += std::to_string(mqviInt);
        }
        auto *mqvi = (unsigned char *) mqviString.c_str();
        size_t mqviSize = mqviString.length();

        if (STREAMOUT) {
            std::cerr << "Step indices" << i << ":" << std::endl;


            std::cerr << mqviString;

            std::cerr << std::endl;
        }
        if (mqviSize > 0) {
            *compressedSizeMapped += this->writeUint8(0x01);

            gabac::Configuration configuration;
            // for now simple config:
            gabac::TransformedSequenceConfiguration tsc;
            tsc.lutTransformationEnabled = 0;
            tsc.diffCodingEnabled = 0;
            tsc.binarizationId = static_cast<gabac::BinarizationId>(0);
            tsc.binarizationParameters.push_back(8);
            tsc.contextSelectionId = static_cast<gabac::ContextSelectionId>(0);
            
            configuration.transformedSequenceConfigurations.push_back(tsc);

            *compressedSizeMapped += this->writeQualBlock(mqvi, mqviSize, configuration);
        } else {
            *compressedSizeMapped += this->writeUint8(0x00);
        }
    }

    return *compressedSizeUnmapped + *compressedSizeMapped;
}

// -----------------------------------------------------------------------------

size_t CQFile::readBlock(calq::DecodingBlock *out,
                         calq::SideInformation *side,
                         std::string *unmapped
){
    gabac::Configuration configuration;
    // for now simple config:
    gabac::TransformedSequenceConfiguration tsc;
    tsc.lutTransformationEnabled = 0;
    tsc.diffCodingEnabled = 0;
    tsc.binarizationId = static_cast<gabac::BinarizationId>(0);
    tsc.binarizationParameters.push_back(8);
    tsc.contextSelectionId = static_cast<gabac::ContextSelectionId>(0);
            
    configuration.transformedSequenceConfigurations.push_back(tsc);


    out->codeBooks.clear();
    out->stepindices.clear();
    out->quantizerIndices.clear();

    unmapped->clear();

    size_t ret = 0;

    std::string buffer;

    // Read block parameters
    ret += this->readUint32(&side->posOffset);
    ret += this->readUint32(reinterpret_cast<uint32_t *>(&side->qualOffset));

    // Read inverse quantization LUTs
    this->readQuantizers(&out->codeBooks);

    // Read unmapped quality values
    uint8_t uqvFlags = 0;
    ret += this->readUint8(&uqvFlags);
    if (uqvFlags & 0x01) { //NOLINT
        ret += this->readQualBlock(unmapped, configuration);
    }

    // Read mapped quantizer indices
    uint8_t mqiFlags = 0;
    ret += this->readUint8(&mqiFlags);
    if (mqiFlags & 0x1) { //NOLINT
        ret += this->readQualBlock(&buffer, configuration);
        std::copy(
                buffer.begin(),
                buffer.end(),
                std::back_inserter(out->quantizerIndices)
        );
        buffer.clear();
    }

    // Read mapped quality value indices
    for (int i = 0; i < static_cast<int>(out->codeBooks.size()); ++i) {
        out->stepindices.emplace_back();
        uint8_t mqviFlags = 0;
        ret += this->readUint8(&mqviFlags);
        if (mqviFlags & 0x1) { //NOLINT
            ret += this->readQualBlock(&buffer, configuration);
            std::copy(
                    buffer.begin(),
                    buffer.end(),
                    std::back_inserter(out->stepindices[i]));
            buffer.clear();
        }
    }

    return ret;
}

// -----------------------------------------------------------------------------

}  // namespace calqapp

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

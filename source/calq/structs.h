#ifndef CALQ_STRUCTS_H
#define CALQ_STRUCTS_H

// -----------------------------------------------------------------------------

#include <string>
#include <vector>

// -----------------------------------------------------------------------------

namespace calq {

// -----------------------------------------------------------------------------

// Side Information for encoding
struct EncodingSideInformation
{
    // Starting positions of reads in respect to genome
    std::vector<uint64_t> positions;
    std::vector<std::string> sequences; // Sequences of reads
    std::vector<std::string> cigars; // CIGARS of reads
    std::string reference; // Reference from positionStart to positionEnd
};

// -----------------------------------------------------------------------------

struct EncodingRead
{
    uint64_t posMin;
    uint64_t posMax;
    std::string qvalues;
    std::string cigar;
    std::string sequence;
    std::string reference;
};

// -----------------------------------------------------------------------------

struct EncodingBlock
{
    std::vector<std::string> qvalues;
};

// -----------------------------------------------------------------------------

struct EncodingOptions
{
    enum struct QuantizerType
    {
        NONE,
        UNIFORM,
        LLOYD_MAX
    };

    enum struct FilterType
    {
        NONE,
        GAUSS,
        RECTANGLE
    };

    enum struct Version
    {
        NONE,
        V1,
        V2
    };

    bool debug = false;
    bool squash = false;
    // EncodingOptions for only compression
    size_t filterSize = 50;
    uint8_t quantizationMin = 2;
    uint8_t quantizationMax = 8;
    uint8_t polyploidy = 2;
    uint8_t qualityValueMax = 93;
    uint8_t qualityValueMin = 0;
    uint8_t qualityValueOffset = 33;
    FilterType filterType = FilterType::GAUSS;
    QuantizerType quantizerType = QuantizerType::UNIFORM;
    Version version = Version::V1;
};

// -----------------------------------------------------------------------------

// Side Information for encoding
struct DecodingSideInformation
{
    // Starting positions of reads in respect to genome
    std::vector<uint64_t> positions;
    std::vector<std::string> cigars; // CIGARS of reads
    uint32_t posOffset;
    uint8_t qualOffset;
};

// -----------------------------------------------------------------------------

struct DecodingRead
{
    uint64_t posMin;
    std::string cigar;
};

// -----------------------------------------------------------------------------

struct DecodingBlock
{
    // Quantizer Index for each genome position in block
    std::vector<uint8_t> quantizerIndices;

    // Step index for each position in each read
    std::vector<std::vector<uint8_t>> stepindices;

    // Code books -> representative value for each step for each quantizer
    std::vector<std::vector<uint8_t>> codeBooks;
};

// -----------------------------------------------------------------------------

}

// -----------------------------------------------------------------------------

#endif //CALQ_STRUCTS_H

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

//
// Created by fabian on 1/1/19.
//

#ifndef CALQ_STRUCTS_H
#define CALQ_STRUCTS_H

#include <string>
#include <vector>

namespace calq {

// Side Information for encoding
struct EncodingSideInformation
{
    std::vector<uint64_t> positions; // Starting positions of reads in respect to genome
    std::vector<std::string> sequences; // Sequences of reads
    std::vector<std::string> cigars; // CIGARS of reads
    std::string reference; // Reference from positionStart to positionEnd
    uint64_t positionStart = 0; // Block starting position
    uint64_t positionEnd = 0; // Block ending position
};

struct EncodingRead
{
    uint64_t posMin;
    uint64_t posMax;
    std::string qvalues;
    std::string cigar;
    std::string sequence;
};

struct EncodingBlock
{
    std::vector<std::string> qvalues;
};

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

// Side Information for encoding
struct DecodingSideInformation
{
    std::vector<uint64_t> positions; // Starting positions of reads in respect to genome
    std::vector<std::string> sequences; // Sequences of reads
    std::vector<std::string> cigars; // CIGARS of reads
    std::string reference; // Reference from positionStart to positionEnd
    uint64_t positionStart = 0; // Block starting position
    uint64_t positionEnd = 0; // Block ending position
};

struct DecodingRead
{
    uint64_t posMin;
    std::string cigar;
    std::vector<uint8_t> stepindices;
};

struct DecodingBlock
{
    std::vector<uint8_t> quantizerIndices; // Quantizer Index for each genome position in block
    std::vector<std::vector<uint8_t>> stepindices; // Step index for each position in each read
    std::vector<std::vector<uint8_t>> codeBooks; // Code books -> representative value for each step for each quantizer
};

}

#endif //CALQ_STRUCTS_H

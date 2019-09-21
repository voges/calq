/**
 * @file sam-file-reader.cc
 */

#include "sam-file-reader.h"
#include <cassert>
#include <string>
#include <vector>
#include "errors.h"
#include "sam-record.h"
#include "string-helpers.h"

namespace calq {

static void parseLine(const std::string &line, std::vector<std::string> *const fields) {
    assert(fields != nullptr);

    fields->clear();
    fields->push_back("");

    int fieldCount = 0;

    for (const auto &c : line) {
        if (c == '\t' && fieldCount <= 10) {
            fieldCount++;
            fields->push_back("");
        }
        (*fields)[fieldCount] += c;
    }

    for (auto &field : *fields) {
        trim(field);
    }
}

SamFileReader::SamFileReader(const std::string &path) : header_(), ifs_() {
    ifs_.open(path, std::ifstream::in | std::ifstream::binary);
    if (!ifs_.is_open()) {
        throwErrorException("Failed to open file: " + path);
    }

    readHeader();
}

SamFileReader::~SamFileReader() { ifs_.close(); }

std::string SamFileReader::header() { return header_; }

size_t SamFileReader::readRecords(const size_t numRecords, std::list<SamRecord> *const records) {
    assert(records != nullptr);

    for (size_t i = 0; i < numRecords; i++) {
        // Read a line
        std::string line;
        getline(ifs_, line);
        if (line.empty()) {
            break;
        }

        // Parse line and construct samRecord
        std::vector<std::string> fields;
        parseLine(line, &fields);
        SamRecord samRecord(fields);
        records->push_back(samRecord);
    }

    return records->size();
}

void SamFileReader::readHeader() {
    // Set file position indicator to the beginning of the file
    auto fpos = 0;
    ifs_.seekg(std::ifstream::beg);

    while (true) {
        // Remember the current file position indicator
        fpos = ifs_.tellg();

        // Read a line
        std::string line;
        getline(ifs_, line);

        // Add the line contents to the header
        if (line[0] == '@') {
            header_ += line;
            header_ += "\n";
        } else {
            break;
        }
    }

    // Rewind file position indicator to the beginning of the the alignment section
    ifs_.seekg(fpos);
}

}  // namespace calq

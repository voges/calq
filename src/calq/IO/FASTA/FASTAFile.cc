/** @file FASTAFile.cc
 *  @brief This file contains the implementation of the FASTAFile class.
 */

// Copyright 2015-2017 Leibniz Universitaet Hannover

#include <cstring>
#include <utility>

#include "IO/FASTA/FASTAFile.h"
#include "Common/ErrorExceptionReporter.h"

namespace calq {

FASTAFile::FASTAFile(const std::string &path, const Mode &mode) : File(path, mode), line_(nullptr) {
    if (path.empty()) {
        throwErrorException("path is empty");
    }
    if (mode != FASTAFile::Mode::MODE_READ) {
        throwErrorException("Currently only MODE_READ supported");
    }

    // Usually, lines in a FASTA file should be limited to 80 chars, so 4 KB
    // should be enough
    line_ = reinterpret_cast<char*>(malloc(LINE_SIZE));
    if (line_ == nullptr) {
        throwErrorException("malloc failed");
    }

    // Parse the complete FASTA file
    parse();
}

FASTAFile::~FASTAFile() {
    free(line_);
}

void FASTAFile::parse() {
    std::string currentHeader;
    std::string currentSequence;

    while (filestream.getline(line_, LINE_SIZE)) {
        // Trim line
        size_t l = strlen(line_) - 1;
        while (l && (line_[l] == '\r' || line_[l] == '\n')) {
            line_[l--] = '\0';
        }

        if (line_[0] == '>') {
            if (!currentSequence.empty()) {
                // We have a sequence, check if we have a header
                if (currentHeader.empty()) {
                    throwErrorException("Found sequence but no header");
                }

                // We have a header, check if it is already present in our
                // references map
                if (references.find(currentHeader) != references.end()) {
                    throwErrorException("Found the same header twice");
                }

                // Everything ok, insert header-sequence pair into our
                // references map
                references.insert(std::pair<std::string, std::string>(currentHeader, currentSequence));
            }

            // Store the header and trim it: do not take the leading '>' and
            // remove everything after the first space
            currentHeader = line_ + 1;
            currentHeader = currentHeader.substr(0, currentHeader.find_first_of(' '));

            // Reset sequence
            currentSequence = "";
        } else {
            currentSequence += line_;
        }
    }

    references.insert(std::pair<std::string, std::string>(currentHeader, currentSequence));
}

}  // namespace calq


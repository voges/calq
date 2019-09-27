/**
 * @file fasta-file-reader.cc
 */

#include "fasta-file-reader.h"
#include <cassert>
#include <string>
#include "exceptions.h"

namespace calq {

FastaFileReader::FastaFileReader(const std::string &path) : ifs_() {
    ifs_.open(path, std::ifstream::in | std::ifstream::binary);
    if (!ifs_.is_open()) {
        throw ErrorException("failed to open file: " + path);
    }
}

FastaFileReader::~FastaFileReader() { ifs_.close(); }

void FastaFileReader::parse(std::vector<FastaRecord> *const fastaRecords) {
    assert(fastaRecords != nullptr);

    // Set file position indicator to the beginning of the file
    size_t fpos = ifs_.tellg();
    ifs_.seekg(std::ifstream::beg);

    std::string currentHeader;
    std::string currentSequence;

    while (true) {
        // Read a line
        std::string line;
        getline(ifs_, line);
        if (line.empty()) {
            break;
        }

        // Process line
        if (line[0] == '>') {
            // Store the previous FASTA record, if there is one
            if (!currentSequence.empty()) {
                // We have a sequence, check if we have a header
                if (currentHeader.empty()) {
                    throw ErrorException("found FASTA sequence, but no header");
                }

                FastaRecord currentFastaRecord(currentHeader, currentSequence);
                fastaRecords->push_back(currentFastaRecord);

                currentHeader = "";
                currentSequence = "";
            }

            // Store the header and trim it: remove everything after the first space
            currentHeader = line;
            currentHeader = currentHeader.substr(0, currentHeader.find_first_of(' '));
        } else {
            currentSequence += line;
        }
    }

    FastaRecord currentFastaRecord(currentHeader, currentSequence);
    fastaRecords->push_back(currentFastaRecord);

    ifs_.seekg(fpos);
}

}  // namespace calq

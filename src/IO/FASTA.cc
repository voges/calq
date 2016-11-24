/** @file FASTA.cc
 *  @brief This file contains the implementation of the FASTAFile class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "IO/FASTA.h"
#include "Common/constants.h"
#include "Common/Exceptions.h"
#include <string.h>

cq::FASTAFile::FASTAFile(const std::string &path, const FASTAFile::Mode &mode)
    : File(path, mode)
    , m_line(NULL)
{
    // Check constructor arguments
    if (path.empty() == true) {
        throwErrorException("path is empty");
    }
    if (mode != FASTAFile::MODE_READ) {
        throwErrorException("Currently only MODE_READ supported");
    }

    // Usually, lines in a FASTA file should be limited to 80 chars, so 4 KB
    // should be enough
    m_line = (char *)malloc(LINE_SIZE);

    // Parse the complete FASTA file
    parse();
}

cq::FASTAFile::~FASTAFile(void)
{
    free(m_line);
}

void cq::FASTAFile::parse(void)
{
    std::string currentHeader("");
    std::string currentSequence("");

    while (fgets(m_line, LINE_SIZE, m_fp) != NULL) {
        // Trim line
        size_t l = strlen(m_line) - 1;
        while (l && (m_line[l] == '\r' || m_line[l] == '\n')) {
            m_line[l--] = '\0';
        }

        //DEBUG("Processed %.2f%%", (double)tell()*100/(double)size());

        if (m_line[0] == '>') {
            if (currentSequence.empty() == false) {
                // We have a sequence, check if we have a header
                if (currentHeader.empty() == true) {
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
            currentHeader = m_line + 1;
            currentHeader = currentHeader.substr(0, currentHeader.find_first_of(" "));

            // Reset sequence
            currentSequence = "";
        }
        else {
            currentSequence += m_line;
        }
    }

    references.insert(std::pair<std::string, std::string>(currentHeader, currentSequence));
}


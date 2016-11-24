/** @file FASTAFile.cc
 *  @brief This file contains the implementation of the FASTAFile class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "IO/FASTAFile.h"
#include "Common/constants.h"
#include "Common/Exceptions.h"
#include <string.h>

FASTAFile::FASTAFile(const std::string &path,
                     const FASTAFile::Mode &mode)
    : File(path, mode)
    , line(NULL)
    , lineSize(sizeof(char) * (4*KB))
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
    line = (char *)malloc(lineSize);

    // Parse the complete FASTA file
    parse();
}

FASTAFile::~FASTAFile(void)
{
    free(line);
}

void FASTAFile::parse(void)
{
    std::string currentHeader("");
    std::string currentSequence("");

    while (fgets(line, lineSize, fp) != NULL) {
        // Trim line
        size_t l = strlen(line) - 1;
        while (l && (line[l] == '\r' || line[l] == '\n')) { line[l--] = '\0'; }

        //DEBUG("Processed %.2f%%", (double)tell()*100/(double)size());

        if (line[0] == '>') {
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
            currentHeader = line + 1;
            currentHeader = currentHeader.substr(0, currentHeader.find_first_of(" "));

            // Reset sequence
            currentSequence = "";
        }
        else {
            currentSequence += line;
        }
    }

    references.insert(std::pair<std::string, std::string>(currentHeader, currentSequence));
}


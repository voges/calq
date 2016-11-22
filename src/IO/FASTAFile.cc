/** @file FASTAFile.cc
 *  @brief This file contains the implementation of the FASTAFile class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#include "FASTAFile.h"
#include "Common/constants.h"
#include "Common/Exceptions.h"
#include <algorithm>
#include <stdio.h>

FASTAFile::FASTAFile(const std::string &path, const char *mode)
    : File(path, mode)
    , line(NULL)
    , lineSize(sizeof(char) * (4 * KB))
{
    // Usually, lines in a FASTA file should be limited to 80 chars, so 4 KB
    // should be enough
    line = (char *)malloc(lineSize);

    std::string currentHeader("");
    std::string currentSequence("");

    while (fgets(line, lineSize, fp) != NULL) {
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

            // Store the header and trim it: remove the leading '>',
            // newline(s), and everything after the first space
            currentHeader = line;
            currentHeader.erase(0, 1);
            currentHeader.erase(std::remove(currentHeader.begin(), currentHeader.end(), '\n'), currentHeader.end());
            currentHeader = currentHeader.substr(0, currentHeader.find_first_of("\n"));

            // Reset sequence
            currentSequence = "";
        }
        else {
            // Append sequence part and remove newline(s)
            currentSequence += line;
            currentSequence.erase(std::remove(currentSequence.begin(), currentSequence.end(), '\n'), currentSequence.end());
        }
    }

    references.insert(std::pair<std::string, std::string>(currentHeader, currentSequence));
}

FASTAFile::~FASTAFile(void)
{
    free(line);
}


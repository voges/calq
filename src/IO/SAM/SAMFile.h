/** @file SAMFile.h
 *  @brief This file contains the definition of the SAMFile class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#ifndef CALQ_IO_SAM_SAMFILE_H_
#define CALQ_IO_SAM_SAMFILE_H_

#include <string>

#include "Common/constants.h"
#include "IO/File.h"
#include "IO/SAM/SAMBlock.h"

namespace calq {

class SAMFile : public File {
public:
    explicit SAMFile(const std::string &path, const Mode &mode = MODE_READ);
    ~SAMFile(void);

    size_t nrBlocksRead(void) const;
    size_t nrMappedRecordsRead(void) const;
    size_t nrUnmappedRecordsRead(void) const;
    size_t nrRecordsRead(void) const;
    size_t readBlock(const size_t &blockSize);

    SAMBlock currentBlock;
    std::string header;

private:
    static const size_t LINE_SIZE = sizeof(char) * (1*MB);

    char *line_;
    size_t nrBlocksRead_;
    size_t nrMappedRecordsRead_;
    size_t nrUnmappedRecordsRead_;
};

} // namespace calq

#endif // CALQ_IO_SAM_SAMFILE_H_


#ifndef CALQ_SAM_FILE_H_
#define CALQ_SAM_FILE_H_

#include <chrono>
#include <string>

#include "calq/constants.h"
#include "calq/file.h"
#include "calq/sam_block.h"

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

    std::chrono::steady_clock::time_point startTime_;
};

}  // namespace calq

#endif  // CALQ_SAM_FILE_H_

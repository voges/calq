#ifndef CALQAPP_SAM_FILE_H_
#define CALQAPP_SAM_FILE_H_

// -----------------------------------------------------------------------------

#include <chrono>
#include <memory>
#include <string>

// -----------------------------------------------------------------------------

#include "calqapp/file.h"
#include "calqapp/sam_block.h"

// -----------------------------------------------------------------------------

namespace calqapp {

// -----------------------------------------------------------------------------

class SAMFile : public File
{
 public:
    explicit SAMFile(const std::string& path,
                     const Mode& mode = Mode::MODE_READ
    );
    ~SAMFile() override;

    size_t nrBlocksRead() const;
    size_t nrMappedRecordsRead() const;
    size_t nrUnmappedRecordsRead() const;
    size_t nrRecordsRead() const;
    size_t readBlock(const size_t& blockSize);

    SAMBlock currentBlock;
    std::string header;

 private:
    static const size_t LINE_SIZE = sizeof(char) * (1 * 1000000);  // 1 MB

    std::unique_ptr<char[]> line_;
    size_t nrBlocksRead_;
    size_t nrMappedRecordsRead_;
    size_t nrUnmappedRecordsRead_;

    std::chrono::steady_clock::time_point startTime_;
};

// -----------------------------------------------------------------------------

}  // namespace calqapp

// -----------------------------------------------------------------------------

#endif  // CALQAPP_SAM_FILE_H_

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

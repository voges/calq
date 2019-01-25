#ifndef CALQ_CQ_FILE_H_
#define CALQ_CQ_FILE_H_

// -----------------------------------------------------------------------------

#include <map>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------

#include "calqapp/file.h"

// -----------------------------------------------------------------------------

namespace calqapp {

// -----------------------------------------------------------------------------

class CQFile : public File
{
 public:
    CQFile(const std::string& path,
           const Mode& mode
    );
    ~CQFile() override;

    size_t nrReadFileFormatBytes() const;
    size_t nrWrittenFileFormatBytes() const;

    size_t readHeader(size_t *blockSize);
    size_t readQuantizers(std::vector<std::vector<uint8_t>> *quantizers);
    size_t readQualBlock(std::string *block);

    size_t writeHeader(const size_t& blockSize);
    size_t writeQuantizers(const std::vector<std::vector<uint8_t>>& quantizers);
    size_t writeQualBlock(unsigned char *block,
                          const size_t& blockSize
    );

 private:
    static constexpr const char *MAGIC = "CQ";
    const size_t MAGIC_LEN = 3;

    size_t nrReadFileFormatBytes_;
    size_t nrWrittenFileFormatBytes_;
};

// -----------------------------------------------------------------------------

}  // namespace calq

// -----------------------------------------------------------------------------

#endif  // CALQ_CQ_FILE_H_

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

#ifndef CALQ_CQ_FILE_H_
#define CALQ_CQ_FILE_H_

#include <map>
#include <string>

#include "calq/file.h"
#include "calq/quantizer.h"

namespace calq {

class CQFile : public File {
 public:
    CQFile(const std::string &path, const Mode &mode);
    ~CQFile(void);

    size_t nrReadFileFormatBytes(void) const;
    size_t nrWrittenFileFormatBytes(void) const;

    size_t readHeader(size_t *blockSize);
    size_t readQuantizers(std::map<int, Quantizer> *quantizers);
    size_t readQualBlock(std::string *block);

    size_t writeHeader(const size_t &blockSize);
    size_t writeQuantizers(const std::map<int, Quantizer> &quantizers);
    size_t writeQualBlock(unsigned char *block, const size_t &blockSize);

 private:
    size_t nrReadFileFormatBytes_;
    size_t nrWrittenFileFormatBytes_;
};

}  // namespace calq

#endif  // CALQ_CQ_FILE_H_

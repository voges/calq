/** @file CQFile.h
 *  @brief This file contains the definition of the CQFile class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

#ifndef CALQ_IO_CQ_CQFILE_H_
#define CALQ_IO_CQ_CQFILE_H_

#include "IO/File.h"

namespace calq {

class CQFile : public File {
public:
    CQFile(const std::string &path, const Mode &mode);
    ~CQFile(void);

    size_t nrReadHeaderBytes(void) const;
    size_t nrWrittenHeaderBytes(void) const;

    size_t readHeader(size_t *blockSize);
    size_t writeHeader(const size_t &blockSize);
    size_t readBuffer(std::string *buffer);
    size_t writeBuffer(unsigned char *buffer, const size_t &bufferSize);

private:
    static constexpr const char *MAGIC = "CQ";
    const size_t MAGIC_LEN = 3;

    size_t nrReadHeaderBytes_;
    size_t nrWrittenHeaderBytes_;
};

} // namespace calq

#endif // CALQ_IO_CQ_CQFILE_H_


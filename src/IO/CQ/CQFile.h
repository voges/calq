/** @file CQFile.h
 *  @brief This file contains the definition of the CQFile class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef CQ_CQFILE_H
#define CQ_CQFILE_H

#include "IO/File.h"

namespace cq {

class CQFile : public File {
public:
    CQFile(const std::string &path, const Mode &mode);
    ~CQFile(void);

    size_t readHeader(void);
    size_t writeHeader(void);

private:
    static constexpr const char *MAGIC = "CQ";
    const size_t MAGIC_LEN = 3;
};

}

#endif // CQ_CQFILE_H


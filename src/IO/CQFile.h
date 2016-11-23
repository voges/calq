/** @file CQFile.h
 *  @brief This file contains the definition of the CQFile class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef CQFILE_H
#define CQFILE_H

#include "IO/File.h"

class CQFile : public File {
public:
    CQFile(const std::string &path, const CQFile::Mode &mode);
    ~CQFile(void);

public:
    size_t readHeader(void);
    size_t writeHeader(void);

private:
    static constexpr const char *MAGIC = "CQ";
};

#endif // CQFILE_H


/** @file FASTAFile.h
 *  @brief This file contains the definition of the FASTAFile class.
 *  @author Jan Voges (voges)
 *  @bug No known bugs
 */

/*
 *  Changelog
 *  YYYY-MM-DD: What (who)
 */

#ifndef CQ_FASTAFILE_H
#define CQ_FASTAFILE_H

#include "Common/constants.h"
#include "IO/File.h"
#include <map>

namespace cq {

class FASTAFile : public File {
public:
    FASTAFile(const std::string &path, const Mode &mode = MODE_READ);
    ~FASTAFile(void);

public:
    std::map<std::string, std::string> references;

private:
    static const size_t LINE_SIZE = sizeof(char) * (4*KB);

private:
    void parse(void);

private:
    char *line_;
};

}

#endif // CQ_FASTAFILE_H


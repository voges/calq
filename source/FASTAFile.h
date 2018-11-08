/** @file FASTAFile.h
 *  @brief This file contains the definition of the FASTAFile class.
 */

// Copyright 2015-2017 Leibniz Universitaet Hannover

#ifndef CALQ_IO_FASTA_FASTAFILE_H_
#define CALQ_IO_FASTA_FASTAFILE_H_

#include <map>
#include <string>
#include <memory>

#include "constants.h"
#include "File.h"

namespace calq {

class FASTAFile : public File {
 public:
    explicit FASTAFile(const std::string &path, const Mode &mode = Mode::MODE_READ);
    ~FASTAFile() override;

    std::map<std::string, std::string> references;

 private:
    static const size_t LINE_SIZE = sizeof(char) * (4*KB);

    void parse();

    std::unique_ptr<char[]> line_;
};

}  // namespace calq

#endif  // CALQ_IO_FASTA_FASTAFILE_H_

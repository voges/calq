#ifndef CALQ_FASTA_FILE_H_
#define CALQ_FASTA_FILE_H_

#include <map>
#include <string>
#include <memory>

#include "calq/constants.h"
#include "calq/file.h"

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

#endif  // CALQ_FASTA_FILE_H_

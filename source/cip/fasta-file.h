#ifndef CALQAPP_FASTA_FILE_H_
#define CALQAPP_FASTA_FILE_H_

// -----------------------------------------------------------------------------

#include <map>
#include <memory>
#include <string>

// -----------------------------------------------------------------------------

#include "file.h"

// -----------------------------------------------------------------------------

namespace cip {

// -----------------------------------------------------------------------------

class FASTAFile : public File {
   public:
    explicit FASTAFile(const std::string& path,
                       const Mode& mode = Mode::MODE_READ);
    ~FASTAFile() override;

    std::map<std::string, std::string> references;

    std::string getReferencesInRange(const std::string& header,
                                     const size_t& start, const size_t& end);

   private:
    static const size_t LINE_SIZE = sizeof(char) * (4 * 1000);  // 4 KB

    void parse();

    std::unique_ptr<char[]> line_;
};

// -----------------------------------------------------------------------------

}  // namespace cip

// -----------------------------------------------------------------------------

#endif  // CALQAPP_FASTA_FILE_H_

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

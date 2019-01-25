#include <cstdint>

// -----------------------------------------------------------------------------

#include <memory>
#include <vector>
#include <string>

// -----------------------------------------------------------------------------

#include "calq/calq_coder.h"

// -----------------------------------------------------------------------------

namespace calqapp {

// -----------------------------------------------------------------------------

class SAMFile;

class FASTAFile;

// -----------------------------------------------------------------------------

struct UnmappedInformation
{
    std::vector<bool> mappedFlags;
    std::vector<std::string> unmappedQualityScores;
};

// -----------------------------------------------------------------------------

class SAMFileHandler
{
 public:

    SAMFileHandler(const std::string& inputFileName,
                   const std::string& referenceFileName
    );
    ~SAMFileHandler();

    size_t readBlock(const size_t& blocksize);
    void getMappedBlock(calq::EncodingBlock *var);
    void getUnmappedBlock(UnmappedInformation *var);
    void getSideInformation(calq::SideInformation *var);
    size_t nrBlocksRead() const;
    size_t nrMappedRecordsRead() const;
    size_t nrUnmappedRecordsRead() const;
    size_t nrRecordsRead() const;

 private:
    std::unique_ptr<SAMFile> samFile_;
    std::unique_ptr<FASTAFile> fastaFile;

    calq::SideInformation side;
    calq::EncodingBlock encBlock;
    UnmappedInformation unmapped;
    size_t refStart;
    size_t refEnd;
    std::string rname;
};

// -----------------------------------------------------------------------------

} // namespace calq

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
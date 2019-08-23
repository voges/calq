#ifndef CALQAPP_SAMFILEHANDLER_H_
#define CALQAPP_SAMFILEHANDLER_H_

// -----------------------------------------------------------------------------

#include <cstdint>

// -----------------------------------------------------------------------------

#include <memory>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------

#include "calq/calq-codec.h"

// -----------------------------------------------------------------------------

namespace cip {

// -----------------------------------------------------------------------------

class SAMFile;
class FASTAFile;

// -----------------------------------------------------------------------------

struct UnmappedInformation {
    std::vector<bool> mappedFlags;
    std::vector<std::string> unmappedQualityScores;
};

// -----------------------------------------------------------------------------

class SAMFileHandler {
   public:
    SAMFileHandler(const std::string& inputFileName,
                   const std::string& referenceFileName);
    ~SAMFileHandler();

    size_t readBlock(const size_t& blocksize);
    void getMappedBlock(calq::EncodingBlock* var);
    void getUnmappedBlock(UnmappedInformation* var);
    void getSideInformation(calq::SideInformation* var);
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

}  // namespace cip

// -----------------------------------------------------------------------------

#endif  // CALQAPP_SAMFILEHANDLER_H_

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

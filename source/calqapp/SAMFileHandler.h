#include <cstdint>

// -----------------------------------------------------------------------------

#include <memory>
#include <vector>
#include <string>

// -----------------------------------------------------------------------------

namespace calq {

// -----------------------------------------------------------------------------

class SAMFile;

// -----------------------------------------------------------------------------

class SAMFileHandler
{
 public:

    SAMFileHandler(const std::string& inputFileName);
    ~SAMFileHandler();

    size_t readBlock(const size_t& blocksize);
    void getPositions(std::vector<uint64_t> *var);
    void getMappedFlags(std::vector<bool> *var);
    void getSequences(std::vector<std::string> *var);
    void getCigars(std::vector<std::string> *var);
    void getMappedQualityScores(std::vector<std::string> *var);
    void getUnmappedQualityScores(std::vector<std::string> *var);
    size_t getRefStart();
    size_t getRefEnd();
    void getRname(std::string *var);
    // paramStruct& getOtherParams();


 private:
    std::unique_ptr<SAMFile> samFile_;

    std::vector<uint64_t> positions;
    std::vector<std::string> sequences;
    std::vector<std::string> cigars;
    std::vector<bool> mappedFlags;
    std::vector<std::string> mappedQualityScores;
    std::vector<std::string> unmappedQualityScores;
    size_t refStart;
    size_t refEnd;
    std::string rname;
    //paramStruct& otherParams;
};

// -----------------------------------------------------------------------------

} // namespace calq

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
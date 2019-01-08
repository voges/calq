/** @file samFileReader.h
 *  @brief This file constains the definition of the SAMFileHandler class.
 */

#include <memory>
#include <vector>
#include <stdint.h>
#include <string>

namespace calq {

class SAMFile;

class SAMFileHandler {
  public:
  	
  	SAMFileHandler(const std::string inputFileName);
  	~SAMFileHandler();

    size_t readBlock(const size_t &blocksize);
  	std::vector<uint64_t> getPositions();
  	std::vector<std::string> getSequences();
  	std::vector<std::string> getCigars();
  	std::vector<std::string> getMappedQualityScores();
    std::string getUnmappedQualityScores();

  	// paramStruct& getOtherParams();

  
  private:
  	std::unique_ptr<SAMFile> samFile_;

  	std::vector<uint64_t> positions;
  	std::vector<std::string> sequences;
  	std::vector<std::string> cigars;
  	std::vector<std::string> mappedQualityScores;
    std::string unmappedQualityScores;
  	//paramStruct& otherParams;
};

} // namespace calq
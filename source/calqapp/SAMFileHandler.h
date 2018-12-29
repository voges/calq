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
  	std::vector<std::string> getQualityScores();
  	// paramStruct& getOtherParams();
  	std::vector<uint8_t> getQuantizerIndices();
  	std::vector<std::vector<uint8_t>> getStepIndices();
  	std::vector<std::vector<uint8_t>> getCodeBooks();
  
  private:
  	std::unique_ptr<SAMFile> samFile_;

  	std::vector<uint64_t> positions;
  	std::vector<std::string> sequences;
  	std::vector<std::string> cigars;
  	std::vector<std::string> qualityScores;
  	//paramStruct& otherParams;
};

} // namespace calq
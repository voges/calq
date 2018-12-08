/** @file samFileReader.h
 *  @brief This file constains the definition of the SAMFileHandler class.
 */

#include <memory>
#include <vector>
#include <stdint.h>
#include <string>

class SamFileHandler {
  public:
  	
  	SamFileHandler(const std::string inputFileName);
  	~SamFileHandler();

  	std::vector<uint64_t> getPositions();
  	std::vector<std::string> getSequences();
  	std::vector<uint64_t> getCigars();
  	std::vector<uint64_t> getQualityScores();
  	// paramStruct& getOtherParams();
  	std::vector<uint8_t> getQuantizerIndices();
  	std::vector<std::vector<uint8_t>> getStepIndices();
  	std::vector<std::vector<uint8_t>> getCodeBooks();
  
  private:
  	// std::unique_ptr<SAMFile> samFile_;

  	std::vector<uint64_t> positions;
  	std::vector<std::string> sequences;
  	std::vector<uint64_t> cigars;
  	std::vector<uint64_t> qualityScores;
  	//paramStruct& otherParams;
  	std::vector<uint8_t> quantizerIndices;
  	std::vector<std::vector<uint8_t>> stepIndices;
  	std::vector<std::vector<uint8_t>> codeBooks;
};

/** @file SamFileHandler.cc
 *  @brief This file contains the implementation of the SAMFileHandler class->
 */

#include "calqapp/SAMFileHandler.h"
#include "calq/helpers.h"
#include "calq/sam_file.h"

namespace calq {

SAMFileHandler::SAMFileHandler(const std::string inputFileName) :  samFile_(nullptr),
																   positions(),
																   sequences(),
																   cigars(),
																   mappedQualityScores(),
																   unmappedQualityScores()
																// otherParams(),
{
	samFile_ = calq::make_unique<SAMFile>(inputFileName);
}

SAMFileHandler::~SAMFileHandler() = default;

size_t SAMFileHandler::readBlock(const size_t &blockSize) {
	size_t returnCode = samFile_->readBlock(blockSize);


	for (auto const &samRecord : samFile_->currentBlock.records) {
		// TO-DO:	differentiate between mapped and unmapped
		//			generate one qualityScore stream for unmapped/mapped
		// 			mapped -> directly to calq lib
		//			unmapped -> to gabac (?) - possibly quantized
		if (samRecord.isMapped() == true) {
			positions.push_back(samRecord.pos);
			sequences.push_back(samRecord.seq);
			cigars.push_back(samRecord.cigar);
			mappedQualityScores.push_back(samRecord.qual);
		} else {
			unmappedQualityScores.push_back(samRecord.qual);
		}
	}
	return returnCode;
}

std::vector<uint64_t> SAMFileHandler::getPositions() {
	return this->positions;
}

std::vector<std::string> SAMFileHandler::getSequences() {
	return this->sequences;
}

std::vector<std::string> SAMFileHandler::getCigars() {
	return this->cigars;
}

std::vector<std::string> SAMFileHandler::getMappedQualityScores() {
	return this->mappedQualityScores;
}

std::vector<std::string> SAMFileHandler::getUnmappedQualityScores() {
	return this->unmappedQualityScores;
}
/*
paramStruct& SamFileHandler::getOtherParams() {
	return this->otherParams;
}
*/

} // namespace calq

/** @file SamFileHandler.cc
 *  @brief This file contains the implementation of the SAMFileHandler class->
 */

#include "calqapp/SAMFileHandler.h"
#include "calq/helpers.h"
#include "calq/sam_file.h"

namespace calq {

SAMFileHandler::SAMFileHandler(const std::string inputFileName) :  samFile_(nullptr)
																// positions(),
																// sequences(),
																// cigars(),
																// qualityScores(),
																// otherParams(),
																// quantizerIndices(),
																// stepIndices(),
																// codeBooks() 
{
	samFile_ = calq::make_unique<SAMFile>(inputFileName);

}

SAMFileHandler::~SAMFileHandler() = default;

size_t SAMFileHandler::readBlock(const size_t &blockSize) {
	size_t returnCode = samFile_->readBlock(blockSize);
	for (auto const &samRecord : samFile_->currentBlock.records) {
		positions.push_back(samRecord.pos);
		sequences.push_back(samRecord.seq);
		cigars.push_back(samRecord.cigar);
		qualityScores.push_back(samRecord.qual);
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

std::vector<std::string> SAMFileHandler::getQualityScores() {
	return this->qualityScores;
}
/*
paramStruct& SamFileHandler::getOtherParams() {
	return this->otherParams;
}
*/

} // namespace calq

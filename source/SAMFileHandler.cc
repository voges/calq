/** @file SamFileHandler.cc
 *  @brief This file contains the implementation of the SAMFileHandler class->
 */

#include "SAMFileHandler.h"

SamFileHandler::SamFileHandler(const std::string inputFileName) //:  samFile_(nullptr)
																// positions(),
																// sequences(),
																// cigars(),
																// qualityScores(),
																// otherParams(),
																// quantizerIndices(),
																// stepIndices(),
																// codeBooks()
																 {}

std::vector<uint64_t> SamFileHandler::getPositions() {
	return this->positions;
}

std::vector<std::string> SamFileHandler::getSequences() {
	return this->sequences;
}

std::vector<uint64_t> SamFileHandler::getCigars() {
	return this->cigars;
}

std::vector<uint64_t> SamFileHandler::getQualityScores() {
	return this->qualityScores;
}
/*
paramStruct& SamFileHandler::getOtherParams() {
	return this->otherParams;
}
*/
std::vector<uint8_t> SamFileHandler::getQuantizerIndices() {
	return this->quantizerIndices;
}

std::vector<std::vector<uint8_t>> SamFileHandler::getStepIndices() {
	return this->stepIndices;
}

std::vector<std::vector<uint8_t>> SamFileHandler::getCodeBooks() {
	return this->codeBooks;
}

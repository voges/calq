/** @file SamFileHandler.cc
 *  @brief This file contains the implementation of the SAMFileHandler class->
 */

#include "calqapp/SAMFileHandler.h"
#include "calq/helpers.h"
#include "calq/sam_file.h"
#include "calq/calq_encoder.h"

namespace calq {

SAMFileHandler::SAMFileHandler(const std::string inputFileName) :  samFile_(nullptr),
																   positions(),
																   sequences(),
																   cigars(),
																   mappedQualityScores(),
																   unmappedQualityScores(),
																   refStart(),
																   refEnd(),
																   rname()
																// otherParams(),
{
	samFile_ = calq::make_unique<SAMFile>(inputFileName);
}

SAMFileHandler::~SAMFileHandler() = default;

uint32_t computeRefLength(const std::string& cigar){
    // Compute 0-based first position and 0-based last position this record
    // is mapped to on the reference used for alignment
    uint32_t posMax = 0;

    size_t cigarIdx = 0;
    size_t cigarLen = cigar.length();
    uint32_t opLen = 0;  // length of current CIGAR operation

    for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++)
    {
        if (isdigit(cigar[cigarIdx]))
        {
            opLen = opLen * 10 + (uint32_t) cigar[cigarIdx] - (uint32_t) '0';
            continue;
        }
        switch (cigar[cigarIdx])
        {
            case 'M':
            case '=':
            case 'X':
                posMax += opLen;
                break;
            case 'I':
            case 'S':
                break;
            case 'D':
            case 'N':
                posMax += opLen;
                break;
            case 'H':
            case 'P':
                break;  // these have been clipped
            default:
                throwErrorException("Bad CIGAR string");
        }
        opLen = 0;
    }
    return posMax;
}

size_t SAMFileHandler::readBlock(const size_t &blockSize) {
	size_t returnCode = samFile_->readBlock(blockSize);

	refStart = samFile_->currentBlock.records[0].pos;
	refEnd = 0;
	rname = samFile_->currentBlock.records[0].rname;
	for (auto const &samRecord : samFile_->currentBlock.records) {
		// TO-DO:	differentiate between mapped and unmapped
		//			generate one qualityScore stream for unmapped/mapped
		// 			mapped -> directly to calq lib
		//			unmapped -> to gabac (?) - possibly quantized
		if (samRecord.isMapped()) {
			positions.push_back(samRecord.pos);
			sequences.push_back(samRecord.seq);
			cigars.push_back(samRecord.cigar);
			mappedQualityScores.push_back(samRecord.qual);
			if ((samRecord.pos + computeRefLength(samRecord.cigar)) > refEnd) {
				refEnd = samRecord.pos + computeRefLength(samRecord.cigar);
			}
		} else {
			unmappedQualityScores+=samRecord.qual;
		}
	}
	return returnCode;
}

void SAMFileHandler::getPositions(std::vector<uint64_t>* var) {
	var->clear();
	var->swap(this->positions);
}

void SAMFileHandler::getSequences(std::vector<std::string>* var) {
	var->clear();
	var->swap(this->sequences);
}

void SAMFileHandler::getCigars(std::vector<std::string>* var) {
	var->clear();
	var->swap(this->cigars);
}

void SAMFileHandler::getMappedQualityScores(std::vector<std::string>* var) {
	var->clear();
	var->swap(this->mappedQualityScores);
}

void SAMFileHandler::getUnmappedQualityScores(std::string* var) {
	var->clear();
	var->swap(this->unmappedQualityScores);
}

size_t SAMFileHandler::getRefStart() {
	return this->refStart;
}

size_t SAMFileHandler::getRefEnd() {
	return this->refEnd;
}

void SAMFileHandler::getRname(std::string* var) {
	var->clear();
	var->swap(this->rname);
}
/*
paramStruct& SamFileHandler::getOtherParams() {
	return this->otherParams;
}
*/

} // namespace calq

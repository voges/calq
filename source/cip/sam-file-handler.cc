#include "sam-file-handler.h"
#include <algorithm>

// -----------------------------------------------------------------------------

#include "fasta-file.h"
#include "sam-file.h"

// -----------------------------------------------------------------------------

namespace cip {

// -----------------------------------------------------------------------------

SAMFileHandler::SAMFileHandler(const std::string& inputFileName,
                               const std::string& referenceFileName)
    : samFile_(nullptr),
      fastaFile(nullptr),
      side(),
      encBlock(),
      unmapped(),
      refStart(),
      refEnd(),
      rname() {
    samFile_ = std::unique_ptr<SAMFile>(new SAMFile(inputFileName));
    if (!referenceFileName.empty()) {
        fastaFile =
            std::unique_ptr<FASTAFile>(new FASTAFile(referenceFileName));
    }
}

// -----------------------------------------------------------------------------

SAMFileHandler::~SAMFileHandler() = default;

// -----------------------------------------------------------------------------

uint32_t computeRefLength(const std::string& cigar) {
    // Compute 0-based first position and 0-based last position this record
    // is mapped to on the reference used for alignment
    uint32_t posMax = 0;

    size_t cigarIdx = 0;
    size_t cigarLen = cigar.length();
    uint32_t opLen = 0;  // length of current CIGAR operation

    for (cigarIdx = 0; cigarIdx < cigarLen; cigarIdx++) {
        if (isdigit(cigar[cigarIdx])) {
            opLen = opLen * 10 + (uint32_t)cigar[cigarIdx] - (uint32_t)'0';
            continue;
        }
        switch (cigar[cigarIdx]) {
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

// -----------------------------------------------------------------------------

size_t SAMFileHandler::readBlock(const size_t& blockSize) {
    this->unmapped.mappedFlags.clear();
    this->unmapped.unmappedQualityScores.clear();
    this->encBlock.qvalues.clear();
    side.positions.clear();
    side.cigars.clear();
    side.sequences.clear();

    size_t returnCode = samFile_->readBlock(blockSize);
    refEnd = 0;
    rname = samFile_->currentBlock.records[0].rname;
    for (auto const& samRecord : samFile_->currentBlock.records) {
        unmapped.mappedFlags.push_back(samRecord.isMapped());
        if (samRecord.isMapped()) {
            if (side.positions.empty()) {
                refStart = samRecord.pos;
            }
            side.positions.push_back(samRecord.pos);
            side.sequences.push_back(samRecord.seq);
            side.cigars.push_back(samRecord.cigar);
            encBlock.qvalues.push_back(samRecord.qual);
            auto endTmp = samRecord.pos + computeRefLength(samRecord.cigar);
            if (endTmp > refEnd) {
                refEnd = endTmp;
            }
        } else {
            unmapped.unmappedQualityScores.push_back(samRecord.qual);
        }
    }

    if (!fastaFile) {
        return returnCode;
    }

    side.reference = fastaFile->getReferencesInRange(rname, refStart, refEnd);
    std::transform(side.reference.begin(), side.reference.end(),
                   side.reference.begin(), ::toupper);

    if (side.positions.empty()) {
        side.posOffset = 0;
    } else {
        side.posOffset = side.positions[0];
    }

    return returnCode;
}

// -----------------------------------------------------------------------------

void SAMFileHandler::getMappedBlock(calq::EncodingBlock* var) {
    var->qvalues.swap(encBlock.qvalues);
}

// -----------------------------------------------------------------------------

void SAMFileHandler::getUnmappedBlock(UnmappedInformation* var) {
    var->mappedFlags.swap(unmapped.mappedFlags);
    var->unmappedQualityScores.swap(unmapped.unmappedQualityScores);
}

// -----------------------------------------------------------------------------

void SAMFileHandler::getSideInformation(calq::SideInformation* var) {
    side.cigars.swap(var->cigars);
    side.sequences.swap(var->sequences);
    side.positions.swap(var->positions);
    side.reference.swap(var->reference);
    var->posOffset = side.posOffset;
    var->qualOffset = side.qualOffset;
}

// -----------------------------------------------------------------------------

size_t SAMFileHandler::nrBlocksRead() const {
    return this->samFile_->nrBlocksRead();
}

// -----------------------------------------------------------------------------

size_t SAMFileHandler::nrMappedRecordsRead() const {
    return this->samFile_->nrMappedRecordsRead();
}

// -----------------------------------------------------------------------------

size_t SAMFileHandler::nrUnmappedRecordsRead() const {
    return this->samFile_->nrUnmappedRecordsRead();
}

// -----------------------------------------------------------------------------

size_t SAMFileHandler::nrRecordsRead() const {
    return this->samFile_->nrRecordsRead();
}

// -----------------------------------------------------------------------------

}  // namespace cip

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
#include "sam-file.h"
#include <cstring>

namespace cip {

static void parseLine(char *fields[SamRecord::NUM_FIELDS], char *line) {
    char *c = line;
    char *pc = c;
    int f = 0;

    while (true) {
        if (*c == '\t' || *c == '\0') {
            *c = '\0';
            if (f == 0) {
                fields[f] = pc;
            }
            if (f == 1) {
                fields[f] = pc;
            }
            if (f == 2) {
                fields[f] = pc;
            }
            if (f == 3) {
                fields[f] = pc;
            }
            if (f == 4) {
                fields[f] = pc;
            }
            if (f == 5) {
                fields[f] = pc;
            }
            if (f == 6) {
                fields[f] = pc;
            }
            if (f == 7) {
                fields[f] = pc;
            }
            if (f == 8) {
                fields[f] = pc;
            }
            if (f == 9) {
                fields[f] = pc;
            }
            if (f == 10) {
                fields[f] = pc;
            }
            if (f == 11) {
                fields[f] = pc;
            }
            f++;
            if (f == 12) {
                break;
            }
            pc = c + 1;
        }
        c++;
    }

    // if (f == 11) { fields[f] = pc; }
}

SAMFile::SAMFile(const std::string &path, const Mode &mode)
    : File(path, mode),
      currentBlock(),
      header(""),
      line_(nullptr),
      nrBlocksRead_(0),
      nrMappedRecordsRead_(0),
      nrUnmappedRecordsRead_(0),
      startTime_(std::chrono::steady_clock::now()) {
    if (path.empty()) {
        throwErrorException("path is empty");
    }
    if (mode != Mode::MODE_READ) {
        throwErrorException("Currently only MODE_READ supported");
    }

    try {
        // 1 million chars should be enough
        line_ = std::unique_ptr<char[]>(new char[LINE_SIZE]);
    } catch (std::exception &e) {
        throwErrorException(std::string("New failed: ") + e.what());
    }

    // Read SAM header
    size_t fpos = 0;
    for (;;) {
        fpos = tell();
        if (readLine(line_.get(), LINE_SIZE)) {
            // Trim line
            size_t l = strlen(line_.get()) - 1;
            while (l && (line_[l] == '\r' || line_[l] == '\n')) {
                line_[l--] = '\0';
            }

            if (line_[0] == '@') {
                header += line_.get();
                header += "\n";
            } else {
                break;
            }
        } else {
            throwErrorException("Could not read SAM header");
        }
    }
    seek(fpos);  // rewind to the begin of the alignment section
    if (header.empty()) {
        std::cout << "CALQ: Warning: No SAM header found" << std::endl;
    }
}

SAMFile::~SAMFile() = default;

size_t SAMFile::nrBlocksRead() const { return nrBlocksRead_; }

size_t SAMFile::nrMappedRecordsRead() const { return nrMappedRecordsRead_; }

size_t SAMFile::nrUnmappedRecordsRead() const { return nrUnmappedRecordsRead_; }

size_t SAMFile::nrRecordsRead() const { return (nrMappedRecordsRead_ + nrUnmappedRecordsRead_); }

size_t SAMFile::readBlock(const size_t &blockSize) {
    if (blockSize < 1) {
        throwErrorException("blockSize must be greater than zero");
    }

    currentBlock.reset();

    std::string rnamePrev;
    uint32_t posPrev = 0;

    for (size_t i = 0; i < blockSize; i++) {
        size_t fpos = tell();
        if (readLine(line_.get(), LINE_SIZE)) {
            // Trim line
            size_t l = strlen(line_.get()) - 1;
            while (l && (line_[l] == '\r' || line_[l] == '\n')) {
                line_[l--] = '\0';
            }

            // Parse line and construct samRecord
            char *fields[SamRecord::NUM_FIELDS];
            parseLine(fields, line_.get());
            SamRecord samRecord(fields);

            if (samRecord.isMapped()) {
                if (rnamePrev.empty()) {
                    // This is the first mapped record in this block; just store
                    // its RNAME and POS and add it to the current block
                    rnamePrev = samRecord.rname;
                    posPrev = samRecord.pos;
                    currentBlock.records.push_back(samRecord);
                    currentBlock.nrMappedRecords_++;
                } else {
                    // We already have a mapped record in this block
                    if (rnamePrev == samRecord.rname) {
                        // RNAME didn't change, check POS
                        if (samRecord.pos >= posPrev) {
                            // Everything fits, just update posPrev and push
                            // the samRecord to the current block
                            posPrev = samRecord.pos;
                            currentBlock.records.push_back(samRecord);
                            currentBlock.nrMappedRecords_++;
                        } else {
                            throwErrorException("SAM file is not sorted");
                        }
                    } else {
                        // RNAME changed, seek back and break
                        seek(fpos);
                        std::cout << "CALQ: Warning: RNAME changed - read only " << currentBlock.nrRecords() << "record(s) (" << blockSize << " record(s) requested)" << std::endl;
                        break;
                    }
                }
            } else {
                currentBlock.records.push_back(samRecord);
                currentBlock.nrUnmappedRecords_++;
            }
        } else {
            std::cout << "CALQ: Warning: Truncated block - read only " << currentBlock.nrRecords() << " record(s) (" << blockSize << " record(s) requested) - reached EOF" << std::endl;
            break;
        }
    }

    if (currentBlock.nrRecords() > 0) {
        nrBlocksRead_++;
        nrMappedRecordsRead_ += currentBlock.nrMappedRecords();
        nrUnmappedRecordsRead_ += currentBlock.nrUnmappedRecords();
    }

    auto elapsedTime = std::chrono::steady_clock::now() - startTime_;
    auto elapsedTimeS = std::chrono::duration_cast<std::chrono::seconds>(elapsedTime).count();
    double elapsedTimeM = static_cast<double>(elapsedTimeS) / 60.0;
    double processedPercentage = (static_cast<double>(tell()) / static_cast<double>(size())) * 100.0;
    auto remainingPercentage = 100 - processedPercentage;
    std::cout << "Processed " << processedPercentage << "% (elapsed: " << elapsedTimeM << " m), remaining: " << remainingPercentage << "% (~" << elapsedTimeM * (remainingPercentage / processedPercentage) << "m)" << std::endl;

    return currentBlock.nrRecords();
}

}  // namespace cip

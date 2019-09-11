/**
 * @file sam-file-reader.h
 */

#ifndef CALQ_SAM_FILE_READER_H_
#define CALQ_SAM_FILE_READER_H_

#include <list>
#include <string>
#include "file-line-reader.h"
#include "sam-record.h"

namespace calq {

class SamFileReader : public FileLineReader {
   public:
    explicit SamFileReader(const std::string &path);
    std::string header();
    size_t readRecords(size_t numRecords, std::list<SamRecord> *records);

   private:
    void readHeader();
    std::string header_;
};

}  // namespace calq

#endif  // CALQ_SAM_FILE_READER_H_

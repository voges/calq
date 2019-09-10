/**
 * @file sam-file-reader.h
 */

#ifndef CALQ_SAM_FILE_READER_H_
#define CALQ_SAM_FILE_READER_H_

#include <list>
#include <string>

#include "file-reader.h"
#include "sam-record.h"

namespace calq {

class SamFileReader : public FileReader {
   public:
    explicit SamFileReader(const std::string &path);

    ~SamFileReader() override;

    size_t readRecords(size_t numRecords, std::list<SamRecord> *records);

   public:
    std::string header;

   private:
    void readHeader();
};

}  // namespace calq

#endif  // CALQ_SAM_FILE_READER_H_

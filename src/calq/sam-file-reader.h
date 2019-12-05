#ifndef CALQ_SAM_FILE_READER_H_
#define CALQ_SAM_FILE_READER_H_

#include <fstream>
#include <list>
#include <string>
#include "sam-record.h"

namespace calq {

class SamFileReader {
   public:
    SamFileReader() = delete;
    explicit SamFileReader(const std::string &path);
    SamFileReader(const SamFileReader &) = delete;
    SamFileReader &operator=(const SamFileReader &) = delete;
    SamFileReader(SamFileReader &&) = delete;
    SamFileReader &operator=(SamFileReader &&) = delete;
    ~SamFileReader();

    std::string header();
    size_t readRecords(size_t numRecords, std::list<SamRecord> *records);

   private:
    void readHeader();

    std::string header_;
    std::ifstream ifs_;
};

}  // namespace calq

#endif  // CALQ_SAM_FILE_READER_H_

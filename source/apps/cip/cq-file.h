#ifndef CALQAPP_CQ_FILE_H_
#define CALQAPP_CQ_FILE_H_

#include <map>
#include <string>
#include <vector>
#include "file.h"
//#include "gabac/configuration.h"
//#include "gabac/gabac.h"

namespace calq {
struct EncodingOptions;
struct DecodingBlock;
struct SideInformation;
}  // namespace calq

namespace cip {

class CQFile : public File {
   public:
    CQFile(const std::string& path, const Mode& mode);
    ~CQFile() override;

//    size_t readBlock(calq::DecodingBlock* out, calq::SideInformation* side, std::string* unmapped,
//                     const gabac::EncodingConfiguration& configuration);
//
//    size_t writeBlock(const calq::EncodingOptions& opts, const gabac::EncodingConfiguration& configuration,
//                      const calq::DecodingBlock& block, const calq::SideInformation& side,
//                      const std::string& unmappedQualityValues_, bool STREAMOUT, size_t* compressedSizeMapped,
//                      size_t* compressedSizeUnmapped);

    size_t readHeader(size_t* blockSize);

    size_t writeHeader(const size_t& blockSize);

    size_t nrWrittenFileFormatBytes() const;

   private:
    static constexpr const char* MAGIC = "CQ";
    const size_t MAGIC_LEN = 3;

    size_t nrReadFileFormatBytes_;
    size_t nrWrittenFileFormatBytes_;

    size_t nrReadFileFormatBytes() const;

    size_t readQuantizers(std::vector<std::vector<uint8_t>>* quantizers);
//    size_t readQualBlock(std::string* block, const gabac::EncodingConfiguration& configuration);

    size_t writeQuantizers(const std::vector<std::vector<uint8_t>>& quantizers);
//    size_t writeQualBlock(unsigned char* block, const size_t& blockSize,
//                          const gabac::EncodingConfiguration& configuration);
};

}  // namespace cip

#endif  // CALQAPP_CQ_FILE_H_

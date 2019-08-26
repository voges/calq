#ifndef CIP_HELPERS_H_
#define CIP_HELPERS_H_

#include <string>

namespace cip {

bool fileExists(const std::string& path);

std::string fileNameExtension(const std::string& path);

}  // namespace cip

#endif  // CIP_HELPERS_H_

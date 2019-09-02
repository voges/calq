#ifndef UTIL_HELPERS_H_
#define UTIL_HELPERS_H_

#include <string>

namespace util {

bool fileExists(const std::string& path);

std::string fileNameExtension(const std::string& path);

}  // namespace util

#endif  // UTIL_HELPERS_H_

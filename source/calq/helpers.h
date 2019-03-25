#ifndef CALQ_HELPERS_H_
#define CALQ_HELPERS_H_

#include <string>

namespace calq {

std::string currentDateAndTime(void);
bool fileExists(const std::string &path);
std::string fileBaseName(const std::string &path);
std::string fileNameExtension(const std::string &path);
std::string removeFileNameExtension(const std::string &path);

}  // namespace calq

#endif  // CALQ_HELPERS_H_

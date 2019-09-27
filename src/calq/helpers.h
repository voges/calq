/**
 * @file helpers.h
 */

#ifndef CALQ_HELPERS_H_
#define CALQ_HELPERS_H_

#include <string>

namespace calq {

bool fileExists(const std::string& path);

std::string fileNameExtension(const std::string& path);

}  // namespace calq

#endif  // CALQ_HELPERS_H_

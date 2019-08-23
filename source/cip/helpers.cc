#include "helpers.h"
#include <fstream>

namespace cip {

bool fileExists(const std::string& path) {
    std::ifstream ifs(path.c_str());
    return ifs.good();
}

std::string fileNameExtension(const std::string& path) {
    if (path.find_last_of('.') != std::string::npos) {
        return path.substr(path.find_last_of('.') + 1);
    }
    return "";
}

}  // namespace cip

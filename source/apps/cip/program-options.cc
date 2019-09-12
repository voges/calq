#include "program-options.h"

#include <cli11@13becad/CLI11.hpp>
#include <sstream>

namespace cip {

ProgramOptions::ProgramOptions(int argc, char *argv[]) : decompress(), force() { processCommandLine(argc, argv); }

void ProgramOptions::processCommandLine(int argc, char *argv[]) {
    CLI::App app("CALQ");

    app.add_flag("-d,--decompress", decompress, "Decompress");
    app.add_flag("-f,--force", force, "Force overwriting of output file(s)");

    try {
        app.parse(argc, argv);
        validate();
    } catch (const CLI::ParseError &e) {
        if (app.exit(e) != 0) {
            throw std::runtime_error{"command line parsing failed: " + std::to_string(app.exit(e))};
        }
    }
}

void ProgramOptions::validate() {
    std::cout << "cip: validating command line" << std::endl;
    validateCommon();
}

void ProgramOptions::validateCommon() {
    // decompress
    if (decompress) {
        std::cout << "cip: decompressing" << std::endl;
    } else {
        std::cout << "cip: compressing" << std::endl;
    }

    // force
    if (force) {
        std::cout << "cip: force switch set - overwriting output file(s)" << std::endl;
    }
}

}  // namespace cip
